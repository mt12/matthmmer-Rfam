#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use File::Copy;
use File::Path;
use DateTime;

use DBI;
use Rfam;

package FileHandling;

#creates temporary working dir, and cwd output dir.
sub setupDirs{
	my $blast = shift;	
	my $altTempDir = shift;
	my $secondaryStructure = shift;

	#sets the temp dir if specified by the user, otherwise use scratch109
	my $tempdir;
	if ($altTempDir){
		$tempdir = "$altTempDir/$$";
	} else {
		#gets the user id
		my $user = getpwuid($<);
		$tempdir = "/lustre/scratch109/sanger/$user/$$";
	}

	#get users current working dir
	my $cwd = Cwd::getcwd();
	#creates a directory on users scratch109 to store all the temporary files created
	unless(mkdir $tempdir){
		die "Unable to create a temporary dir in your scratch109.\n";
	}
	#creates a directory in cwd to return the SEEDs
	unless(mkdir "$cwd/matthmmer$$"){
		die "Unable to create a dir in your cwd for the results.\n";
	}
	#creates a directory inside the one just made for the nhmmer.tblout/blast.out files.
	my $out = "nhmmertblout";
	if ($blast){
		$out = "blastout";
	}
	#creates a dir for the nhmmer/blast output
	unless(mkdir "$cwd/matthmmer$$/$out"){
		die "Unable to create a dir in your cwd for the $out results.\n";
	}

	#creates a dir for the secondary structure information if the user has specified secondary structure
	if($secondaryStructure){
		unless(mkdir "$cwd/matthmmer$$/secondaryStructure"){
			die "Unable to create a dir in your cwd for the secondary structure results.\n";
		}
	}

	return ($tempdir, $cwd);
}

#copies the original SEED into the temporary working dir
sub copyOrigSEED{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	File::Copy::copy("$dirVars[1]/SEED", "$dirVars[0]/SEED0");	
}

sub copySEED{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $iter = shift;
	my $prevIter = $iter-1;

	File::Copy::copy("$dirVars[0]/SEED$prevIter", "$dirVars[0]/SEED$iter");	
}

#copys the SEED from the temporary dir to the output dir in the users cwd
sub copySEEDToOutput{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $iter = shift;

	File::Copy::copy("$dirVars[0]/SEED$iter-99", "$dirVars[1]/matthmmer$$/SEED$iter");
}

#Copys a file from one dir to another (used for moving SEED files from the temp dir to the final results dir)
sub copyOut{
	my $path = shift;
	my $outDir = shift;
	
	File::Copy::copy($path, $outDir);
}

sub formatAlignmentSEED{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	&copyOut("$dirVars[1]/SEED", "$dirVars[0]/SEED0.stk");
	
	if(system("sreformat fasta $dirVars[0]/SEED0.stk > $dirVars[0]/SEED0") != 0){
		&cleanDie($dirVars[0], $dirVars[1], "Unable to create the report file.\n");
	}
}

#Deletes the temp dir (that was on users scratch)
sub removeTempDir{
	my $rmdir = shift;
	my $message = shift;
	File::Path::rmtree("$rmdir") or die "$message and can't clean up";
}

#sub that cleans up (removes the temp dirs) and then displays the message before dieing
sub cleanDie{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $message = shift;
	#close the db connection
	&closeDBConnection(\@dirVars);
	#removes the tempdir in the users scratch and the output dir in the cwd
	if ($dirVars[4]){
		print "Kept temporary files as dirty was specified in $$.\n";
	} else {
		&removeTempDir($dirVars[0], $message);
		&removeTempDir("$dirVars[1]/matthmmer$$", $message);
	}
	die $message;
}

#***REPORT TEXT CODE***#

sub createReport{
	my $startTime = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $optionsRef = shift;
	my %options = %$optionsRef;
	my $defaultGroup = shift;
	my $defaultNoofIters = shift;
	my $defaultMinBitScore = shift;
	my $defaultMinPercentageLength = shift;

	my $user = getpwuid($<);

	#creates the report.txt file
	open (REPORT, ">$dirVars[1]/matthmmer$$/report.txt") or &cleanDie($dirVars[0], $dirVars[1], "Unable to create the report file.\n");
		
		print REPORT	 "#matthmmer Report Automatically Generated.\n"
						."#******************************************************************************#\n"
						."Produced On:                \t".$startTime->dmy('-') . " " . $startTime->strftime('%H:%M:%S') . "\n"
						."User:                       \t$user\n"
						."ProcessID:                  \t$$\n"
						."Database Directory:         \t$dirVars[8]\n"
						."Temporary Working Directory:\t$dirVars[0]\n"
						."\n#***PARAMETERS:****************************************************************#\n"
						."SEED:				 \n";

		open (SEED, "$dirVars[1]/SEED") or &cleanDie($dirVars[0], $dirVars[1], "Unable to open SEED.\n");
			
			while(<SEED>){
				print REPORT $_;
			}
		
		close(SEED);

		#prints the parameters
		print REPORT	 "\nAccession No:                     \t";
		if (defined $options{a}){
			print REPORT $options{a} . "	(Given)\n"
						."Taxonomy ID:                      \t$options{i}	(Converted)\n";
		} else {
			print REPORT "Not defined\n"
						."Taxonomy ID:                      \t$options{i}	(Given)\n";			
		}
		print REPORT	"Parent?:                          \t";
		if (defined $options{p}){
			print REPORT "Yes\n";
		} else {
			print REPORT "No\n";
		}
		print REPORT	"Depth:                            \t";
		if (defined $options{l}){
			print REPORT $options{l} . "\n";
		} else {
			print REPORT "Not defined\n";
		}		
		print REPORT	"Group:                            \t";
		if (!defined $options{p} && !defined $options{l} && !defined $options{g}){
			$options{g} = 'family (Default)';
		}
		if (defined $options{g}){
			print REPORT $options{g} . "\n";	
		} else {
			print REPORT "Not defined\n";
		}
		print REPORT	"Search:                           \t";
		if (defined $options{b}){
			print REPORT "BLAST\n";
		} else {
			print REPORT "NHMMER\n";
		}
		print REPORT	"Number of Iterations:             \t";
		if (defined $options{n}){
			print REPORT $options{n} . "\n";
		} else {
			print REPORT $defaultNoofIters . " (Default maximum before terminating)\n";
		}
		print REPORT	 "Check for convergence?:           \t";
		if (defined $options{c}){
			print REPORT "Yes\n";
		} else {
			print REPORT "No\n";
		}		
		print REPORT	 "Use e-values instead of bit Score:\t";
		if (defined $options{e}){
			print REPORT "Yes\n";
			print REPORT "Maximum e-value:                  \t";
		} else {
			print REPORT "No\n";
			print REPORT "Minimum bit score:                \t";
		}
			if (defined $options{s}){
				print REPORT $options{s} . "\n";
			} else {
				print REPORT $defaultMinBitScore . " (Default)\n";
			}
		print REPORT	"Minimum Percentage Length:        \t";
		if (defined $options{z}){
			print REPORT $options{z} . "\n";
		} else {
			print REPORT $defaultMinPercentageLength . " (Default)\n";
		}
		print REPORT	"Keep Temporary Files:             \t";
		if (defined $options{d}){
			print REPORT "Yes\n";
		} else {
			print REPORT "No\n";
		}
		print REPORT "\n#***TAXONOMY TREE INFO:********************************************************#\n";

	close(REPORT);	
}

#adds the tree data to the report
sub addTreeInfo{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $parentID = shift;
	my $parentGroup = shift;	

	my $taxonStr = Taxon::getTaxonStr($parentID, \@dirVars);

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		print REPORT "Parent Node: Taxonomy ID:    \t$parentID\n"
					."             Taxonomy String:\t$taxonStr\n"
					."             Group:          \t$parentGroup\n";
		
	close(REPORT);
}

#adds the total number of siblings to the report
sub addNoofSiblings{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $noofSiblings = shift;

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		print REPORT "Number of siblings found:    \t$noofSiblings\n"
					."\n#***MINIDB INFO:***************************************************************#\n";
		
	close(REPORT);		
}

#adds the seqstat info for the minidbs
sub addMiniDBStats{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $miniDBIndex = shift;
	my $finalMiniDB = shift;
	my $seqstatRef = shift;
	my @seqstat = @$seqstatRef;

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		print REPORT 	 "MINIDB$miniDBIndex: Number of sequences:\t$seqstat[0]\n"
						."         Total # residues:   \t$seqstat[1]\n";
		if ($finalMiniDB){
			print REPORT "Total number of minidbs:     \t$miniDBIndex\n";
		}					
		
	close(REPORT);	
}

#adds the info about the original SEED to the report
sub addSEEDInfo{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $metaRef = shift;
	my @meta = @$metaRef;
	my $blast = shift;

	#gets the taxon string from the accession number
	my $taxonStr = "";
	if (@meta){
		my $accession = Taxon::stripAccessionVersion($meta[0]);
		my %accessionHash = Taxon::getTaxonIDsByAccessionNo([($accession)], \@dirVars);
		$taxonStr = Taxon::getTaxonStr($accessionHash{$accession}[0], \@dirVars);
	}

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		if (@meta){
			print REPORT "\n#***ORIGINAL SEED INFO:********************************************************#\n"
						."Accession No:   \t$meta[0]\n"
						."Taxonomy String:\t$taxonStr\n"
						."Sequence Length:\t$meta[3]\n"
						."Start Co-ord:   \t$meta[1]\n"
						."End Co-ord:     \t$meta[2]\n";
		}
		if ($blast){
			print REPORT "\n#***BLAST RESULTS INFO:********************************************************#\n";
		} else {
			print REPORT "\n#***NHMMER RESULTS INFO:*******************************************************#\n";
		}
	
	close(REPORT);	
}

#creates the seq_summary.txt file
sub createSeqSum{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $noofIters = shift;
	my $useEValue = shift;

	open (SEQSUM, ">$dirVars[1]/matthmmer$$/seq_summary.txt") or &cleanDie(\@dirVars, "Unable to create the sequence summary file.\n");

		print SEQSUM	 "#matthmmer Sequence Summary File Automatically Generated.\n"
						."#******************************************************************************#\n";
		#				."                                                                           Co-Ords      ";
		#if ($useEValue){
		#	print SEQSUM "E-Values\n";
		#} else {
		#	print SEQSUM "Bit Scores\n";
		#}
		print SEQSUM	 "Accession \tSeqLabel Level                                     \tLength\t   Start\t     End";
		for (my $i = 1; $i <= $noofIters; $i++){
			print SEQSUM "\tIter$i";
		}
		print SEQSUM "\n";

	close (SEQSUM);
}

#Makes the sort conditions as a string (as the sort params are variable length) so it can be 'eval'ed later.
sub makeSortStr{
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $noofIters = shift;

	my $sortStr = " ";

	#sorts by each bitscore starting with the final bitscore and working up to the first iteration's bit score as the least significant sort
	for (my $s = $noofIters; $s >= 0; $s--){
		$sortStr .= "sort {\$b->{bitScore}[$s] <=> \$a->{bitScore}[$s]} ";
	}

	return $sortStr . "values %seqs";
}

#Adds the output of nhmmer/blast results info to the report
sub addSearchInfo{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $noofIters = shift;
	my $blast = shift;
	my $useEValue = shift;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;

	use constant MAX_ACCESSION_LENGTH => 9;
	use constant MAX_LEVEL_LENGTH => 45;
	use constant MAX_LENGTH_LENGTH => 4;
	use constant MAX_COORDS_LENGTH => 8;
	my $maxScoreLength;
	if ($useEValue){
		$maxScoreLength = 7;
	} else {
		$maxScoreLength = 5;
	}

	open (SEQSUM, ">>$dirVars[1]/matthmmer$$/seq_summary.txt") or &cleanDie(\@dirVars, "Unable to ammend the sequence summary file.\n");
		
		my $sortStr = &makeSortStr(\%seqs, $noofIters-1);
		foreach my $seq (eval $sortStr){
			#sort {$b->{bitScore}[$noofIters-1] <=> $a->{bitScore}[$noofIters-1]}
			#sort {$b->{bitScore}[1] <=> $a->{bitScore}[1]}
			#sort {$b->{bitScore}[0] <=> $a->{bitScore}[0]}
			print SEQSUM "$seq->{accession}";
			my $i;
			for ($i = 0; $i < (MAX_ACCESSION_LENGTH - length($seq->{accession})); $i++){
				print SEQSUM " ";
			}
			if ($seq->{inSEED}){
				print SEQSUM "\tSEED \t";
			} else {
				print SEQSUM "\tALIGN\t";
			}
			my $accession = Taxon::stripAccessionVersion($seq->{accession});
			my %accessionHash = Taxon::getTaxonIDsByAccessionNo([($accession)], \@dirVars);
			my $taxonStr = Taxon::getTaxonStr($accessionHash{$accession}[0], \@dirVars);
			if (!$taxonStr){
				$taxonStr = " ";
			}
			print SEQSUM "$taxonStr";
			if (length($taxonStr) > (MAX_LEVEL_LENGTH + (MAX_LENGTH_LENGTH - length($seq->{realLength})))){
				print SEQSUM "\n                                                                ";
				for ($i = 0; $i < (MAX_LENGTH_LENGTH - length($seq->{realLength})); $i++){
					print SEQSUM " ";
				}
			} else {
				for ($i = 0; $i < ((MAX_LENGTH_LENGTH+MAX_LEVEL_LENGTH) - (length($taxonStr)+length($seq->{realLength}))); $i++){
					print SEQSUM " ";
				}
			}		
			print SEQSUM "\t$seq->{realLength}\t";
			for ($i = 0; $i < (MAX_COORDS_LENGTH - length($seq->{start})); $i++){
				print SEQSUM " ";
			}	
			print SEQSUM "$seq->{start}\t";
			for ($i = 0; $i < (MAX_COORDS_LENGTH - length($seq->{end})); $i++){
				print SEQSUM " ";
			}	
			print SEQSUM "$seq->{end}\t";
			for ($i = 0; $i < $noofIters; $i++){
				for (my $j = 0; $j < ($maxScoreLength - length($seq->{bitScore}[$i])); $j++){
					print SEQSUM " ";
				}				
				if ($seq->{bitScore}[$i] > 0){
					print SEQSUM "$seq->{bitScore}[$i]\t";
				} else {
					print SEQSUM "\t";
				}
			}
			print SEQSUM "\n";
		}	
		
	close(SEQSUM);	
}

#adds info on the Final SEED to the report
sub addFinalSEEDInfo{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $noofSEEDs = shift;

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		for (my $i = 1; $i <= $noofSEEDs; $i++){
			my $noofSeqs = SEEDHMM::getNoofSeqsFromFaSEED(\@dirVars, $i);
			print REPORT "Iteration $i - Number of sequences:\t$noofSeqs";
		}

		print REPORT 	"\n#***SEED RESULTS INFO:*********************************************************#\n";
	
		my @avSeqLength = ();
		for (my $i = 1; $i <= $noofSEEDs; $i++){
			my @seqstat = CreateMiniDB::getResidues("$dirVars[1]/matthmmer$$/SEED$i", \@dirVars);
			push (@avSeqLength, $seqstat[2]);
			print REPORT "SEED $i - Number of sequences:\t\t$seqstat[0]\n";
		}
		print REPORT	 "\n";
		for (my $i = 1; $i <= $noofSEEDs; $i++){			
			print REPORT "SEED $i - Average sequence length:\t" . $avSeqLength[$i-1] . "\n";
		}

	close(REPORT);		
}

#finishes the report
sub finishReport{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $endTime = shift;

	open (REPORT, ">>$dirVars[1]/matthmmer$$/report.txt") or &cleanDie(\@dirVars, "Unable to ammend the report file.\n");

		print REPORT "\n#***FINAL REPORT STATISTICS:***************************************************#\n"
					."Total time to run: $endTime seconds\n"
					."#                                                                              #\n"
					."#******************************************************************************#\n";
		
	close(REPORT);		
}

sub createLevelList{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $level = shift;
	my $noofIters = shift;
	my $levelsRef = shift;
	my %levels = %$levelsRef;

	use constant MAX_MATCH_LENGTH => 5;

	open (LIST, ">$dirVars[1]/matthmmer$$/$level-list.txt") or &cleanDie(\@dirVars, "Unable to create the $level-list file.\n");
		
		if ($level eq "genus"){
			print LIST "#***LIST OF RESULTS BY GENUS FAMILY:*******************************************#\n\n";
		} else {
			print LIST "#***LIST OF RESULTS BY SPECIES FAMILY:*****************************************#\n\n";
		}
		for (my $i = 1; $i < $noofIters+1; $i++){
			print LIST "Iter$i ";
		}
		print LIST " Family\n";

		my $sortStr = " ";

		#sorts by each bitscore starting with the final bitscore and working up to the first iteration's bit score as the least significant sort
		for (my $s = $noofIters-1; $s >= 0; $s--){
			$sortStr .= "sort {\$levels{\$b}[$s] <=> \$levels{\$a}[$s]} ";
		}

		$sortStr .= "keys %levels";

		foreach my $level (eval $sortStr){
			my $taxonStr = Taxon::getTaxonStr($level, \@dirVars);
			for (my $i = 0; $i < $noofIters; $i++){
				for (my $j = 0; $j < (MAX_MATCH_LENGTH - length($levels{$level}[$i])); $j++){
					print LIST " ";
				}
				print LIST $levels{$level}[$i] . " ";
			}
			print LIST " " . $taxonStr . "\n";
		}		

	close (LIST);
}

sub countNoOfBigMiniDBs{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my $count = 1;
	while (-e "$dirVars[8]/bigMiniDBs/miniDB$count.fa"){
		$count++;
	}
	$count--;
	if ($count > 0){
		return $count;
	} else {
		&cleanDie(\@dirVars, "No bigMiniDBs found...\n");
	}
}

#opens the db connection to rfam_11_0, and returns the db params
sub openDBConnection{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#RDB on dev connection stuff
	my $dbHost=$Rfam::rdbHostDev;
	my $dbUser=$Rfam::rdbUserDev;
	my $dbPort=$Rfam::rdbPortDev;
	my $dbPass=$Rfam::rdbPassDev;
	my $dbName=$Rfam::rdbNameDev;

	my $dsn    = "dbi:mysql:$dbName:$dbHost:$dbPort";
	my $dbAttr = { RaiseError => 1, PrintError => 1 };

	#connect
	my $dbh = DBI->connect( $dsn, $dbUser, $dbPass, $dbAttr )
	  or &cleanDie(\@dirVars, "(EE) ERROR: couldn't connect to database: $! \n");

	#prepare the statement
	my $asth = $dbh->prepare( 'SELECT level from taxonomy_websearch where ncbi_id=?' )
	  or &cleanDie(\@dirVars, "(EE) ERROR: couldn't prepare query to retrieve taxonomy string: " . $dbh->errstr . "\n");

	#return the connection and the prepared statment in an array
	return ($dbh, $asth);
}

#finishes the db query and closes the connection to the db
sub closeDBConnection{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	$dirVars[3]->finish();

	$dirVars[2]->disconnect();

	if (!$dirVars[5]){
		print "Rfam Database connection closed.\n";
	}
}

sub leaveTestHarnessID{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	
	open (FILE, ">$dirVars[1]/testHarnessID.txt") or &cleanDie(\@dirVars, "Cannot create file to store process id.\n");

		print FILE "$$";

	close (FILE);
}

#prints the help page
sub help{
    print STDERR <<EOF;

matthmmer.pl:   matthmmer traverses the taxonomy tree from a given started point and a given range building a subset from this
                using a given SEED file in the cwd matthmmer iteratively runs nhmmer on the SEED and HMM created from the subset
                creating a larger SEED that grows until all relevant sequences are found.

matthmmer.pl <options>
GENERAL OPTIONS:        -h|--help       Shows this help page
TAXONOMY TREE OPTIONS:  -a|--accession  Provide the accession number for starting the sibling tree search on
                        -t|--taxonomy   Alternatively, specify a taxonomy id for starting the sibling tree search on
                        -d|--depth      Depth to traverse up the tree
                        -g|--group      Alternatively, specify a group to traverse up the tree until
                        -p|--parent     The accession/tax id given is a parent node, so find all the children of this
NHMMER/BLAST OPTIONS:   -b|--blast      Use blast instead of nhmmer
                        -n|--iterations Number of iterations to nhmmer/blast on
                        --bitscore      Set the bitscore threshold
                        -e              Use evalue instead of bitscore
                        --evalue        Set the evalue threshold
                        --seqlength     Set the sequence length percentage threshold
FILE OPTIONS:           -s|--secondary  Produces the secondary structure of the final output
                        --tempdir       Changes the temporary dir to a user specified dir
                        --dirty         Keep temporary files that were in the users temporary dir
OTHER OPTIONS:          --testHarness   Set if running a test harness
                        -q|--quiet      Supresses the output progress

For more information on running matthmmer see also the README file located at:
/warehouse/pfam01/rfam/Users/mt12/README

EOF
}

1;
