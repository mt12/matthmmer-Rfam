#!/usr/bin/perl

use strict;
use warnings;
#use Getopt::Std;
use Getopt::Long;
use DateTime;
#use Data::Dumper;

push (@INC, "/nfs/users/nfs_m/mt12/scripts");

require Taxon;
require TaxonTree;
require FileHandling;
require CreateMiniDB;
require SEEDHMM;
require ExtractSeqs;

#Dir constants for location of files and dbs
use constant SCRATCH_DIR => "/warehouse/pfam01/rfam/Users/mt12";
use constant DBS_DIR => "/lustre/scratch101/blastdb/Rfam/rfamseq/EMBL_Bacterial_Genomes";

#All of the following constants are only used if the user does not specify them in the flags:
use constant DEFAULT_GROUP => 'family';					#The default group matthmmer traverses up to in creating a subset for the minidbs
use constant MAX_ITERATIONS_BEFORE_CONVERGENCE => 10;	#The default maximum number of iterations that matthmmer is limited to before stopping
use constant MIN_BIT_SCORE => 20;						#The default bit score threshold
use constant MIN_PERCENTAGE_LENGTH => 0.95;				#The default minimum length a nhmmer match can be as a percentage of the starting seq length
use constant MIN_EVALUE => 0.01;						#The default evalue threshold (if used)

my %options=();
#getopts("a:i:l:g:n:pbds:z:", \%options); #short opts
&GetOptions(	#cjkmruw (x)
	"a|accession=s" => \$options{a},	#The accession number for starting the sibling tree search on
	"t|taxonomy=i" => \$options{i},		#Alternatively, specifying a taxonomy id for starting the sibling tree search on
	"d|depth=i" => \$options{l},		#Depth to traverse up the tree
	"g|group=s" => \$options{g},		#Alternatively, specify a group to find up the tree
	"n|iterations=i" => \$options{n},	#Number of iterations to nhmmer/blast on
	"p|parent" => \$options{p},			#The accession/tax id given is a parent node, so find all the children of this
	"b|blast" => \$options{b},			#Use blast instead of nhmmer
	"s|secondary" => \$options{s},		#Create the secondary structure
	"dirty" => \$options{d},			#Keep temporary files that were on users scratch109
	"bitscore=i" => \$options{o},		#Set the bitscore threshold
	"e" => \$options{e},				#Use evalue instead of bitscore
	"evalue=f" => \$options{v},			#Set the evalue threshold
	"seqlength=f" => \$options{z},		#Set the sequence length percentage threshold
	"testHarness" => \$options{t},		#if testharness is being run
	"tempdir=s" => \$options{f},		#Changes the temporary dir to a user specified dir
	"q|quiet" => \$options{q},			#Supresses the output progress
	"v|verbose" => \$options{y},		#Verbosely runs the code with output on everything that is going on
	"h|help" => \$options{h},			#Shows the help page
);

#if help is requested, print the help page and exit
if(defined $options{h}) {
    FileHandling::help();
    exit(1);
}

#stupidity check...
if (defined $options{y} && defined $options{q}){
	die "So you want the script to run quietly AND verbosely??? I'm just going to die instead...\n";
}

#get times for the report.txt file
my $startTime = DateTime->now();
my $start = time;

#sets up the temp and output directories

my @dirVars = FileHandling::setupDirs($options{b}, $options{f}, $options{s});
#open connection to rfam_11_0 db:
push (@dirVars, FileHandling::openDBConnection(\@dirVars));
push (@dirVars, $options{d});
push (@dirVars, $options{q});
push (@dirVars, $options{y});
push (@dirVars, SCRATCH_DIR);
push (@dirVars, DBS_DIR);
#IMPORTANT - @dirVars is an array used throughout matthmmer it contains the following setup:
#$dirVars[0] = temporary working dir (probs users scratch109)
#$dirVars[1] = output dir in users cwd
#$dirVars[2] = db connection
#$dirVars[3] = db prepared statement
#$dirVars[4] = dirty
#$dirVars[5] = quiet
#$dirVars[6] = verbose
#$dirVars[7] = the scratch_dir (actually warehouse, the location of all files needed such as the nodes.dmp file)
#$dirVars[8] = the dbs_dir (the scratch dir where the dbs and bigMiniDBs are located)
#The majority of these variables are needed if at any point matthmmer needs to die, it needs to die gracefully, i.e. remove all files (unless dirty is specified, and close the db connection off.
if (!defined $options{q}){
	print "Directories and database connection to rfam_11_0 setup.\n";
}

my @meta;
my $seqLength;
my $alignment = 0;
#Checks a file called SEED exists in cwd
if (SEEDHMM::checkSEEDExists(\@dirVars)){
	if (defined $options{y}) { print "Checked that SEED exists\n"; }
	#gets info on the SEED based on the header and works out what type of SEED has been provided. Will then proceed to get meta on the SEED if valid
	my $isSEEDValid = SEEDHMM::checkForValidSEEDHeader(\@dirVars);
	if ($isSEEDValid == 1){
		if (defined $options{y}) { print "Checked that SEED header is in the format acc/start-end.\n"; }
		@meta = ExtractSeqs::getMetaFromSEED(1, 0, \@dirVars);
		if (defined $options{y}) { print "SEED accession is in database and SEED sequence matches that in the db.\n"; }
	} elsif ($isSEEDValid == 2){
		if (defined $options{y}) { print "SEED header is in stockholm format.\n"; }
		$alignment = 1;
		@meta = ExtractSeqs::getMetaFromSEED(1, 1, \@dirVars);
	} else {
		if (defined $options{y}) { print "SEED header is not in expected format, but could be valid, so continuing anyway.\n"; }
		$seqLength = ExtractSeqs::getMetaFromSEED(0, 0, \@dirVars);
	}
} else {
	#dies if SEED does not exist
	FileHandling::cleanDie(\@dirVars, "SEED file does not exist in your current working directory\n");
}

#if neither -a or -t flags are set, then set them -a to the default of 'SEED'
if (!defined $options{a} && !defined $options{i}){
	$options{a} = 'SEED';
}

#If an accession no is provided instead of a taxonomy id need to convert this into a taxonomy id
if (defined $options{a}){
	if (defined $options{i}){
		FileHandling::cleanDie(\@dirVars, "Can only have one of flag -a (accession no) or flag -t (taxonomy id) defined, not both.\n");		
	}
	#if the -a flag is SEED then sets the accession no to be the accession in the SEED
	if ($options{a} eq 'SEED'){
		if ($alignment){
			FileHandling::cleanDie(\@dirVars, "Cannot give -a SEED when the SEED file is an alignment file.\n");
		}
		if (!@meta){
			FileHandling::cleanDie(\@dirVars, "The SEED file provided doesn't have a real accession. This is ok, but you can't provide SEED as the -a flag variable.\n");					
		}
		#Assigns the accession number found in the SEED to the accession no option.
		$options{a} = $meta[0];
		if (defined $options{y}) { print "\'SEED\' has been specified as the accession the flag. The accession in the SEED is $options{a}.\n"; }
	}
	#removes the accession version if this is given with the accession
	$options{a} = Taxon::stripAccessionVersion($options{a});
	#converts the accession into equivalent taxonomy ids
	my %accessionHash = Taxon::getTaxonIDsByAccessionNo([($options{a})], \@dirVars);
	#takes the first tax id returned for the accession (even if there is multiple tax ids)
	if (defined $accessionHash{$options{a}}[0]){
		$options{i} = $accessionHash{$options{a}}[0];
		if (defined $options{y}) {
			if (scalar(@{$accessionHash{$options{a}}}) > 1){
				print "There is more than one equivalent taxonomy id for the given accession. Only the one shown below is being used.\n";
			}
			print "Accession provided is valid and the equivalent taxonomy id is: $options{i}.\n";
		}
	} else {
		FileHandling::cleanDie(\@dirVars, "$options{a} is not a valid accession number.\n");
	}	
}

#creates the report
if (defined $options{e}){
	FileHandling::createReport($startTime, \@dirVars, \%options, DEFAULT_GROUP, MAX_ITERATIONS_BEFORE_CONVERGENCE, MIN_EVALUE, MIN_PERCENTAGE_LENGTH);
} else {
	FileHandling::createReport($startTime, \@dirVars, \%options, DEFAULT_GROUP, MAX_ITERATIONS_BEFORE_CONVERGENCE, MIN_BIT_SCORE, MIN_PERCENTAGE_LENGTH);
}

#checks that taxon id given is valid
if (Taxon::taxonIDValid($options{i}, $options{p}, \@dirVars)){
	if (!defined $options{q}){
		print "Valid taxonomy id/accession number given and SEED exists.\n";
	}
	#gets the siblings for the given taxon id to the required depth/group
	my @siblings;
	if (defined $options{l}){
		@siblings = TaxonTree::getSiblings($options{i}, 0, $options{l}, 0, \@dirVars);
	} elsif (defined $options{p}){
		@siblings = TaxonTree::getSiblings($options{i}, 1, 0, 0, \@dirVars);
	} elsif (defined $options{g}){
		@siblings = TaxonTree::getSiblings($options{i}, 0, 0, $options{g}, \@dirVars);
	} else {
		#default
		@siblings = TaxonTree::getSiblings($options{i}, 0, 0, DEFAULT_GROUP, \@dirVars);
	}

	my $miniDBCount;
	my @miniDBStats;

	#If an array consisting of 1 element being -1 then this signifies that all the siblings would have been returned, so rather than create all the minidbs, use the bigMiniDBs!!
	if (scalar(@siblings) == 1 && $siblings[0] == -1){
		FileHandling::addNoofSiblings(\@dirVars, "All");
		if (!defined $options{q}){
			print "The whole taxonomy tree is being used, minidbs have been pre-created for this.\n";
		}
		$miniDBCount = FileHandling::countNoOfBigMiniDBs(\@dirVars);
		@miniDBStats = CreateMiniDB::getBigMiniDBStats(\@dirVars, $miniDBCount);
		#inverting the minidbcount saves a flag for future to inform scripts that bigMiniDBs are being used
		$miniDBCount *= -1;
	} else {
		#not all the siblings returned so need to create the minidbs!!
		if (!defined $options{q}){
			print scalar(@siblings) . " siblings have been found.\n";
		}
		FileHandling::addNoofSiblings(\@dirVars, scalar(@siblings));
		my @accessions;

		#converts each siblings taxonomy id to the accession no(s) (with version no).
		my %accessionsHash = Taxon::getAccessionNosByTaxonID(\@siblings, 1, \@dirVars);	
		foreach my $accessionNos (keys %accessionsHash){
			push (@accessions, @{$accessionsHash{$accessionNos}});
		}

		#finds which database each accession no is in and returns this as a hash
		my %dbs = CreateMiniDB::findDatabases(\@accessions, 0, \@dirVars);

		#creates the minidb(s) from the accession nos
		@miniDBStats = CreateMiniDB::fetchSeqs(\%dbs, \@dirVars);
		$miniDBCount = scalar(@miniDBStats);
		if (!defined $options{q}){
			print "$miniDBCount MiniDB(s) have been created.\n";
		}
	}

	FileHandling::addSEEDInfo(\@dirVars, \@meta, $options{b});

	#set params to default if not specified:
	if (!defined $options{n}){
		$options{n} = MAX_ITERATIONS_BEFORE_CONVERGENCE;
		if (defined $options{y}) { print "Script will iterate until convergence or until the number of iterations reaches the maximum iteration limit of: " . MAX_ITERATIONS_BEFORE_CONVERGENCE . ".\n"; }
	}
	if (!defined $options{z}){
		$options{z} = MIN_PERCENTAGE_LENGTH;
		if (defined $options{y}) { print "The minimum length percentage for a sequence to be accepted has been set to the default value of: " . MIN_PERCENTAGE_LENGTH . ".\n"; }
	}
	if (defined $options{e}){
		if (defined $options{v}){
			$options{o} = $options{v};
		} else {
			$options{o} = MIN_EVALUE;
			if (defined $options{y}) { print "The evalue threshold has been set to the default value of: " . MIN_EVALUE . ".\n"; }
		}
	} else {
		if (!defined $options{o}){
			$options{o} = MIN_BIT_SCORE;
			if (defined $options{y}) { print "The bit score threshold has been set to the default value of: " . MIN_BIT_SCORE . ".\n"; }
		}
	}

	FileHandling::copyOrigSEED(\@dirVars);

	my %seqs;

	if (defined $options{b}){ #use blast

		if ($miniDBCount > 0){
			CreateMiniDB::blastXDFormat($miniDBCount, \@dirVars);
		}

		for (my $i = 0; $i < $options{n}; $i++){
			#Runs blast on the SEED and the minidb(s)
			SEEDHMM::blastSEED($i, $miniDBCount, \@miniDBStats, \@dirVars);

			if (@meta){
				#Finds the closest matching sequences from the nhmmer output
				%seqs = ExtractSeqs::getGoodBlastSeqs($i, 1, \@meta, $miniDBCount, $options{o}, $options{z}, $options{e}, \%seqs, $options{n}, \@dirVars);
			} else {
				%seqs = ExtractSeqs::getGoodBlastSeqs($i, 0, $seqLength, $miniDBCount, $options{o}, $options{z}, $options{e}, \%seqs, $options{n}, \@dirVars);
			}

			#Takes these seqs and puts them in a new SEED file
			ExtractSeqs::blastFetchSeqs($i+1, \%seqs, \@dirVars);

			SEEDHMM::blastAlign($i+1, \@dirVars);

			#Copies the newly created SEED from the temp rough working dir on users scratch to the output dir in users cwd
			FileHandling::copySEEDToOutput(\@dirVars, $i+1);

			if (!defined $options{q}){
				print "Using BLAST, Iteration $i has been completed successfully.\n";
			}

			#check to see that the number of sequences is greater than last time, if not then finish blasting
			if (SEEDHMM::checkForConvergence(\@dirVars, $i)){
				$options{n} = $i+1;
				last;
			}
		}

	} else { #use nhmmer (default)

		if (defined $alignment){
			#if the SEED file is in stockholm format, then makes a copy of the SEED in fasta format for use later on
			FileHandling::formatAlignmentSEED(\@dirVars);
		} else {
			#converts SEED to stockholm format to allow the HMM to be built
			SEEDHMM::convertSEEDToStk(\@dirVars);
		}

		for (my $i = 0; $i < $options{n}; $i++){
			#Build HMM from SEED file
			SEEDHMM::buildHMM($i, \@dirVars);

			#Runs nhmmer on the HMM and the minidb(s)
			SEEDHMM::nhmmerSEED($i, $miniDBCount, \@dirVars);

			if (@meta){
				#Finds the closest matching sequences from the nhmmer output
				%seqs = ExtractSeqs::getGoodNhmmerSeqs($i, 1, \@meta, $miniDBCount, $options{o}, $options{z}, $options{e}, \%seqs, $options{n}, \@dirVars);
			} else {
				%seqs = ExtractSeqs::getGoodNhmmerSeqs($i, 0, $seqLength, $miniDBCount, $options{o}, $options{z}, $options{e}, \%seqs, $options{n}, \@dirVars);
			}

			#Takes these seqs and puts them in a new SEED file
			ExtractSeqs::nhmmerFetchSeqs($i+1, \%seqs, \@dirVars);

			#Aligns the new SEED with the old HMM and removes redundant seqs
			SEEDHMM::alignSEED($i, \@dirVars);

			#Copies the newly created SEED from the temp rough working dir on users scratch to the output dir in users cwd
			FileHandling::copySEEDToOutput(\@dirVars, $i+1);

			if (!defined $options{q}){
				print "Using NHMMER, Iteration $i has been completed successfully.\n";
			}

			#check to see that the number of sequences is greater than last time, if not then finish nhmmering
			if (SEEDHMM::checkForConvergence(\@dirVars, $i)){
				$options{n} = $i+1;
				last;
			}
		}

	}

	#creates the secondary structure, if requested
	if (defined $options{s}){
		SEEDHMM::createSecondaryStructure(\@dirVars, $options{n});
		if (!defined $options{q}){
			print "Secondary Structure Created.\n";
		}
	}

	#creates the species-list and genus-list reports
	TaxonTree::getAllSeqsParents(\@dirVars, \%seqs, $options{n}, 'species');
	TaxonTree::getAllSeqsParents(\@dirVars, \%seqs, $options{n}, 'genus');

	#writes all the data to seq_summary.txt
	%seqs = SEEDHMM::checkWhichSeqsInSEED(\@dirVars, \%seqs, $options{n});
	FileHandling::createSeqSum(\@dirVars, $options{n}, $options{e});
	FileHandling::addSearchInfo(\@dirVars, $options{n}, $options{b}, $options{e}, \%seqs);

	#adds SEED stats to the report
	FileHandling::addFinalSEEDInfo(\@dirVars, $options{n});

	#if dirty has not been specified then remove the temp dir and files
	if (!$options{d}){
		FileHandling::removeTempDir($dirVars[0], "");
		if (!defined $options{q}){
			print "Clean up completed.\n";
		}
	}

	#finishes making the report.txt file
	FileHandling::finishReport(\@dirVars, time - $start);

	#closes the database connection
	FileHandling::closeDBConnection(\@dirVars);

	#leaves a file behind with the process id for the test harness to read, if the test harness is being run
	if ($options{t}){
		FileHandling::leaveTestHarnessID(\@dirVars);
	}

} else {
	FileHandling::cleanDie(\@dirVars, "$options{i} is not a valid taxonomy ID.\n");
}

if (!defined $options{q}){
	print "Well that all worked! Check you current directory for a new directory called matthmmer$$ for your files.\n";
}

1;
