#!/usr/bin/perl

use strict;
use warnings;

package SEEDHMM;

#nhmmer location
use constant NHMMER_DIR => "/nfs/users/nfs_s/sb30/bin/hmmer3.1_alpha_0.30/src/";

#returns true if the SEED exists in the users cwd else dies (cleanly)
sub checkSEEDExists{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#checks SEED file exists
	if (!-e "$dirVars[1]/SEED"){
		FileHandling::cleanDie(\@dirVars, "SEED file does not exist in current directory.\n");
	}

	return 1;
}

#converts the original SEED to a stk file into the temporary working dir
sub convertSEEDToStk{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	
	#converts SEED to stockholm format
	if (system("sreformat -d stockholm $dirVars[1]/SEED > $dirVars[0]/SEED0.stk") != 0){
		FileHandling::cleanDie(\@dirVars, "Error converting SEED to Stockholm format. SEED may not be correctly formatted\n");
	}
}

#nhmmer hmmbuild
sub buildHMM{
	my $iter = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#Builds HMM
	if (system(NHMMER_DIR . "hmmbuild $dirVars[0]/SEED$iter.hmm $dirVars[0]/SEED$iter.stk > $dirVars[0]/hmmbuild$iter") != 0){
		FileHandling::cleanDie(\@dirVars, "hmmbuild failed.\n");
	}
}

#nhmmer hmmalign
sub alignSEED{
	my $iter = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;	

	my $seedIter = $iter+1;

	#Aligns HMM
	if (system(NHMMER_DIR . "hmmalign $dirVars[0]/SEED$iter.hmm $dirVars[0]/SEED$seedIter > $dirVars[0]/SEED$seedIter.stk") != 0){
		FileHandling::cleanDie(\@dirVars, "hmmalign failed.\n");
	}

	#removes repetitive sequences
	if (system("esl-weight -f --idf 0.99 $dirVars[0]/SEED$seedIter.stk > $dirVars[0]/SEED$seedIter-99") != 0){
		FileHandling::cleanDie(\@dirVars, "Removing repetitive sequences failed.\n");
	}
}

#for blast, aligns the file using mafft and then converts output SEED to stk and then remove repetitive seqs
sub blastAlign{
	my $iter = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my @seqStat = CreateMiniDB::getResidues("$dirVars[0]/SEED$iter", \@dirVars);

	#aligns the blast output
	if ($seqStat[0] > 1){
		if (system("mafft --quiet $dirVars[0]/SEED$iter > $dirVars[0]/SEEDm$iter") != 0){
			FileHandling::cleanDie(\@dirVars, "Aligning SEED failed.\n");
		}
		if (system("sreformat -d stockholm $dirVars[0]/SEEDm$iter > $dirVars[0]/SEED$iter.stk") != 0){
			FileHandling::cleanDie(\@dirVars, "Converting SEED to stockholm format failed.\n");
		}
	} else {
		if (system("sreformat -d stockholm $dirVars[0]/SEED$iter > $dirVars[0]/SEED$iter.stk") != 0){
			FileHandling::cleanDie(\@dirVars, "Converting SEED to stockholm format failed.\n");
		}		
	}

	#removes repetitive sequences
	if (system("esl-weight -f --idf 0.99 $dirVars[0]/SEED$iter.stk > $dirVars[0]/SEED$iter-99") != 0){
		FileHandling::cleanDie(\@dirVars, "Removing repetitive sequences failed.\n");
	}
	
}

#runs nhmmer
sub nhmmerSEED{
	my $iter = shift;
	my $miniDBCount = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my $miniDBPath;
	if ($miniDBCount < 0){
		$miniDBPath = "$dirVars[8]/bigMiniDBs";
		$miniDBCount *= -1;
	} else {
		$miniDBPath = $dirVars[0];
	}

	#generate table output for each minidb
	for (my $i = 1; $i <= $miniDBCount; $i++){
		my $miniDB = "$miniDBPath/miniDB$i.fa";
		#checks minidb exists
		if (!-e $miniDB){
			FileHandling::cleanDie(\@dirVars, "Can't find MiniDB.\n");
		} elsif (system(NHMMER_DIR . "nhmmer --tblout $dirVars[0]/SEED$iter-$i.nhmmer.tblout $dirVars[0]/SEED$iter.hmm $miniDB > $dirVars[0]/nhmmer$iter-$i") != 0) {
			FileHandling::cleanDie(\@dirVars, "nhmmer error.\n");
		}
		FileHandling::copyOut("$dirVars[0]/SEED$iter-$i.nhmmer.tblout", "$dirVars[1]/matthmmer$$/nhmmertblout");
	}
}

#blast
sub blastSEED{
	my $iter = shift;
	my $miniDBCount = shift;
	my $miniDBStatsRef = shift;
	my @miniDBStats = @$miniDBStatsRef;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#converts SEED to fa
	if ($iter > 0){
		if (system("sreformat -d fasta $dirVars[0]/SEED$iter-99 > $dirVars[0]/SEED$iter.fa") != 0){
			FileHandling::cleanDie(\@dirVars, "Converting SEED to fasta format failed.\n");
		}
	} else {
		#might be able to do this more efficiently...
		FileHandling::copyOut("$dirVars[0]/SEED$iter", "$dirVars[0]/SEED$iter.fa");		
	}

	my $miniDBPath;
	if ($miniDBCount < 0){
		$miniDBPath = "$dirVars[8]/bigMiniDBs";
		$miniDBCount *= -1;
	} else {
		$miniDBPath = $dirVars[0];
	}

	#generate table output for each minidb
	for (my $i = 1; $i <= $miniDBCount; $i++){
		my $miniDB = "$miniDBPath/miniDB$i.fa";
		#checks minidb exists
		if (!-e $miniDB){
			FileHandling::cleanDie(\@dirVars, "Can't find MiniDB.\n");
		} elsif (system("wublastn $miniDB $dirVars[0]/SEED$iter.fa -mformat 3 --e 0.01 -Z " . $miniDBStats[$i-1] . " -hspsepsmax 200 > $dirVars[0]/SEED$iter-$i.blast.out") != 0) {
			FileHandling::cleanDie(\@dirVars, "blast error.\n");
		}
		FileHandling::copyOut("$dirVars[0]/SEED$iter-$i.blast.out", "$dirVars[1]/matthmmer$$/blastout");
	}
}

#returns what type of SEED has been provided based on the SEED header
sub checkForValidSEEDHeader{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#gets the header of the SEED
	my $header = `head $dirVars[1]/SEED -n 1` or FileHandling::cleanDie(\@dirVars, "Could not take the first line of the SEED\n");
	chomp $header;

	if ($header =~ m/^>[A-Z0-9]+\.[0-9]\/[0-9]+-[0-9]+[\s]*$/){
		#format fasta header: ">acc.v/start-end"
		return 1;
	} elsif ($header =~ m/^# STOCKHOLM.*$/){
		#format stockholm
		return 2;
	} elsif ($header =~ m/^>.*$/){
		#format fasta ">junk"
		return 0;
	} else {
		#format not recognised so die
		FileHandling::cleanDie(\@dirVars, "SEED header is not valid.\n");
	}
}

#counts number of seqs in a fasta formatted SEED
sub getNoofSeqsFromFaSEED{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $index = shift;

	return `grep ">" $dirVars[0]/SEED$index | wc -l` or FileHandling::cleanDie(\@dirVars, "Cannot grep SEED$index.\n");
}

#checks which seqs are in the final SEED and which were repeated seqs (used to generate the seq_summary report)
sub checkWhichSeqsInSEED{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $finalIter = shift;

	foreach my $seq (keys %seqs){
		if(`grep $seq $dirVars[0]/SEED$finalIter-99` ne ""){
			$seqs{$seq}->{inSEED} = 1;
		}
	}

	return %seqs;	
}

#checks that the number of sequences has increased since last iteration
sub checkForConvergence{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $prevIter = shift;
	my $newIter = $prevIter+1;
	
	if (&getNoofSeqsFromFaSEED(\@dirVars, $newIter) eq &getNoofSeqsFromFaSEED(\@dirVars, $prevIter)){
		#converged
		return 1;
	} else {
		#new seqs added
		return 0;
	}
}

#creates the secondary structure for the final SEED
sub createSecondaryStructure{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $iter = shift;

	#converts final SEED in fasta format in secondary structure dir
	if (system("sreformat fasta $dirVars[1]/matthmmer$$/SEED$iter > $dirVars[1]/matthmmer$$/secondaryStructure/fastaSEED") != 0){
		FileHandling::cleanDie(\@dirVars, "Converting the final SEED to fasta format failed.\n");		
	}

	#secondary structure command
	if (system("runMAFFT_RNAalifold.sh $dirVars[1]/matthmmer$$/secondaryStructure/fastaSEED $dirVars[1]/matthmmer$$/secondaryStructure/ >& $dirVars[1]/matthmmer$$/secondaryStructure/output.txt") != 0){
		FileHandling::cleanDie(\@dirVars, "Creating the secondary structure failed.\n");		
	}	
}

1;
