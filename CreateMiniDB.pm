#!/usr/bin/perl

use strict;
use warnings;
#use Data::Dumper;

package CreateMiniDB;

#returns a hash of databases with the accession nos in them. Takes an array of accession numbers as the only arg
sub findDatabases{
    my $accessionsRef = shift;
    my @accessions = @$accessionsRef;
	my $singleAcc = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

    my %databases;

	#opens the file containing the list of dbs with the corresponding first and last seq
	if ($dirVars[0]){
    	open (FILE, "$dirVars[8]/seqIndex/seqIndex.txt") or FileHandling::cleanDie(\@dirVars, "Can't find/open seqIndex file for reading.\n");
	} else {
    	open (FILE, "$dirVars[8]/seqIndex/seqIndex.txt") or die "Can't find/open seqIndex file for reading.\n";
	}

        while(<FILE>){
            chomp $_;
            my @line = split(/\t/, $_);
			#loops through each accession for each line in the seqIndex
            foreach my $accession (@accessions){
				#if the accession no is >= the first accession in this db and <= the last accession in this db
				#then this accession no is appended to the hash of databases for the db in question
                if ($accession ge $line[1] && $accession le $line[2]){
					if ($singleAcc){
						return $line[0];
					} else {
		                if (exists $databases{$line[0]}) {
		                    push (@{$databases{$line[0]}}, $accession);
		                } else {
		                    $databases{$line[0]}[0] = $accession;
		                }
					}
                }
            }
        }

    close (FILE);
	
	#If only one accession was provided it should have returned in the loop, so if it has not then there has been a problem
	if ($singleAcc && $dirVars[0]){
		FileHandling::cleanDie(\@dirVars, "$accessions[0] not found in seqIndex.\n");
	} else {
    	return %databases;
	}
}

#returns stats on the bigMiniDBs!
sub getBigMiniDBStats{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $noofBigMiniDBs = shift;

	my @miniDBStats = ();

	#loops through each miniDB and finds stats on each, pushing to an array
	for (my $i = 1; $i <= $noofBigMiniDBs; $i++){
		my @dbStats = &getResidues("$dirVars[8]/bigMiniDBs/miniDB$i.fa", \@dirVars);
		push (@miniDBStats, $dbStats[1]);
	}

	return @miniDBStats;
}

#takes a fasta file and returns the number of sequences and total number of residues from it by using seqstat
sub getResidues{
    my $file = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

    my @seqstat = `seqstat --quiet $file` or FileHandling::cleanDie(\@dirVars, "Can't find $file to perform seqstat on\n");
	my $noofSeqs = $seqstat[2];
    my $residues = $seqstat[3];
	my $largest  = $seqstat[5];
	my $avlength = $seqstat[6];

    $noofSeqs =~ s/.*:\s*//;
    chomp $noofSeqs;
    $residues =~ s/.*:\s*//;
    chomp $residues;
    $largest =~ s/.*:\s*//;
    chomp $largest;
    $avlength =~ s/.*:\s*//;
    chomp $avlength;
	my @seqstatRet = ($noofSeqs, $residues, $avlength, $largest);

	#returned as an array
    return @seqstatRet;
}

#makes minidb files containing all of the sequences required. (Takes the returned hash from &findDatabases as the arg
sub fetchSeqs{
	my $databasesRef = shift;
	my %databases = %$databasesRef;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#biggest size of a minidb
	use constant MAX_FILE_SIZE => 3500000000;
	#biggest number of sequences xdget can return in one go
	use constant MAX_XDGET_SIZE => 3000;

	#Closes the MINIDB file, and checks that the number of sequences and total number of residues match up with what is expected
	sub closeAndCheckMiniDB{
		my $miniDBIndex = shift;
		my $miniDBSize = shift;
		my $seqCount = shift;
		my $finalMiniDB = shift;
		my $notUpdate = shift;
		my $dirVarsRef = shift;
		my @dirVars = @$dirVarsRef;

		close(MINIDB);
		my @seqstat = &getResidues("$dirVars[0]/miniDB$miniDBIndex.fa", \@dirVars);
		#if either are different throw an error
		if ($seqstat[0] != $seqCount || $seqstat[1] != $miniDBSize){
			warn "miniDB$miniDBIndex.fa has been created but doesn't pass seqstat validation checks. You may want to rerun this script.\n";
		}
		if ($notUpdate){
			FileHandling::addMiniDBStats(\@dirVars, $miniDBIndex, $finalMiniDB, \@seqstat);	
		}

		return $seqstat[1];
	}

	my @miniDBStats = ();
	my $miniDBSize = 0;
	my $seqCount = 0;
	#opens the first MINIDB to write to
	open(MINIDB, ">$dirVars[0]/miniDB1.fa") or FileHandling::cleanDie(\@dirVars, "Could not create miniDB1 to write to\n");

	#loops through each db in the hash
	foreach my $db (keys %databases){
		while(@{$databases{$db}}){
			#creates a long string of the accession nos, as required for xdget
			my $accessionList = join(" ", splice(@{$databases{$db}}, 0, MAX_XDGET_SIZE));
			#opens xdget on the db, with the accession list
		    open( XD, "xdget -n $dirVars[8]/$db $accessionList |") or FileHandling::cleanDie(\@dirVars, "Could not open xdget pipe:[$!]\n");
				while(<XD>){
					next if not defined($_);
					#if header line
				    if(/^>(\S+)/){;
						#extracts the accession number from the header line
						my $acc = $_;
						$acc =~ s/\ .*//;
						chomp $acc;
						substr($acc, 0, 1) = "";
						#finds the accession number file size from the hashed file
						my $grep = "grep \"$acc\" $dirVars[8]/bacteria_vs_sv.txt | head -n 1";
						my $sv = `$grep` or FileHandling::cleanDie(\@dirVars, "Can't find/open hash file.\n");
						my @svline = split(/\ /, $sv);
						#if the total size of the MINIDB with the new accession is < the biggest file size then can append the seq to file MINIDB
						if ($miniDBSize+$svline[1] < MAX_FILE_SIZE){
							$seqCount++;
							$miniDBSize += $svline[1];
						} else { #if not then close close this MINIDB and open a new MINIDB
							#checks to see if the sub has been called by matthmmer (if is true) or the updatedb script (if is false)
							if ($dirVars[1]){
								push (@miniDBStats, &closeAndCheckMiniDB(scalar(@miniDBStats)+1, $miniDBSize, $seqCount, 0, 1, \@dirVars));
							} else {
								push (@miniDBStats, &closeAndCheckMiniDB(scalar(@miniDBStats)+1, $miniDBSize, $seqCount, 0, 0, \@dirVars));
							}
							$miniDBSize = $svline[1];
							$seqCount = 1;
							#$miniDBCount++;
							my $miniDBIndex = scalar(@miniDBStats)+1;
							open(MINIDB, ">$dirVars[0]/miniDB" . $miniDBIndex . ".fa") or FileHandling::cleanDie(\@dirVars, "Could not create miniDB" . $miniDBIndex . " to write to.\n");
						}
					}
					#appends the line from xdget to the open MINIDB
					print MINIDB $_;
				}
			close(XD);
		}
	}
	#closes the last file, as all the seqs have been written
	if ($dirVars[1]){
		push (@miniDBStats, &closeAndCheckMiniDB(scalar(@miniDBStats)+1, $miniDBSize, $seqCount, 1, 1, \@dirVars));
	} else {
		push (@miniDBStats, &closeAndCheckMiniDB(scalar(@miniDBStats)+1, $miniDBSize, $seqCount, 1, 0, \@dirVars));
	}
	return @miniDBStats;
}

sub blastXDFormat{
	my $miniDBCount = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	
	for (my $i = 1; $i <= $miniDBCount; $i++){
		#runs xdformat on the miniDB
		if (system("xdformat -n -I $dirVars[0]/miniDB$i.fa > $dirVars[0]/xdformatMiniDB$i 2>&1") != 0){
			FileHandling::cleanDie(\@dirVars, "Cannot perform xdformat on the miniDB.\n");
		}
	}
}

1;
