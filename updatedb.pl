#!/usr/bin/perl

use strict;
use warnings;

#push (@INC, "/nfs/users/nfs_m/mt12/scripts");

require Taxon;
require TaxonTree;
require CreateMiniDB;

use constant SCRATCH_DIR => "/warehouse/pfam01/rfam/Users/mt12";
use constant DBS_DIR => "/lustre/scratch101/blastdb/Rfam/rfamseq/EMBL_Bacterial_Genomes";
use constant MAX_DB_SIZE => 3500000000;

my $newCount = 0;

#opens for reading the new list of accessions
open (NEWACCS, "accession_nos") or die "Can't find the new accession numbers file\n";

#opens for appending the list of accessions and size
open (BACTERIA_VS_SV, ">>" . DBS_DIR . "/bacteria_vs_sv.txt") or die "Cant open bacteria_vs_sv file for ammending to\n";

#opens for appending the accession-taxid link file
open (PAIRS, ">>" . SCRATCH_DIR . "/files/pairs") or die "Cant open pairs file for ammending to\n";

#opens for appending the list of taxonomy ids
open (TAXONIDS, ">>" . SCRATCH_DIR . "/files/taxonIDs") or die "Cant open taxonIDs file for ammending to\n";

	while(<NEWACCS>){
		my $newAcc = $_;
		chomp $newAcc;
		if (!-e DBS_DIR."/fastavfiles/$newAcc.fa"){	#for testing using results, should change to fastavfiles AND .embl to .fa
			print $newAcc . "\n";
			$newCount++;

			my $dir = DBS_DIR."/fastavfiles";
			if(system("pfetch $newAcc > $dir/$newAcc.fa") != 0){
				die "Couldn't pfetch $newAcc\n";
			}
            my $accv = `head -n 1 $dir/$newAcc.fa`;
            $accv =~ s/\ .*//;
            chomp $accv;
            substr($accv, 0, 1) = "";
			print BACTERIA_VS_SV $accv . " " . &getResidues("$dir/$newAcc.fa") . "\n";
			
			my $pfetchID = `pfetch -F $newAcc | grep '/db_xref="taxon:'`;
			my @taxonids = split("\n", $pfetchID);
			foreach my $taxonid (@taxonids){
				chomp $taxonid;
				$taxonid =~ s/.*\/db_xref="taxon://;
				$taxonid =~ s/"//;
				print PAIRS $newAcc . " => " . $taxonid . "\n";
				if(system("grep $taxonid " . SCRATCH_DIR . "/files/taxonIDs") != 0){
					print TAXONIDS $taxonid . "\n";
				}
			}
		}
	}

close (TAXONIDS);

close (PAIRS);

close (BACTERIA_VS_SV);

close (NEWACCS);

#sorts the 3 files appended
system("sort " . DBS_DIR . "/bacteria_vs_sv.txt > " . DBS_DIR . "/bacteria_vs_sv_sorted.txt");
system("sort " . SCRATCH_DIR . "/files/pairs > " . SCRATCH_DIR . "/files/newPairs");
system("sort " . SCRATCH_DIR . "/files/taxonIDs > " . SCRATCH_DIR . "/files/newTaxonIDs");

system("mv " . DBS_DIR . "/bacteria_vs_sv_sorted.txt " . DBS_DIR . "/bacteria_vs_sv.txt");
system("mv " . SCRATCH_DIR . "/files/newPairs " . SCRATCH_DIR . "/files/pairs");
system("mv " . SCRATCH_DIR . "/files/newTaxonIDs " . SCRATCH_DIR . "/files/taxonIDs");

my @taxonIDs;

#opens list of taxon id file
open (TAXONIDS, SCRATCH_DIR . "/files/taxonIDs") or die "Can't open taxonids file\n";

	while(<TAXONIDS>){
		my $taxID = $_;
		chomp $taxID;
		push (@taxonIDs, $taxID);
	}

close (TAXONIDS);

my %pairs;
my %lines;

#opens the massive nodes.dmp file
open (NODES, SCRATCH_DIR . "/files/nodes.dmp") or die "Can't open taxonomy node file for reading!!!!\n";

    while (<NODES>) {
		    my $linein = $_;
		    chomp $linein;
	    	$linein =~ s/\s+//g;
		    my @line = split(/\|/, $linein);
			$lines{$line[0]} = $_;
			$pairs{$line[0]} = $line[1];
	}

close (NODES);

my %linesToKeep;

sub recTree{
	my $nodeID = shift;
	if (!exists $linesToKeep{$nodeID} && exists $pairs{$nodeID}){
		#print $nodeID . "\n";
		$linesToKeep{$nodeID} = $lines{$nodeID};
		#print $pairs{$nodeID} . "\n";
		if ($nodeID ne 1){	#not been tested with 1 yet
			&recTree($pairs{$nodeID});
		}
	}
}

foreach (@taxonIDs){
	&recTree($_);	
}

#write to newNodes
open (NEWNODES, ">" . SCRATCH_DIR . "/files/newNodes.dmp") or die "Cant create newnodes file\n";

	foreach my $line (sort {$a <=> $b} keys %linesToKeep){
		print NEWNODES $linesToKeep{$line};
	}

close (NEWNODES);

#open for reading taxon list
open (TAXONIDS, SCRATCH_DIR . "/files/taxonIDs") or die "Can't open taxonids file\n";

my $firstTaxon = 0;

#write a new list of taxon ids not in dmp
open (TAXIDSNOTINDUMP, ">" . SCRATCH_DIR . "/files/taxIdsNotInNodesdmp") or die "Can't create taxIdsNotInNodesdmp file\n";

	while(<TAXONIDS>){
		my $taxID = $_;
		chomp $taxID;

		if ($firstTaxon != 0){
			$firstTaxon = $taxID;
		}

		my $grstr = "\"^$taxID\t|\"";

		if(system("grep $grstr " . SCRATCH_DIR . "/files/newNodes.dmp") != 0){
    	    print TAXIDSNOTINDUMP $taxID . "\n";
		}
	}

close (TAXIDSNOTINDUMP);

close (TAXONIDS);


open (BACTERIA_VS_SV, DBS_DIR . "/bacteria_vs_sv.txt") or die "Cant open bacteria_vs_sv file\n";

open (SEQINDEX, ">" . DBS_DIR . "/seqIndex/seqIndex.txt") or die "Cant create seqIndex file\n";

	my $dbSize = 0;
	my $dbCount = 1;
	my $first = 1;
	my $lastAcc = "";
	my @files;

	print SEQINDEX "DB1.fa\t";

	while(<BACTERIA_VS_SV>){
		my $thisAcc = $_;
		chomp $thisAcc;
		$thisAcc =~ s/[^A-Z0-9].*$//;
		if ($first){
			print SEQINDEX "$thisAcc\t";
			$first = 0;
		}
		my $residues = $_;
		chomp $residues;
		$residues =~ s/.*\s//;
		if ($dbSize + $residues > MAX_DB_SIZE){
			#&catFiles(\@files, $dbCount);
			$dbCount++;
			print SEQINDEX 	 "$lastAcc\n"
							."DB$dbCount.fa\t$thisAcc\t";
			$dbSize = $residues;
			@files = ();
			push(@files, $thisAcc);
		} else {
			push(@files, $thisAcc);
			$dbSize += $residues;
		}
		$lastAcc = $thisAcc;		
	}

close (BACTERIA_VS_SV);

print SEQINDEX "$lastAcc\n";
#&catFiles(\@files, $dbCount);

close (SEQINDEX);

print $newCount . "\n";

sub getResidues{
    my $file = shift;
    my @seqstat = `seqstat --quiet $file`;
    my $residues = $seqstat[3];
    $residues =~ s/.*:\s*//;
    chomp $residues;
    return $residues;
}

sub catFiles{
	my $filesRef = shift;
    my @files = @$filesRef;
	my $dbCount = shift;
	my $catStr = "cat ";
    foreach (@files){
        $catStr .= DBS_DIR . "/fastavfiles/" . $_ . ".fa ";
    }
    $catStr .= "> " . DBS_DIR . "/DB" . $dbCount . ".fa";
    system($catStr);
	print "DB $dbCount has been written\n";
	system("xdformat -n -I " . DBS_DIR . "/DB$dbCount.fa");
}

#remakes the bigMiniDBs (effectively runs code from matthmmer


my @dirVars = (0, 0, 0, 0, 0, 0, 0, SCRATCH_DIR, DBS_DIR);

my @siblings = TaxonTree::getSiblings(2, 1, 0, 0, \@dirVars);

my @accessions;

#converts each siblings taxonomy id to the accession no(s) (with version no).
my %accessionsHash = Taxon::getAccessionNosByTaxonID(\@siblings, 1, \@dirVars);	
foreach my $accessionNos (keys %accessionsHash){
	push (@accessions, @{$accessionsHash{$accessionNos}});
}

#finds which database each accession no is in and returns this as a hash
my %dbs = CreateMiniDB::findDatabases(\@accessions, 0, \@dirVars);

#creates the minidb(s) from the accession nos
$dirVars[0] = DBS_DIR . "/bigMiniDBs";
my @miniDBStats = CreateMiniDB::fetchSeqs(\%dbs, \@dirVars);

for (my $i = 1; $i <= scalar(@miniDBStats); $i++){
	#runs xdformat on the miniDB
	if (system("xdformat -n -I $dirVars[0]/miniDB$i.fa > $dirVars[0]/xdformatMiniDB$i 2>&1") != 0){
		die "Cannot perform xdformat on the miniDB.\n";
	}
}

system("cp accession_nos " . SCRATCH_DIR . "/files/");

print "All done!!\n";

1;
