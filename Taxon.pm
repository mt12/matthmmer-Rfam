#!/usr/bin/perl

use strict;
use warnings;

package Taxon;

#Returns all the taxonomy ids for an array of accession nos (or all the accession nos for an array of taxonomy ids, depending on arguments given)
sub createHash{
	my $indexesRef = shift;
	my @indexes = @$indexesRef;
	#if flag is -a then a taxonomy id has been given, if flag is -t then an accession no has be given
	my $flag = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	open FILE, "$dirVars[7]/files/pairs" or die "Can't find/open pairs file for reading\n";;

    	my %accessionHash;

	    while(<FILE>){
	        my $accessionNo = $_;
	        chomp $accessionNo;
	        my $taxonID = $accessionNo;
			#splits the pairs line into the accessionNo => taxonID
	        $accessionNo =~ s/[^A-Z0-9].*$//;
	        $taxonID =~ s/.*=>\ \s*//;
			#hash built on taxonIDs
	        if ($flag eq '-a'){
				#if the taxon id is in the list of taxon ids required
				if (grep(/^$taxonID$/, @indexes)) {
					#if the taxon id does already exist in the hash then push corresponding accession no to the array
					if (exists $accessionHash{$taxonID}) {
			            push (@{$accessionHash{$taxonID}}, $accessionNo);
			        } else { #taxon id doesnt exist in the hash then add a new entry to the hash with the corresponding accession no
			            $accessionHash{$taxonID}[0] = $accessionNo;
			        }
				}
	        } elsif ($flag eq '-t') { #hash built on accessionNos
				#if the taxon id is in the list of accession nos required
				if (grep(/^$accessionNo$/, @indexes)) {	
					#if the accession no does already exist in the hash then push corresponding accession no to the array			
			        if (exists $accessionHash{$accessionNo}) {
			            push (@{$accessionHash{$accessionNo}}, $taxonID);
			        } else { #taxon id doesnt exist in the hash then add a new entry to the hash with the corresponding accession no
						$accessionHash{$accessionNo}[0] = $taxonID;
					}
				}
			} else {
				die "Neither -a or -t was set.\n";
			}
		}

	close (FILE);

	return %accessionHash;
}

#gets all accession nos that correspond to a given taxonomy id(s)
sub getAccessionNosByTaxonID{
	my $taxonIDsRef = shift;
	my @taxonIDs = @$taxonIDsRef;
	my $wantVersionNo = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	#converts the accession nos into the accessionNo.version
	if ($wantVersionNo){
		my %accessionNos = &createHash(\@taxonIDs, '-a', \@dirVars);
		foreach my $accessionNo (keys %accessionNos){
			foreach my $accession (@{$accessionNos{$accessionNo}}){
				#finds the accession in the hash file
				my $grep = "grep \"$accession\" $dirVars[8]/bacteria_vs_sv.txt | head -n 1";
				my $sv = `$grep` or die "Hash file does not exist\n";
				my @svline = split(/\ /, $sv);
				$accession = $svline[0];
			}
		}
		return %accessionNos;
	} else {
		return &createHash(\@taxonIDs, '-a', \@dirVars);
	}
}

#gets all taxonomy ids that correspond to a given accession no(s)
sub getTaxonIDsByAccessionNo{
	my $accessionNosRef = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my @accessionNos = @$accessionNosRef;
	return &createHash(\@accessionNos, '-t', \@dirVars);
}

#checks the pairs file to see if the taxonomy id exists.
sub taxonIDExists{
	my $taxonID = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my $grep = "grep \"=> " . $taxonID . "\$\" $dirVars[7]/files/pairs | wc -l";
    my $lines = `$grep` or FileHandling::cleanDie(\@dirVars, "Can't find/open pairs file\n");
	chomp $lines;
	return $lines;
}

#checks that the taxonomy id is present in newNodes.dmp
sub taxIDInNodesdmp{
	my $taxonID = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my $grep = "grep \"" . $taxonID . "\$\" $dirVars[7]/files/taxIdsNotInNodesdmp | wc -l";
    my $lines = `$grep` or FileHandling::cleanDie(\@dirVars, "Can't find/open taxIdsNotInNodesdmp file \n");
	chomp $lines;
	return $lines;
}

#check parent id is in newNodes.dmp
sub parentIDInNodesdmp{
	my $taxonID = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my $grep = "grep \"^".$taxonID."[[:blank:]]*|\" $dirVars[7]/files/newNodes.dmp | wc -l";
	my $lines = `$grep` or FileHandling::cleanDie(\@dirVars, "Can't find/open taxIdsNotInNodesdmp file \n");
	chomp $lines;
	return $lines;
}

#checks that the taxonomy id given is a valid id
sub taxonIDValid{
	my $taxonID = shift;
	my $isParent = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	if ($isParent){
		if (&parentIDInNodesdmp($taxonID, \@dirVars)){
			return 1;
		} else {
			return 0;
		}
	} else {
		if (&taxonIDExists($taxonID, \@dirVars) && !&taxIDInNodesdmp($taxonID, \@dirVars)){
			return 1;
		} else {
			return 0;
		}
	}
}

#removes the version number from a given accession no (returns the original accession no if no version id is given
sub stripAccessionVersion{
	my $accession = shift;
	$accession =~ s/[^A-Z0-9].*$//;
	return $accession;
}

#gets the taxonomy string the rfam_11_0 database for a given ncbi id
sub getTaxonStr{
	my $ncbiID = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;	

	#executes the query
	$dirVars[3]->execute($ncbiID);
	if( $DBI::err ) {
		print STDERR "(WW) WARNING: error executing query asth to get data: "
		    . $dirVars[2]->errstr . "\n";
	}

	#gets the result of the querys
	my $result = $dirVars[3]->fetchrow_array();
	if( $dirVars[3]->err ) {
		print STDERR "(WW) WARNING: error whilst retrieving query asth"
		    . $dirVars[2]->errstr . "\n";
	}

	return $result;
}

1;
