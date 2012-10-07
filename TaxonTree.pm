#!/usr/bin/perl

use strict;
use warnings;

#require FileHandling;

package TaxonTree;

#Object for all nodes, with a parent and group attribute.
sub new {
    my $class = shift;
    my %arguments = @_;
    my $self = bless ( {}, $class);

	$self->{parent} = "";
	$self->{group} = "";

    while (my ($arg, $value) = each %arguments) {
        if ($arg eq "parent") {
            $self->{parent} = $value;
        } elsif ($arg eq "group") {
			$self->{group} = $value;
		}
	}

	#print $self->{parent} . " " . $self->{group} . "\n";
	return $self;
}

#reads in all the bacteria nodes and returns a hash of taxonomy ID -> parent taxonomy ID
sub readNodesDump{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my %pairs;
	my $nodefile = "$dirVars[7]/files/newNodes.dmp";
	if ($dirVars[0]){
		open (NODES, "<$nodefile") or FileHandling::cleanDie(\@dirVars, "Can't open/find newNodes.dmp file for reading.\n");
	} else {
		open (NODES, "<$nodefile") or die "Can't open/find newNodes.dmp file for reading.\n";
	}

	    while (<NODES>) {
			my $linein = $_;
			chomp $linein;
			$linein =~ s/\s+//g;
			my @line = split(/\|/, $linein);
			#print $line[0] . " => ". $line[1] . "\n";
			$pairs{$line[0]} = TaxonTree->new(parent => $line[1], group => $line[2]);
	    }

	close (NODES);
	
	return %pairs;
}

#recursive subroutine to find the parent of a given taxonomy id by the parents group
sub getParentByGroup{
	my $child = shift;
	my $group = shift;
	my $forReport = shift;
	my $pairsRef = shift;
	my %pairs = %$pairsRef;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	if ($pairs{$child}->{group} ne $group && $child != 1){
	    #get the parent of this parent
		&getParentByGroup($pairs{$child}->{parent}, $group, $forReport, \%pairs, \@dirVars);
	} else {
	    if ($child == 1){
			if ($forReport){
				return 0;
			} else {
				FileHandling::cleanDie(\@dirVars, "Group entered does not exist as a parent of the given child.\n");
			}
	    } else {
		    #at the nth parent so return this taxonomy id
		    return $child;
		}
	}
}

#given a taxonomy id and depth, returns all the siblings/cousins for taxonomy id from the depth given.
sub getSiblings{
	my $node = shift;
	my $nodeIsParent = shift;
	my $depth = shift;
	my $group = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#recursive subroutine to retrieve all the children (leafs only) of a given parent taxonomy id;
    sub getChildren{
		my $parent = shift;
		my $pairsRef = shift;
		my %pairs = %$pairsRef;
        my @leaves;
        my $noofChildren = 0;
        foreach my $child (keys %pairs){
            if ($pairs{$child}->{parent} == $parent){
                $noofChildren++;
				#returns the grandchildren of the parent taxonomy id, by calling this recursively, this retrieves all of the leaves from the branch
                my @leaf = &getChildren($child, $pairsRef);
				#adds the leaves returned (could be multiple despite the variable name!) to the overall collection of leaves
                push(@leaves, @leaf);
            }
        }
        if ($noofChildren == 0){
			#if this taxonomy id has no children, then return this taxonomy id, as this is a leaf
            return $parent;
        } else {
			#if this taxonomy id has children, then this is not a leaf, so return the leaves from this branch (the leaves are filtering up the tree)
            return @leaves;
        }
    }

	#recursive subroutine to find the nth parent of a given taxonomy id
	sub getParentByDepth{
	    my $child = shift;
	    my $count = shift;
		my $pairsRef = shift;
		my %pairs = %$pairsRef;
	    if ($count > 0 && $child != 2){
			#get the parent of this parent
			&getParentByDepth($pairs{$child}->{parent}, --$count, $pairsRef);
	    } else {
			if ($child == 2){
				warn "Depth given is larger than the full tree depth. All children will be returned...\n";
			}
			#at the nth parent so return this taxonomy id
			#print $child . "\n";

			return $child;
	    }
	}

	#reads in the file
	my %pairs = &readNodesDump(\@dirVars);
	
	#finds all of its children from the given parent and returns them in an array
	if ($nodeIsParent){
		return &getChildren($node, \%pairs);
	} else { #finds the parent and then finds all of its children and returns them in an array
		my $parent;
		if ($depth){
			$parent = &getParentByDepth($node, $depth, \%pairs);
		} else {
			$parent = &getParentByGroup($node, $group, 0, \%pairs, \@dirVars);
		}
		if ($parent == 2){
			return (-1);
		} else {
			if ($dirVars[0]){
				FileHandling::addTreeInfo(\@dirVars, $parent, $pairs{$parent}->{group});		
			}
			return &getChildren($parent, \%pairs);
		}
	}
}

#gets the genus or species for each sequence (used for creating reports)
sub getAllSeqsParents{
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $noofIters = shift;
	my $level = shift;

	my %pairs = &readNodesDump(\@dirVars);
	my %levels;
	
	foreach my $seq (keys %seqs){
		my $accession = Taxon::stripAccessionVersion($seqs{$seq}->{accession});
		my %accessionHash = Taxon::getTaxonIDsByAccessionNo([($accession)], \@dirVars);
		#takes the first tax id returned for the accession (even if there is multiple tax ids)
		if (defined $accessionHash{$accession}[0]){
			$seqs{$seq}->{$level} = &getParentByGroup($accessionHash{$accession}[0], $level, 1, \%pairs, \@dirVars);
			if ($seqs{$seq}->{$level} ne 0){
				if (exists $levels{$seqs{$seq}->{$level}}){
					for (my $i = 0;$i < $noofIters;$i++){
						if ($seqs{$seq}->{bitScore}[$i] > 0){
							$levels{$seqs{$seq}->{$level}}[$i]++;
						}
					}
				} else {
					for (my $i = 0;$i < $noofIters;$i++){
						if ($seqs{$seq}->{bitScore}[$i] > 0){
							$levels{$seqs{$seq}->{$level}}[$i] = 1;
						} else {
							$levels{$seqs{$seq}->{$level}}[$i] = 0;
						}
					}
				}
			}
		} else {
			FileHandling::cleanDie(\@dirVars, "Accession found from search does not have a valid taxonomy id.\n");
		}
	}

	FileHandling::createLevelList(\@dirVars, $level, $noofIters, \%levels);
	
}

1;
