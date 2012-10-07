#!/usr/bin/perl

use strict;
use warnings;

package ExtractSeqs;

use constant OVERLAP_SLACK => 5;

#constants detailing the column indeces of table results
use constant NHMMER_TBL_EVALUE_INDEX = 13;
use constant NHMMER_TBL_BITSCORE_INDEX = 14;
use constant NHMMER_TBL_HMMFROM_INDEX = 4;
use constant NHMMER_TBL_HMMTO_INDEX = 5;
use constant NHMMER_TBL_STRAND_INDEX = 12;
use constant NHMMER_TBL_ACCESSION_INDEX = 0;
use constant NHMMER_TBL_ALIFROM_INDEX = 6;
use constant NHMMER_TBL_ALITO_INDEX = 7;

#Object for top matching seqs, with start and end attributes for co-ords in the accession.
sub new {
    my $class = shift;
    my %arguments = @_;
    my $self = bless ( {}, $class);

	$self->{accession} = 0;
	$self->{start} = 0;
	$self->{end} = 0;
	$self->{realLength} = 0;
	$self->{bitScore} = ();
	$self->{species} = "";
	$self->{genus} = "";
	$self->{inSEED} = 0;

    while (my ($arg, $value) = each %arguments) {
        if ($arg eq "accession") {
            $self->{accession} = $value;
        } elsif ($arg eq "start") {
            $self->{start} = $value;
        } elsif ($arg eq "end") {
			$self->{end} = $value;
		} elsif ($arg eq "realLength") {
			$self->{realLength} = $value;
		} elsif ($arg eq "bitScore") {
			for (my $i = 0; $i < $value; $i++){
				$self->{bitScore}[$i] = 0.0;
			}			
		}
	}

	#print $self->{start} . " " . $self->{end} . "\n";
	return $self;
}

#returns a hash containing the sequences from the nhmmer output that have a high match against the original SEED sequence.
sub getGoodNhmmerSeqs{
	my $iter = shift;
	my $validAcc = shift;
	my @meta;
	my $seqLength;
	if ($validAcc){
		my $metaRef = shift;
		@meta = @$metaRef;		
		$seqLength = $meta[3];
	} else {
		$seqLength = shift;
	}
	my $noofTables = shift;
	my $minBitScore = shift;
	my $minPercentageLength = shift;
	my $useEvalue = shift;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $noofIters = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#decodes the nooftables var
	if ($noofTables < 0) {
		$noofTables *= -1;
	}

	my $acceptedOneSeq = 0;

	#loops through each table output
	for (my $i = 1; $i <= $noofTables; $i++){

		open FILE, "$dirVars[0]/SEED$iter-$i.nhmmer.tblout" or FileHandling::cleanDie(\@dirVars, "Can't find/open nhmmer results file for reading.\n");

			#loops through each line in the table (could be made more efficient - could terminate when min bit score reached)
			while(<FILE>){
				#excludes lines which start with # (not real rows in the table)
				if (substr($_, 0, 1) ne "#"){
					#strips out multiple spaces
					$_ =~ s/\ \ +/\ /g;
					#split the line into fields/columns
					my @lines = split(/\ /, $_);

					#if the bit score is greater than the bit score threshold and the length of the seq is greater than the theshold length
					my $satisfiesScoreCond;
					my $score;
					if (defined $useEvalue){
						$score = $lines[NHMMER_TBL_EVALUE_INDEX];
						$satisfiesScoreCond = ($score < $minBitScore);
					} else {
						$score = $lines[NHMMER_TBL_BITSCORE_INDEX];
						$satisfiesScoreCond = ($score > $minBitScore);						
					}
					if ($satisfiesScoreCond){
						if (($lines[NHMMER_TBL_HMMTO_INDEX]-$lines[NHMMER_TBL_HMMFROM_INDEX])+1 >= ($seqLength * $minPercentageLength)){
							$acceptedOneSeq = 1;
							#adjusts the co-ords of the sequence if the length is a little shorter than the original SEED
							if ($lines[NHMMER_TBL_HMMFROM_INDEX] > 1){
								if ($lines[NHMMER_TBL_STRAND_INDEX] eq "+"){
									if($dirVars[6]){ print "Adjusting start Co-ords for $lines[".NHMMER_TBL_ACCESSION_INDEX."] from $lines[".NHMMER_TBL_ALIFROM_INDEX."]"; }
									$lines[NHMMER_TBL_ALIFROM_INDEX] -= ($lines[NHMMER_TBL_HMMFROM_INDEX]-1);
									if($dirVars[6]){ print " to $lines[".NHMMER_TBL_ALIFROM_INDEX."], because this sequence length is less than $seqLength.\n"; }
								} else {
									if($dirVars[6]){ print "Adjusting start Co-ords for $lines[".NHMMER_TBL_ACCESSION_INDEX."] from $lines[".NHMMER_TBL_ALIFROM_INDEX."]"; }
									$lines[NHMMER_TBL_ALIFROM_INDEX] += ($lines[NHMMER_TBL_HMMFROM_INDEX]-1);
									if($dirVars[6]){ print " to $lines[".NHMMER_TBL_ALIFROM_INDEX."], because this sequence length is less than $seqLength.\n"; }
								}
							}
							if ($lines[NHMMER_TBL_HMMTO_INDEX] < $seqLength){
								if ($lines[NHMMER_TBL_STRAND_INDEX] eq "+"){
									if($dirVars[6]){ print "Adjusting end Co-ords for $lines[".NHMMER_TBL_ACCESSION_INDEX."] from $lines[".NHMMER_TBL_ALITO_INDEX."]"; }
									$lines[NHMMER_TBL_ALITO_INDEX] += ($seqLength-$lines[NHMMER_TBL_HMMTO_INDEX]);
									if($dirVars[6]){ print " to $lines[".NHMMER_TBL_ALITO_INDEX."], because this sequence length is less than $seqLength.\n"; }
								} else {
									if($dirVars[6]){ print "Adjusting end Co-ords for $lines[".NHMMER_TBL_ACCESSION_INDEX."] from $lines[".NHMMER_TBL_ALITO_INDEX."]"; }
									$lines[NHMMER_TBL_ALITO_INDEX] -= ($seqLength-$lines[NHMMER_TBL_HMMTO_INDEX]);
									if($dirVars[6]){ print " to $lines[".NHMMER_TBL_ALITO_INDEX."], because this sequence length is less than $seqLength.\n"; }
								}
							}
							#Checks that this line is not the original seq (same acc, same starts and ends)
							if (!@meta || (@meta && !($lines[NHMMER_TBL_ACCESSION_INDEX] eq $meta[0] && $lines[NHMMER_TBL_ALIFROM_INDEX] == $meta[1] && $lines[NHMMER_TBL_ALITO_INDEX] == $meta[2]))){
								#Generates a unique key for the seq
								my $key = $lines[NHMMER_TBL_ACCESSION_INDEX]."/".$lines[NHMMER_TBL_ALIFROM_INDEX]."-".$lines[NHMMER_TBL_ALITO_INDEX];
								#if there is not an entry in the hash with the key, then add a new entry to the hash with the data from the line
								if (!exists $seqs{$key}){
									$seqs{$key} = ExtractSeqs->new(accession => $lines[NHMMER_TBL_ACCESSION_INDEX], start => $lines[NHMMER_TBL_ALIFROM_INDEX], end => $lines[NHMMER_TBL_ALITO_INDEX], realLength => ($lines[NHMMER_TBL_HMMTO_INDEX]-$lines[NHMMER_TBL_HMMFROM_INDEX])+1, bitScore => $noofIters);
								}
								#add this iterations bit score to the hash
								$seqs{$key}->{bitScore}[$iter] = $score;
							}
						}
					} else {
						#terminates checking as min bitscore or max evalue has been reached, so no future row will not satisfy the score condition
						last;
					}
				}
			}

		close(FILE);

	}

	if (!$acceptedOneSeq && !$dirVars[5]){
		print "No sequences accepted at iteration $iter, will therefore converge...\n";
	}

	return %seqs;
}

#returns a hash containing the sequences from the blast output that have a high match against the original SEED sequence.
#this sub is very similar to the previous one but with blast rather than nhmmer
sub getGoodBlastSeqs{
	my $iter = shift;
	my $validAcc = shift;
	my @meta;
	my $seqLength;
	if ($validAcc){
		my $metaRef = shift;
		@meta = @$metaRef;		
		$seqLength = $meta[3];
	} else {
		$seqLength = shift;
	}
	my $noofTables = shift;
	my $minBitScore = shift;
	my $minPercentageLength = shift;
	my $useEvalue = shift;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $noofIters = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	if ($noofTables < 0) {
		$noofTables *= -1;
	}
	
	#loops through each table output
	for (my $i = 1; $i <= $noofTables; $i++){

		open FILE, "$dirVars[0]/SEED$iter-$i.blast.out" or FileHandling::cleanDie(\@dirVars, "Can't find/open blast results file for reading.\n");

			#loops through each line in the table
			while(<FILE>){
				#excludes lines which start with # (not real rows in the table) could change this to regexp?
				if (substr($_, 0, 1) ne "#"){
					#split the line into fields/columns
					chomp $_;
					my @lines = split(/\t/, $_);
					my $length = $lines[18] - $lines[17];
					if ($length < 0){
						$length *= -1;
					}
					$length++;

					#if the bit score is greater than the bit score threshold and the length of the seq is greater than the theshold length
					if (($lines[4] > $minBitScore) && ($length >= ($seqLength * $minPercentageLength))) {
						#adjusts the co-ords of the sequence if the length is a little shorter than the original SEED
						if ($lines[19] eq "+1"){
							if ($lines[16] eq "+1"){
								if ($lines[17] > 1){
									if($dirVars[6]){ print "Adjusting start Co-ords for $lines[1] from $lines[20]"; }
									$lines[20] -= ($lines[17]-1);
									if($dirVars[6]){ print " to $lines[20], because this sequence length is less than $seqLength.\n"; }
								}
								if ($lines[18] < $seqLength){
									if($dirVars[6]){ print "Adjusting end Co-ords for $lines[1] from $lines[21]"; }
									$lines[21] += ($seqLength-$lines[18]);
									if($dirVars[6]){ print " to $lines[21], because this sequence length is less than $seqLength.\n"; }
								}
							} else {
								if ($lines[17] < $seqLength){
									if($dirVars[6]){ print "Adjusting start Co-ords for $lines[1] from $lines[20]"; }
									$lines[20] -= ($seqLength-$lines[17]);
									if($dirVars[6]){ print " to $lines[20], because this sequence length is less than $seqLength.\n"; }
								}
								if ($lines[18] > 1){
									if($dirVars[6]){ print "Adjusting end Co-ords for $lines[1] from $lines[21]"; }
									$lines[21] += ($lines[18]-1);
									if($dirVars[6]){ print " to $lines[21], because this sequence length is less than $seqLength.\n"; }
								}
								my $temp = $lines[20];
								$lines[20] = $lines[21];
								$lines[21] = $temp;						
							}
						} else {
							warn "-1 in subject frame, this has not been accounted for and so may cause future errors.\n";
						}
						#Checks that this line is not the original seq (same acc, same starts and ends)
						if (!@meta || (@meta && !($lines[1] eq $meta[0] && $lines[20] == $meta[1] && $lines[21] == $meta[2]))){
							#Generates a unique key for the seq
							my $key = $lines[1]."/".$lines[20]."-".$lines[21];
							#if there is not an entry in the hash with the key, then add a new entry to the hash with the data from the line
							if (!exists $seqs{$key}){
								$seqs{$key} = ExtractSeqs->new(accession => $lines[1], start => $lines[20], end => $lines[21], realLength => $length-1, bitScore => $noofIters);
							}
							#add this iterations bit score to the hash
							$seqs{$key}->{bitScore}[$iter] = $lines[4];
							#checks to see if there is an overlap with other seqs and removes less good seq. (Could do slightly more efficiently?)
							&checkAndResolveOverlaps(\%seqs, $key, $iter, \@dirVars);
						}
					}
				}
			}

		close(FILE);

	}

	return %seqs;	
}

#returns the SEEDs meta data - accession no (with version), start, and end co-ords.
sub getMetaFromSEED{
	my $validAcc = shift;
	my $alignmentFile = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	my @meta;
	my $header;
	my $seedSeq = "";

	my $firstSeq = 1;

	open (SEED, "SEED") or die "Unable to open SEED\n";

		while(<SEED>){
			chomp $_;
			if($alignmentFile){
				#stockholm format, finds first accession line
				if(!/^#/ && !/^$/){
					$header = $_;
					$header =~ s/\s.*//;
					last;
				}
			} else {
				if(/^>/){
					#fasta format finds first accession line
					if ($firstSeq){
						$header = $_;
						$header = substr($header, 1);
						$firstSeq = 0;
					} else {
						last;
					}
				} else {
					#fasta first seq
					$seedSeq .= $_;
				}				
			}
		}

	close (SEED);

	$seedSeq =~ s/\s+//g;

	#if not a junk header (i.e. more validation needs to be performed)
	if ($validAcc){
		#splits the header line up into the appropriate segments.
		my @metas = split('/', $header);
		push (@meta, $metas[0]);
		push (@meta, split('-', $metas[1]));
		push (@meta, $meta[1]-$meta[2]);
		if ($meta[3] < 0){
			$meta[3] *= -1;
		}
		#print $meta[0] . " " . $meta[1] . " " . $meta[2] . " " . $meta[3] . "\n";

		if(!$alignmentFile){
			#removes the accession version if this is given with the accession
			my $accession = Taxon::stripAccessionVersion($meta[0]);
			#converts the accession into equivalent taxonomy ids
			my %accessionHash = Taxon::getTaxonIDsByAccessionNo([($accession)], \@dirVars);
			#takes the first tax id returned for the accession (even if there is multiple tax ids)
			if (!defined $accessionHash{$accession}[0]){
				FileHandling::cleanDie(\@dirVars, "$accession is not a valid accession number in the SEED file header.\n");
			}

			my @accession = ($accession);
			my $db = CreateMiniDB::findDatabases(\@accession, 1, \@dirVars);

			#gets seq from db for acc and co-ords provided
			my $seq = "";
			open (XD, "xdget -n $dirVars[8]/$db $meta[0] -a $meta[1] -b $meta[2] |") or FileHandling::cleanDie(\@dirVars, "Could not open xdget pipe:[$!]\n");
				while(<XD>){
					if (!/^>/){
						chomp $_;
						$seq .= $_;
					}
				}
			close (XD);	

			$seq =~ s/\s+//g;

			$seedSeq =~ tr/[a,c,g,t,u]/[A,C,G,T,T]/;
			$seq =~ tr/[a,c,g,t,u]/[A,C,G,T,T]/;

			$seedSeq =~ s/U/T/g;
			$seq =~ s/U/T/g;

			#checks the seqs are the same, if not die
			if ($seedSeq ne $seq){
				FileHandling::cleanDie(\@dirVars, "SEED file sequence does not match the databases sequence of:\n$seq\n");
			}
		}

		return @meta;
	} else {
		#only meta that can be returned is the first sequence length
		return length($seedSeq);
	}
}

#used only for BLAST
sub checkAndResolveOverlaps{
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $newSeqKey = shift;
	my $iter = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;
	
	foreach my $key (keys %seqs){
		if ($key ne $newSeqKey && $seqs{$key}->{accession} eq $seqs{$newSeqKey}->{accession} && defined($seqs{$key}->{bitscore}[$iter])){
			my $start1 = $seqs{$key}->{start};
			my $end1 = $seqs{$key}->{end};
			my $start2 = $seqs{$newSeqKey}->{start};
			my $end2 = $seqs{$newSeqKey}->{end};

			if ($start1 > $end1){
				if (($start1 - $end1) > (2*OVERLAP_SLACK)){
					$start1 -= OVERLAP_SLACK;
					$end1 += OVERLAP_SLACK;
				}
			} else {
				if (($end1 - $start1) > (2*OVERLAP_SLACK)){
					$start1 += OVERLAP_SLACK;
					$end1 -= OVERLAP_SLACK;
				}
			}

			my $points = 0;

			if ($start1 > $start2){ $points++; }
			if ($start1 > $end2){ $points++; }
			if ($end1 > $end2){ $points++; }
			if ($end1 > $start2){ $points++; }

			if ($points > 0 && $points < 4){
				#overlap has occured
				print "overlap\n";
				if ($seqs{$key}->{bitscore}[$iter] > $seqs{$newSeqKey}->{bitscore}[$iter]){
					if ($iter > 0){
						if ($seqs{$newSeqKey}->{bitscore}[$iter-1] == 0.0){
							delete($seqs{$newSeqKey});
						} else {
							$seqs{$newSeqKey}->{bitscore}[$iter] = -1;
						}
					} else {
						delete($seqs{$newSeqKey});
					}					
				} else {
					if ($iter > 0){
						if ($seqs{$key}->{bitscore}[$iter-1] == 0.0){
							delete($seqs{$key});
						} else {
							$seqs{$key}->{bitscore}[$iter] = -1;
						}
					} else {
						delete($seqs{$key});
					}
				}
				last;
			}
		}
	}
}

#prints seq into the new SEED
sub addSeq{
	my $seq = shift;
	my $iter = shift;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#finds which db the seq is in
	my @accession = ($seq->{accession});
	my $db = CreateMiniDB::findDatabases(\@accession, 1, \@dirVars);
	#xdgets the seq based on the coords provided
	#print "$seq->{accession} $seq->{start} $seq->{end}";
	open (XD, "xdget -n $dirVars[8]/$db $seq->{accession} -a $seq->{start} -b $seq->{end} |") or FileHandling::cleanDie(\@dirVars, "Could not open xdget pipe:[$!]\n");
		while(<XD>){
			if (substr($_, 0, 1) eq ">"){
				print SEED ">$seq->{accession}/$seq->{start}-$seq->{end}\n";
			} else {
				print SEED $_;
			}
		}
	close (XD);
}

#fetches the sequences from the dbs and adds them to the SEED could maybe merge the nhmmerFetchSeqs with blastFetchSeqs??
sub nhmmerFetchSeqs{
	my $iter = shift;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	#copies the original SEED sequence into the new SEED rather than refetching for efficiency
	FileHandling::copySEED(\@dirVars, $iter);
	open (SEED, ">> $dirVars[0]/SEED$iter ") or FileHandling::cleanDie(\@dirVars, "Could not create SEED$iter to write to\n");

		#loops through each seq in hash adding it to the new SEED
		foreach my $seq (keys %seqs){
			#print $seq . " " . $seqs{$seq}->{start} . " " . $seqs{$seq}->{end} . "\n";
			if ($iter > 1){
				#only adds a seq that is new this iter
				if ($seqs{$seq}->{bitScore}[$iter-1] > 0 && $seqs{$seq}->{bitScore}[$iter-2] == 0){
					&addSeq($seqs{$seq}, $iter, \@dirVars);
				}
			} else {
				#adds all seqs as its first iter
				&addSeq($seqs{$seq}, $iter, \@dirVars);
			}
		}

	close(SEED);
}

#fetches the sequences from the dbs and adds them to the SEED
sub blastFetchSeqs{
	my $iter = shift;
	my $seqsRef = shift;
	my %seqs = %$seqsRef;
	my $dirVarsRef = shift;
	my @dirVars = @$dirVarsRef;

	open (SEED, "> $dirVars[0]/SEED$iter ") or FileHandling::cleanDie(\@dirVars, "Could not create SEED$iter to write to\n");

		#loops through each seq in hash adding it to the new SEED
		foreach my $seq (keys %seqs){
			#print $seq . " " . $seqs{$seq}->{start} . " " . $seqs{$seq}->{end} . "\n";
			if ($seqs{$seq}->{bitScore}[$iter-1] >= 0){
				&addSeq($seqs{$seq}, $iter, \@dirVars);
			}
		}

	close(SEED);
}

1;
