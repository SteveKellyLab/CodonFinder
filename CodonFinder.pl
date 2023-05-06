#! /usr/bin/perl -w

use Getopt::Std;

&run_program;

sub run_program	{
				%codons = ();
				&def_param;
				&get_options;
				&set_starting_scores;
				&read_sequences($seq_file);								
				&score_all_genes(0);							
				&optimise_new;			
				}
				
sub exit_screen	{
				print "\n";
				print "Steve's Codon Optimiser\n";
				print "optimise_codons.pl -i [input FASTA file] -w [weights file]\n\n";
				exit;
				}
				
sub optimisation_criteria	{
							my $local_CAI = $_[0];
							my $local_ENT = $_[1];
							my $dist = $_[2];
							my $result = $local_CAI*$local_CAI*$local_ENT*$dist;
							return $result;
							}				
				
sub get_options	{
				&getopts('i:w:',\%parameters);
				if (exists $parameters{'i'})
					{
					$seq_file = $parameters{'i'};
					chomp $seq_file;					
					}				
				else
					{
					&exit_screen;
					}
				if (exists $parameters{'w'})
					{
					$weights_file = $parameters{'w'};
					chomp $weights_file;
					$weights_start = 1;					
					}				
				else
					{
					$weights_start = 0;
					&exit_screen
					}
				}

sub def_param	{
				$mut_rate = 0.1;
				$generations = 20;
				$population_size = 400;
				$founder_size = 20;
				%synon = (TTT, F, TCT, S, TAT, Y, TGT, C, TTC, F, TCC, S, TAC, Y, TGC, C, TTA, L, TCA, S, TAA, STOP, TGA, STOP, TTG, L, TCG, S, TAG, STOP, TGG, W, CTT, L, CCT, P, CAT, H, CGT, R, CTC, L, CCC, P, CAC, H, CGC, R, CTA, L, CCA, P, CAA, Q, CGA, R, CTG, L, CCG, P, CAG, Q, CGG, R, ATT, I, ACT, T, AAT, N, AGT, S, ATC, I, ACC, T, AAC, N, AGC, S, ATA, I, ACA, T, AAA, K, AGA, R, ATG, M, ACG, T, AAG, K, AGG, R, GTT, V, GCT, A, GAT, D, GGT, G, GTC, V, GCC, A, GAC, D, GGC, G, GTA, V, GCA, A, GAA, E, GGA, G, GTG, V, GCG, A, GAG, E, GGG, G);
				%data = ();
				%temp_count = ();				
				
				foreach my $codon (keys(%synon))
					{
					unless ($synon{$codon} eq 'STOP')
						{
						$data{$synon{$codon}}{$codon} = 1;
						$temp_count{$synon{$codon}}{$codon} = 0;
						}
					}	
				@aminos = sort {$a cmp $b} keys %data;					
				}
				
sub reset_data	{
				%temp_sequence = ();
				%temp_codons = ();
				}
				
sub evaluate_sequence	{
						my $local_sequence = $_[0];
						&reset_data;
						&read_seq($local_sequence);						
						my $cai = &compute_nCAI_fast;
						my $ent = &compute_entropy_fast_2;
						my $dist = &calc_dist($local_sequence);
						my $v = &optimisation_criteria($cai, $ent, $dist);	
						if ($v > $best)
							{
							$best = $v;
							$MAX_CAI = $cai;
							$MAX_ENT = $ent;
							%best_sequence = %temp_sequence;
							%best_codons = %temp_codons;							
							$pass1 = 1;
							}
						return $pass1;
						}
						
sub evaluate_sequence_2	{
						my $local_sequence = $_[0];
						&reset_data;
						&read_seq($local_sequence);	
						my $dist = &calc_dist($local_sequence);
						my $cai = &compute_nCAI_fast;
						my $ent = &compute_entropy_fast_2;								
						my $v = &optimisation_criteria($cai, $ent, $dist);	
						return $v;
						}
			
sub generate 	{
				my $local_seq = $_[0];
				my $return_seq = '';
				while ($local_seq =~ m/(\w{3})/gmis)
					{
					my $codon = $1;
					if (($codon eq 'TGA')||($codon eq 'TAG')||($codon eq 'TAA'))
						{
						$return_seq .= $codon;
						}
					else
						{
						$return_seq .= $lookup{$synon{$codon}}{int(rand($draw_max{$synon{$codon}}))};						
						}
					}
				return $return_seq;
				}
				
sub count_codons	{
					my $local_seq = $_[0];
					my $return_seq = '';
					my %local_count = ();
					while ($local_seq =~ m/(\w{3})/gmis)
						{
						my $codon = $1;
						unless (($codon eq 'TGA')||($codon eq 'TAG')||($codon eq 'TAA'))
							{								
							++$temp_count{$synon{$codon}}{$codon};
							}
						}
					open OUTX2, ">Codon_usage_profile.txt";
					print OUTX2 ">AA\tCodon\tCount\tWeight\n";
					foreach my $aa (keys(%temp_count))
						{
						foreach my $codon (keys(%{$temp_count{$aa}}))
							{
							print OUTX2 "$aa\t$codon\t$temp_count{$aa}{$codon}\t$START{$codon}\n";
							}							
						}
					close OUTX2;
					}
				
sub test_sequence	{
					my $local_seq = $_[0];
					my $return = 0;
					#my %banned = (BpiI, GAAGAC, BsaI, GGTCTC, Esp3i, CGTCTC);	
					#DraIII	CACNNNGTG					
					if ($return == 0)
						{
						if (($local_seq =~ m/AAAAA/)||($local_seq =~ m/CCCCC/)||($local_seq =~ m/GGGGG/)||($local_seq =~ m/TTTTT/))
							{
							$return = 1;
							}
						}					
					return $return;
					}
					
sub generate_starting_sequence	{
								my $local_seq = $_[0];
								my $test_result = 1;
								my $start_seq = '';
								while ($test_result == 1)
									{
									$start_seq = &generate($local_seq);
									#$start_seq = &randomise_codons($start_seq);									
									$test_result = &test_sequence($start_seq);									
									}
								return $start_seq;
								}
								
sub mutate_sequence	{
					my $local_seq = $_[0];
					my @local_seq = ($local_seq =~ m/.{3}/g);					
					my $mutations = int((length $local_seq)*($mut_rate));
					my $i = 0;
					while ($i <= $mutations)
						{
						my $rand_1 = int(rand($#local_seq));					
						my $codon = $local_seq[$rand_1];
						my $rand_2 = int(rand($draw_max{$synon{$codon}}));
						my $replacement = $lookup{$synon{$codon}}{$rand_2};					
						$local_seq[$rand_1] = $replacement;
						++$i;
						}
					my $return_seq = join '', @local_seq;
					my $shuff = rand(1);
					if ($shuff > 0.95)
						{
						$return_seq = &randomise_codons($return_seq);
						}
					return $return_seq;					
					}
					
sub recombine_parents	{
						my $parent_1 = $_[0];
						my $parent_2 = $_[1];
						my $length = length($parent_1);
						my $recombination_point = 3*(int(rand(length($parent_1)/3)));						
						my $part_1 = substr($parent_1, 0, $recombination_point);
						my $part_2 = substr($parent_2, $recombination_point, $length - $recombination_point);
						my $hybrid = "$part_1"."$part_2";
						return $hybrid;
						}
					
sub evolve_sequence	{
					my $local_seq = $_[0];
					my $parent = $local_seq;
					my $winner = $parent;					
					my $i = 0;
					my $local_max = 0;	
					my @founders = ();
					my $n = 0;
					my $local_best_seq = $parent;
					my $local_best_score = &evaluate_sequence_2($local_best_seq);
					while ($n < $founder_size)
						{						
						push (@founders, $parent);
						++$n;
						}					
					while ($i < $generations)
						{						
						print "generation $i\: $local_best_score\n";
						my @next_founders = ();
						my %population = ();
						my %scores = ();												
						my $j = 0;							
						while ($j <= $population_size)
							{
							my $local_parent_ID_1 = int(rand($#founders));
							my $local_parent_ID_2 = int(rand($#founders));							
							my $local_parent_1 = $founders[$local_parent_ID_1];
							my $local_parent_2 = $founders[$local_parent_ID_2];
							my $genome = &recombine_parents($local_parent_1, $local_parent_2);
							$population{$j} = &mutate_sequence($genome);							
							my $test_result = &test_sequence($population{$j});							
							if ($test_result == 0)
								{
								$scores{$j} = &evaluate_sequence_2($population{$j});															
								++$j;
								}
							}
						my @sort = sort {$scores{$a} <=> $scores{$b}} keys %scores;
						@sort = reverse @sort;
						my $qq = 0;
						while ($qq < $founder_size)
							{							
							my $selected_founder = $population{$sort[$qq]};							
							if ($qq == 0)
								{								
								if ($local_best_score > $scores{$sort[$qq]})
									{
									$selected_founder = $local_best_seq;									
									}
								else	
									{
									$local_best_seq = $selected_founder;
									$local_best_score = $scores{$sort[$qq]};
									}
								}
							push (@next_founders, $selected_founder);
							++$qq;
							}								
						@founders = @next_founders;							
						++$i;
						}

					my $qq = 0;	
					open OUTX, ">FINAL_FOUNDERS.txt";						
					foreach my $a (@founders)
						{							
						print OUTX ">$qq\n$a\n\n";
						my $score = &evaluate_sequence_2($a);
						if ($score > $local_max)
							{
							$winner = $a;
							$local_max = $score;
							}
						++$qq;
						}						
					close OUTX;						
					return $winner;
					}
					
sub randomise_codons	{
						my $temp_seq = $_[0];
						my %amino_location = ();
						my $i = 0;
						while ($temp_seq =~ m/(.{3})/gmis)
							{
							my $codon = $1;							
							$amino_location{$synon{$codon}}{$i} = $codon;
							++$i;
							}							
						my @new_order = ();
						foreach my $aa (@aminos)
							{
							if (exists $amino_location{$aa})
								{
								my @start1 = keys %{$amino_location{$aa}};
								my @start2 = values %{$amino_location{$aa}};								
								my @tmp1 = @start1;
								my @tmp2 = @start2;
								my $dist_start = &calc_dist_2(\@tmp1, \@tmp2);
								my $k = 0;
								while ($k <= 100)
									{
									&fisher_yates_shuffle(\@tmp1);
									&fisher_yates_shuffle(\@tmp2);
									my $dist_current = &calc_dist_2(\@tmp1, \@tmp2);
									if ($dist_current > $dist_start)
										{										
										$dist_start = $dist_current;
										@start1 = @tmp1;
										@start2 = @tmp2;
										}
									++$k;
									}
								
								my $j = 0;
								while ($j <= $#start1)
									{
									$new_order[$start1[$j]] = $start2[$j];
									++$j;
									}								
								}
							}
						my $shuffled_sequence = join '', @new_order;
						return $shuffled_sequence;
						}

sub calc_dist_2	{
				my $array1 = $_[0];
				my $array2 = $_[1];
				my %loc = ();
				my $i = 0;
				my $cut = 40;
				my $lambda = -10;
				while ($i <= $#$array1)
					{
					$loc{$$array2[$i]}{$$array1[$i]} = 1;
					++$i;
					}
				my $dist = 0;
				my $max_dist = 0;
				foreach my $e (keys(%loc))
					{
					my @tmp = keys (%{$loc{$e}});
					foreach my $j (@tmp)
						{						
						foreach my $k (@tmp)
							{
							unless ($j == $k)
								{
								my $v = sqrt(($k - $j)**2);
								if ($v < $cut)
									{
									my $ld = exp(sqrt(($k - $j)**2)/($lambda));
									$dist += $ld;
									++$max_dist;
									}
								}
							}
						}
					}
				if ($max_dist > 0)
					{
					$dist = 1 - ($dist/$max_dist);
					}
				return $dist;
				}						
				
sub optimise_new	{				
					print "__________________________________________________\n";
					open OUT, ">Optimised_sequences\_$seq_file";
					open OUT2, ">Optimised_sequences\_$seq_file\.log";
					print OUT2 "Accession\tCAI\tEntropy\tGC_content\n";
					print "\nOptimised sequence scores\n";
					print "Accession\tCAI\tnormCAI\tEnt\tScore\tGC\tSpread\n";
					foreach my $acc (@order)
						{						
						&reset_data;
						&read_seq($seqs{$acc});						
						($MAX_CAI, $MAX, $MIN) = &compute_minmaxCAI_fast;
						$MAX_ENT = &compute_entropy_fast_2;					
						my $temp_dist = &calc_dist($seqs{$acc});
						$best = &optimisation_criteria($MAX_CAI, $MAX_ENT, $temp_dist);					
						
						%best_sequence = %temp_sequence;
						%best_codons = %temp_codons;
						
						my $gen_seq = &generate_starting_sequence($seqs{$acc});						
						$gen_seq = &evolve_sequence($gen_seq);
						
						&count_codons($gen_seq);
						&evaluate_sequence($gen_seq);
						
						$optimised{$acc} = $gen_seq;
						my $dist = &calc_dist($gen_seq);					
						my $A = $optimised{$acc} =~ tr/A/A/;
						my $T = $optimised{$acc} =~ tr/T/T/;
						my $C = $optimised{$acc} =~ tr/C/C/;
						my $G = $optimised{$acc} =~ tr/G/G/;
						my $GC = ($G + $C)/($A + $C + $T + $G);
						my $norm_cai = $MAX_CAI/$MAX;
						$dist = sprintf("%.2f", $dist);
						$GC = sprintf("%.2f", $GC);
						$MAX_CAI = sprintf("%.2f", $MAX_CAI);
						$norm_cai = sprintf("%.2f", $norm_cai);
						$MAX_ENT = sprintf("%.2f", $MAX_ENT);
						$best = sprintf("%.2f", $best);
						print "$acc\t$MAX_CAI\t$norm_cai\t$MAX_ENT\t$best\t$GC\t$dist\n";
						print OUT2 "$acc\t$MAX_CAI\t$norm_cai\t$MAX_ENT\t$GC\n";
						print OUT ">$acc\n$optimised{$acc}\n\n";
						}
					close OUT;
					close OUT2;
					}
				
sub read_seq	{
				my $s = $_[0];
				%temp_codons = ();	
				%aacount = ();			
				$codon_count = 0;
				%temp_sequence = ();		
				while ($s =~ m/(\w{3})/gmis)
					{
					my $w1 = $1;
					unless (($w1 eq 'TGA')||($w1 eq 'TAG')||($w1 eq 'TAA'))
						{
						++$temp_codons{$w1};
						$temp_sequence{$codon_count} = $w1;
						++$aacount{$synon{$w1}};
						++$codon_count;							
						}
					}
				}				
				
sub	score_all_genes	{					
					my $num = $_[0];					
					%CAI = ();
					%ENT = ();
					print "\nOriginal sequence scores\n";
					print "Accession\tCAI\tnormCAI\tEnt\tScore\tGC\tSpread\n";
					foreach my $a (@order)
						{						
						&read_seq($seqs{$a});
						my $dist = &calc_dist($seqs{$a});
						$CAI{$a} = &compute_nCAI_fast;
						$ENT{$a} = &compute_entropy_fast_2;
						($MAX_CAI, $MAX, $MIN) = &compute_minmaxCAI_fast;						
						my $local_dist = &calc_dist($seqs{$a});
						my $v = &optimisation_criteria($CAI{$a}, $ENT{$a}, $local_dist);						
						$v = sprintf("%.2f", $v);						
						
						my $A = $seqs{$a} =~ tr/A/A/;
						my $T = $seqs{$a} =~ tr/T/T/;
						my $C = $seqs{$a} =~ tr/C/C/;
						my $G = $seqs{$a} =~ tr/G/G/;
						my $gc = ($G + $C)/($A + $C + $T + $G);
						
						my $temp4 = sprintf("%.2f", $dist);
						my $temp1 = sprintf("%.2f", $CAI{$a});
						my $temp2 = sprintf("%.2f", $CAI{$a}/$MAX);
						my $temp3 = sprintf("%.2f", $ENT{$a});
						$gc = sprintf("%.2f", $gc);
						print "$a\t$temp1\t$temp2\t$temp3\t$v\t$gc\t$temp4\n";						
						}					
					}
				
sub compute_CAI_fast	{
						my $a = $_[0];												
						my $value = 0;						
						foreach my $c (keys(%{$codons{$a}}))
							{						
							$value += ($codons{$a}{$c})*($LOG{$c});								
							}
						$value = exp($value/$number_of_codons{$a});											
						return $value;				
						}
						
sub compute_minmaxCAI_fast	{
							my $value1 = 0;		
							my $value2 = 0;	
							my $value3 = 0;					
							foreach my $c (keys(%temp_codons))
								{		
								$value1 += ($temp_codons{$c})*($LOG{$c});	
								$value2 += ($temp_codons{$c})*($MAX{$c});	
								$value3 += ($temp_codons{$c})*($MIN{$c});									
								}
							$value1 = exp($value1/$codon_count);
							$value2 = exp($value2/$codon_count);
							$value3 = exp($value3/$codon_count);											
							return ($value1, $value2, $value3);							
							}
						
sub compute_nCAI_fast	{
						my $value = 0;						
						foreach my $c (keys(%temp_codons))
							{		
							$value += ($temp_codons{$c})*($LOG{$c});									
							}
						$value = exp($value/$codon_count);											
						return $value;							
						}

sub compute_entropy_fast_2	{															
							my $value = 0;
							my $count = 0;
							my $weighted_score = 0;
							my $length = 0;
							foreach my $aa (keys(%data))
								{								
								unless (($aa eq 'M')||($aa eq 'W'))
									{
									my $value2 = 0;								
									my $q = keys %{$data{$aa}};
									my $na = $aacount{$aa};
									my $count = 0;
									my $count2 = 0;
									foreach my $c (keys(%{$data{$aa}}))
										{
										if (exists $temp_codons{$c})
											{	
											my $na = $aacount{$synon{$c}};	
											$value2 -= ($temp_codons{$c}/$na)*(log($temp_codons{$c}/$na)/log(2));										
											++$count2;
											}									
										++$count;
										}
									if ($count2 > 0)
										{
										$value2 = $value2/$count;
										$weighted_score += $na*$value2;
										$length += $na;
										}									
									}
								}						
							$weighted_score /= $length;	
							return $weighted_score;				
							}
							
sub read_sequences		{						
						my $file = $_[0];
						chomp $file;
						open FILE, $file;
						my @in = <FILE>;
						close FILE;
						%seqs = ();
						my $j = join '', @in;
						$j =~ s/\r//gmis;
						$j =~ s/>/\n\n>/gmis;
						$j = $j."\n\n\n";						
						while ($j =~ m/>(.+?)\n(.+?)\n\n/gmis)
							{
							my %local = ();
							my $a = $1;							
							my $s = $2;							
							$s =~ tr/a-z/A-Z/;
							$s =~ s/[\s\n]//gmis;
							unless (($s =~ m/N/)||($s =~ m/Y/)||($s =~ m/R/)||($s =~ m/K/))
								{
								$s =~ s/TAA$//;
								$s =~ s/TGA$//;
								$s =~ s/TAG$//;
								$seqs{$a} = $s;								
								push(@order, $a);															
								my $i = 0;
								while ($s =~ m/(\w{3})/gmis)
									{
									my $w = $1;
									unless (($w eq 'TAA')||($w eq 'TGA')||($w eq 'TAG'))
										{
										++$codons{$a}{$w};										
										++$number_of_codons{$a};																			
										++$i;
										}
									}								
								}							
							}					
						}


sub set_starting_scores	{
						%START = ();
						%synon = (TTT, F, TCT, S, TAT, Y, TGT, C, TTC, F, TCC, S, TAC, Y, TGC, C, TTA, L, TCA, S, TTG, L, TCG, S, TGG, W, CTT, L, CCT, P, CAT, H, CGT, R, CTC, L, CCC, P, CAC, H, CGC, R, CTA, L, CCA, P, CAA, Q, CGA, R, CTG, L, CCG, P, CAG, Q, CGG, R, ATT, I, ACT, T, AAT, N, AGT, S, ATC, I, ACC, T, AAC, N, AGC, S, ATA, I, ACA, T, AAA, K, AGA, R, ATG, M, ACG, T, AAG, K, AGG, R, GTT, V, GCT, A, GAT, D, GGT, G, GTC, V, GCC, A, GAC, D, GGC, G, GTA, V, GCA, A, GAA, E, GGA, G, GTG, V, GCG, A, GAG, E, GGG, G);
						
						%draw_max = ();
						%lookup = ();
						%aa_groups = ();
						foreach my $c (keys(%synon))
							{
							my $a = $synon{$c};
							$aa_groups{$a}{$c} = 1;
							}
						
						if (($weights_start == 1)&&(-e $weights_file))
							{
							open FILE, $weights_file;
							while (<FILE>)
								{
								my $l = $_;
								chomp $l;
								my @t = split(/\t/, $l);
								$START{$t[0]} = $t[1];
								$LOG{$t[0]} = log($t[1]);
								}
							close FILE;
							}
							
						%MAX = ();
						%MIN = ();
						foreach my $aa (keys(%aa_groups))
							{
							my @sort = sort {$START{$a} <=> $START{$b}} keys %{$aa_groups{$aa}};
							#@sort = reverse @sort;
							my $i = 0;
							foreach my $c (@sort)
								{
								$MAX{$c} = $LOG{$sort[-1]};
								$MIN{$c} = $LOG{$sort[0]};															
								++$i;
								}
							}
							
						foreach my $aa (keys(%aa_groups))
							{
							my $i = 0;
							my $sum = 0;
							foreach my $c (keys(%{$aa_groups{$aa}}))
								{
								my $v = int(100*$START{$c});
								my $j = 0;
								while ($j < $v)
									{
									$lookup{$aa}{$i} = $c;
									++$i;
									++$j;
									}
								$sum += $v;								
								}
							$draw_max{$aa} = $sum;							
							}						
						}	
						
sub calc_dist	{
				my @local_sequence = $_[0] =~ /[ACTG]{3}/g;
				my $i = 0;
				my %loc = ();
				my $cut = 40;
				my $lambda = -10;
				foreach my $e (@local_sequence)
					{
					$loc{$e}{$i} = 1;
					++$i;
					}
				my $dist = 0;
				my $max_dist = 0;
				foreach my $e (keys(%loc))
					{
					my @tmp = keys (%{$loc{$e}});
					foreach my $j (@tmp)
						{
						foreach my $k (@tmp)
							{
							unless ($j == $k)
								{
								my $v = sqrt(($k - $j)**2);
								if ($v < $cut)
									{
									my $ld = exp(sqrt(($k - $j)**2)/($lambda));
									$dist += $ld;
									++$max_dist;
									}
								}
							}
						}
					}
				$dist = 1 - ($dist/$max_dist);
				return $dist;
				}

sub fisher_yates_shuffle	{
							my $array = $_[0];
							my $i;
							for ($i = @$array; --$i;) 
								{
								my $j = int rand ($i+1);
								next if $i == $j;
								@$array[$i,$j] = @$array[$j,$i];
								}
							}
