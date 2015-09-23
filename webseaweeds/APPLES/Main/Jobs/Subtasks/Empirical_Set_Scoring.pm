#!/usr/bin/perl

=head1 Empirical scoring using the PSSM Tool

=cut

use MooseX::Declare;

class Jobs::Subtasks::Empirical_Set_Scoring extends Jobs::Job {

	use Runtime;
	use Configuration::AppleSeeds;

	use Scalar::Util qw(blessed);
	use JSON;
	use File::Temp qw(tempfile tempdir);
	use Cwd;

	use Datatypes::Moose_Types;
    use Data::Dumper;
    
	require Serialization::Serializable_Array;
	require Datatypes::Motifs::PSSM;
	require Datatypes::Motifs::Motif_Presence;
	require Datatypes::Concepts::Transcription_Factor;
	require Datatypes::Sequence::Set;

	require Statistics::Overrepresentation_Test;
	require Statistics::Overrepresentation_Test::Binomial;
	require Statistics::Overrepresentation_Test::PoissonBinomial;
	require Statistics::Overrepresentation_Test::QValue;

	## the pssms to work with
	has 'pssms' => (
		is => "rw",
		isa => 'Datatypes::Motifs::PSSM_Set',
		required => 1,
	);

	## the sequences to work with
	has 'sequences' => (
		is => "rw",
		isa => 'Datatypes::Sequence::Set',
		required => 1,
	);

	## p-value for accepting matches
    has 'pvalue' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 0.001; },
        required => 1,
		documentation => "Required p-value for accepting matches.",
    );

	## p-value for accepting matches
    has 'overrep_threshold' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 0.02; },
        required => 1,
		documentation => "Overrepresentation test threshold for accepting matches.",
    );

	## which overrepresentation test to use
	has 'overrep_test' => (
		is => 'rw',
		isa => 'Overrep_Test',
		default => default_value( 'Overrep_Test' ),
		required => 1,
		documentation => "Which Overrepresentation test to use."
	);

	## The score type
	has 'score_type' => (
        is => 'rw',
        isa => 'PSSM_Scoringfunction',
        default  => default_value( 'PSSM_Scoringfunction' ),
        required => 1,
		documentation => "Which score type to use.|".
		"Acceptable are: bio or mult; if mult: you can specify the pseudocount type as well: mult none, mult linear, mult sqrt, or mult bifa",
	);

	## Use repeat masked sequences
	has 'repeatmasked' => (
		is => 'rw',
		isa => 'Bool',
		default => sub { return 1; },
		required => 1,
		documentation => 'Use repeatmasked sequences.'
	);

=head2 Run Scoring

 Parameters:
 None

 Returns:
 a Sequence_Set with added motif hits

=cut

	method ran () {
		my $pssmtool_executable = find_executable("pssmtool");

		if ( !-e $pssmtool_executable ) {
			die "PSSM Tool executable $pssmtool_executable not found.";
		}

		my %matrices_by_name = ();

        ## write matrices to file
        ## store matrix names in a hash - i.e. matrix by name
		my @matrices = @ { $self->pssms()->pssms };
        my ($fh, $filename) = tempfile(DIR=>getcwd, UNLINK => 0);
        foreach my $matrix ( @matrices ) {
            $matrices_by_name{$matrix->name()} = $matrix;
            print $fh to_json (Serialization::Serializable::to_hash($matrix), {allow_blessed => 1}) . "\n";
        }
        close ($fh);

        #Collect sequences & store them in files
		my %stf = ();

        #Make a temp directory
		my $dir = tempdir( DIR=>getcwd, UNLINK => 0, CLEANUP => 0 );
        #Go through all of the sequences
		foreach my $seq (@{$self->sequences->sequences}) {
            #Make a temporary file
			my ($fh, $filename) = tempfile(DIR => $dir, UNLINK => 0);
            #Put the sequence in this file
			$stf{$seq->id()} = {
				seq => $seq,
				file => $filename,
			};

			if ($self->repeatmasked) {
				if (!defined $seq->{masked_sequence}) {
					die ("No masked sequence available for " . $seq->id());
				}
				print $fh ($seq->{masked_sequence} || $seq->seq());
			} else {
				print $fh $seq->seq();
			}
			close $fh;
		}
        
        
        
		##<--Run the scoring-->##

        #Where is the pssm tool?
		my $pdir = get_config_key('pssmtool_histogram_directory');

		die "We need a writable histogram directory for the PSSM tool, please configure as pssmtool_histogram_directory"
			if !defined ($pdir) || !-d $pdir || !-w $pdir;

		my $st = lc($self->score_type);

		$st =~ s/[^a-z]//g;

        #=Do we use the integrated stuff?=
		## All tests which are not done within PSSM_Tool need the output of all significant matches.
		my $binop = -1;
		## if we use PSSM_Tool's integrated binomial test, we enable it here.
		if ($self->overrep_test eq "binomial") {
			$binop = $self->overrep_threshold;
		} elsif ($self->overrep_test eq 'binomial_and_qvalue') {
			$binop = 1;
		}

        
        
        ##Execute the PSSM tool for all sequences
		my $exec_str = "$pssmtool_executable -a score -f $filename --sequence $dir --profiles $pdir --scoretype $st -p "
			. $self->pvalue . " --binomial_p $binop";
		debug ($exec_str);
		print `$exec_str`;
        

        #Collect global pvals
        my %global_pvals = ();
        my %pvals = ();
        
        #For each sequences
		my @sequences_to_return = ();
        
		while ( my ($k, $v) = each (%stf)) {
            
            #Get the hits for this sequence and all matrices
			my $hits;
			{
				open FIN, "<", "$v->{file}.motifs.json"
					or do {
						warn ("PSSM_Tool didn't output anything for sequence $k ( $v->{file} )");
						goto NEXTONE;
					};
				local $/ = undef;
				my $txt = <FIN>;
				$hits = from_json( $txt );
				close FIN;
			}

            #If no hits
			if (!defined ($hits)) {
				warn ("Unable to read result for sequence $k");
				goto NEXTONE;
			}
			if (ref ($hits) ne "ARRAY") {
				warn( "Unable to read result for sequence $k : " . ref($hits));
				goto NEXTONE;
			}
            
            
            #Get the sequence, its start pos, and strand
			my $seq = $v->{seq};
			my $startpos = $seq->{five_prime_pos};
			my $factor   = 1;
			if ($seq->{strand} eq 'negative') {
				$factor   = -1;
			}

            #Did we find the hits?
			my $hitsfound = 0;
			my $presencefound = 0;

			foreach my $hit (@$hits) {
                
                #If we've got a motif hit
				if ($hit->{SERIAL_VERSIONID} =~ m/^Datatypes\:\:Motifs\:\:Motif_Hit/) {
					$hit = Serialization::Serializable::from_hash($hit);
					# debug ($hit->name . " : pval " . $hit->pvalue);
                    
                    #Add to our pvals hash, which is matrixname_accession format
					if (defined ($pvals{$hit->{'name'}. "_".$hit->{accession}})) {
						push @{$pvals{$hit->{'name'}. "_".$hit->{accession}}->{hits}}, $hit;
					} else {
						$pvals{$hit->{'name'} . "_" . $hit->{accession}} = {
							hits => [ $hit ]
						};
                    }
					
                    #Add it to our global pvals
                    if (defined ($global_pvals{$k}->{$hit->{'name'}. "_".$hit->{accession}})) {
						push @{$global_pvals{$k}->{$hit->{'name'}. "_".$hit->{accession}}->{hits}}, $hit;
					} else {
						$global_pvals{$k}->{$hit->{'name'} . "_" . $hit->{accession}} = {
							hits => [ $hit ]
						};
                    }
                
                    #This is how many hits we found
					++$hitsfound;
                    
                #Or if we've got a motif presence
				} elsif ($hit->{SERIAL_VERSIONID} =~ m/^Datatypes\:\:Motifs\:\:Motif_Presence/) {
                    
                    
					my $mp = Serialization::Serializable::from_hash($hit);
					if (defined ($pvals{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}})) {
						$pvals{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}}->{motifpresence} = $mp;
					} else {
						$pvals{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}} = {
							motifpresence => $mp,
							hits => [],
						};
                    }
                    
                    if (defined ($global_pvals{$k}->{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}})) {
						$global_pvals{$k}->{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}}->{motifpresence} = $mp;
					} else {
						$global_pvals{$k}->{$mp->{'pssm_name'}. "_".$mp->{pssm_accession}} = {
							motifpresence => $mp,
							hits => [],
						};
					}
                    
				}
			}
            
            #Ignore all this
			if ($self->overrep_test =~ m/qvalue/) {
				my @ppvals = ();
				my @ppval_hits = ();
				my %motifpresences = ();
				
				my $i = 0;				
				while (my ($xk, $xv) = each (%pvals)) {
					for (my $i = 0; $i < scalar @{$xv->{hits}}; ++$i) {
						push @ppvals, $xv->{hits}->[$i]->{pvalue};
						push @ppval_hits, $xv->{hits}->[$i++];
					}
					$motifpresences{$xk} = $xv->{motifpresence};
				}
	
				my $qtest = Statistics::Overrepresentation_Test::QValue->new (
					pvalues => \@ppvals,
					threshold => $self->overrep_threshold
				);
				
				my $testresult = $qtest->test();
				
				%pvals = ();
				
				for my $ix (@{$testresult->{indexes}}) {
					my $hit = $ppval_hits[$ix];
					$hit->qvalue($testresult->{all_qvalues}->[$ix]);
					
					if (defined ($pvals{$hit->{'name'}. "_".$hit->{accession}})) {
						push @{$pvals{$hit->{'name'}. "_".$hit->{accession}}->{hits}}, $hit;
					} else {
						$pvals{$hit->{'name'} . "_" . $hit->{accession}} = { 
							hits => [ $hit ] 
						};
					}
					print STDERR $hit->{'name'} . " : " . $hit->pvalue . " / " . $hit->qvalue . "\n";
				}
				
				while (my ($xk, $xv) = each (%pvals)) {
					my $pv = 0;
					my $qv = 0;
					for (my $i = 0; $i < scalar @{$xv->{hits}}; ++$i) {
						if ($xv->{hits}->[$i]->{pvalue} > $pv) {
							$pv = $xv->{hits}->[$i]->{pvalue};
						}
						if ($xv->{hits}->[$i]->{qvalue} > $qv) {
							$qv = $xv->{hits}->[$i]->{qvalue};
						}
					}
					$xv->{motifpresence} = $motifpresences{$xk};
					$xv->{motifpresence}->pvalue($pv);
					$xv->{motifpresence}->{qvalue} = $qv;
				}
			}
            #//Ignore
			
                        

			NEXTONE:
		}

        my %all_pvals_fixed = ();
        
        #Now for each motif...
        #Tomorrow, you need to put all the pvals together, add a sequence field, and fix the PSSM motifpresence
        foreach my $gene (keys %global_pvals)
        {
            my %pvals = %{$global_pvals{$gene}};
            
            while (my ($m_id, $vals) = each (%pvals))
            {
                my $new_count = $vals->{"motifpresence"}->{"count"};
                my $new_sequence_name = $vals->{"motifpresence"}->{"sequence_name"};
                my $new_collected_scores = $vals->{"motifpresence"}->{"collected_scores"};

                my $pssm_name = $vals->{"motifpresence"}->{"pssm_name"};
                my $pssm_accession = $vals->{"motifpresence"}->{"pssm_accession"};
                my $serial_id = $vals->{"motifpresence"}->{"SERIAL_VERSIONID"};

                #Go through hits
                foreach my $hit (@{$vals->{hits}})
                {
                    $hit->{"gene"} = $gene;
                    push(@{$all_pvals_fixed{$m_id}->{"hits"}}, $hit);
                }
                
                #Sort out PSSM 
                if(defined($all_pvals_fixed{$m_id}->{"motifpresence"}))
                {
                    $all_pvals_fixed{$m_id}->{"motifpresence"}->{"count"} = $all_pvals_fixed{$m_id}->{"motifpresence"}->{"count"} + $new_count;
                    $all_pvals_fixed{$m_id}->{"motifpresence"}->{"collected_scores"} = $all_pvals_fixed{$m_id}->{"motifpresence"}->{"count"} + $new_collected_scores;
                }
                else
                {
                    $all_pvals_fixed{$m_id}->{"motifpresence"}->{"count"} = $new_count;
                    $all_pvals_fixed{$m_id}->{"motifpresence"}->{"collected_scores"} = $new_collected_scores;
                }
                
                push(@{$all_pvals_fixed{$m_id}->{"motifpresence"}->{"sequence_names"}}, $new_sequence_name);
                push(@{$all_pvals_fixed{$m_id}->{"motifpresence"}->{"sequence_accessions"}}, $gene);
                
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"pssm_name"} = $pssm_name;
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"pssm_accession"} = $pssm_accession;
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"SERIAL_VERSIONID"} = $serial_id;
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"links"} = {};
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"test"} = "none";
                $all_pvals_fixed{$m_id}->{"motifpresence"}->{"pvalue"} = 1;
            }
        }

        my %pvals = %all_pvals_fixed;
        
        my @motifs_to_return = ();
        
            while (my ($xk, $xv) = each (%pvals)) {
                
                #Empty motif collection
                my %cur_motif = ();
                
                my @pvals = map {$_->{pvalue}} @{$xv->{hits}};
                
                #<==Sort pvals and take at most 100==>#
                @pvals = sort {$a <=> $b} @pvals;
                my $limit=99;
                if(scalar(@pvals) <= $limit)
                {
                    $limit = scalar(@pvals)-1;
                }
                @pvals = @pvals[0..$limit];
                
                if (scalar @pvals > 0 && defined $xv->{motifpresence}) {
                    # debug($xv->{motifpresence}->collected_scores . " samples, " . Dumper(\@pvals));
                    #print "\n" . $xv->{motifpresence}->{collected_scores} . " resulted in " . $xv->{motifpresence}->{count} .  " samples, " . Dumper(\@pvals);
                    $cur_motif{"pssm_name"} = $xv->{motifpresence}->{pssm_name};
                    $cur_motif{"count"} = $xv->{motifpresence}->{count};
                    $cur_motif{"sequence_accessions"} = $xv->{motifpresence}->{sequence_accessions};
                    $cur_motif{"collected_scores"} = $xv->{motifpresence}->{collected_scores};
                    $cur_motif{"pssm_accession"} = $xv->{motifpresence}->{pssm_accession};
                    
                    
                    my $test;
                    my $testresult;
                    

                    
                    if ($self->overrep_test eq "none"
                        || $self->overrep_test eq "binomial"
                        ) {
                            $test = Statistics::Overrepresentation_Test->new (
                            pvalues => \@pvals,
                            threshold => $self->overrep_threshold
                            );
                            $testresult = $test->test();
                        } elsif ($self->overrep_test =~ m/qvalue/) {
                            $test = Statistics::Overrepresentation_Test->new (
                            pvalues => \@pvals,
                            threshold => 1.0
                            );
                            $testresult = $test->test();
                        } elsif ($self->overrep_test eq "binomial_R") {
                            $test = Statistics::Overrepresentation_Test::Binomial->new (
                            pvalues => \@pvals,
                            threshold => $self->overrep_threshold
                            );
                            $testresult = $test->test($xv->{motifpresence}->{collected_scores});
                        } elsif ($self->overrep_test eq "poissonbinomial_R") {
                            $test = Statistics::Overrepresentation_Test::PoissonBinomial->new (
                            pvalues => \@pvals,
                            threshold => $self->overrep_threshold
                            );
                            $testresult = $test->test($xv->{motifpresence}->{collected_scores});
                        } else {
                            die "Unknown overrepresentation test " . $self->overrep_test;
                        }
                
                    print "\nXK: $xk" . "\t" . $testresult->{pvalue};
                    ## add the hits we accepted
                    
                        if (scalar @{$testresult->{indexes}} > 0)
                        {
                            $cur_motif{"pvalue"} = $testresult->{pvalue};
                            $cur_motif{"result_count"} = scalar @{$testresult->{indexes}};
                            
                            push(@motifs_to_return, \%cur_motif);
                        }
                    }
                }
        return \@motifs_to_return;
    }

}
