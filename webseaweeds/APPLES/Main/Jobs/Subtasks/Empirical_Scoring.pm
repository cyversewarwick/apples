#!/usr/bin/perl

=head1 Empirical scoring using the PSSM Tool

=cut

use MooseX::Declare;

class Jobs::Subtasks::Empirical_Scoring extends Jobs::Job {

	use Runtime;
	use Configuration::AppleSeeds;

	use Scalar::Util qw(blessed);
	use JSON;
	use File::Temp qw(tempfile tempdir);
	use Cwd;

	use Datatypes::Moose_Types;

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

	method _run () {
		my $pssmtool_executable = find_executable("pssmtool");

		if ( !-e $pssmtool_executable ) {
			die "PSSM Tool executable $pssmtool_executable not found.";
		}

		my %matrices_by_name = ();

        ## write matrices to file
		my @matrices = @ { $self->pssms()->pssms };
        my ($fh, $filename) = tempfile(DIR=>getcwd, UNLINK => 0);
        foreach my $matrix ( @matrices ) {
            $matrices_by_name{$matrix->name()} = $matrix;
            print $fh to_json (Serialization::Serializable::to_hash($matrix), {allow_blessed => 1}) . "\n";
        }
        close ($fh);

		my %stf = ();

		my $dir = tempdir( DIR=>getcwd, UNLINK => 0, CLEANUP => 0 );
		foreach my $seq (@{$self->sequences->sequences}) {
			my ($fh, $filename) = tempfile(DIR => $dir, UNLINK => 0);
			$stf{$seq->id()} = {
				seq => $seq,
				file => $filename,
			};

			if ($self->repeatmasked) {
				if (!defined $seq->{masked_sequence}) {
					warn ("No masked sequence available for " . $seq->id());
				}
				print $fh ($seq->{masked_sequence} || $seq->seq());
			} else {
				print $fh $seq->seq();
			}
			close $fh;
		}

		## Run the scorer

		my $pdir = get_config_key('pssmtool_histogram_directory');

		die "We need a writable histogram directory for the PSSM tool, please configure as pssmtool_histogram_directory"
			if !defined ($pdir) || !-d $pdir || !-w $pdir;

		my $st = lc($self->score_type);

		$st =~ s/[^a-z]//g;

		## All tests which are not done within PSSM_Tool need the output of all significant matches.
		my $binop = -1;
		## if we use PSSM_Tool's integrated binomial test, we enable it here.
		if ($self->overrep_test eq "binomial") {
			$binop = $self->overrep_threshold;
		} elsif ($self->overrep_test eq 'binomial_and_qvalue') {
			$binop = 1;
		}

		my $exec_str = "$pssmtool_executable -a score -f $filename --sequence $dir --profiles $pdir --scoretype $st -p "
			. $self->pvalue . " --binomial_p $binop";

		debug ($exec_str);

		print `$exec_str`;

		my @sequences_to_return = ();

		while ( my ($k, $v) = each (%stf)) {
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

			if (!defined ($hits)) {
				warn ("Unable to read result for sequence $k");
				goto NEXTONE;
			}

			if (ref ($hits) ne "ARRAY") {
				warn( "Unable to read result for sequence $k : " . ref($hits));
				goto NEXTONE;
			}

			my $seq = $v->{seq};
			my $startpos = $seq->{five_prime_pos};
			my $factor   = 1;
			if ($seq->{strand} eq 'negative') {
				$factor   = -1;
			}

			my $hitsfound = 0;
			my $presencefound = 0;
			
			my %pvals = ();

			foreach my $hit (@$hits) {
				if ($hit->{SERIAL_VERSIONID} =~ m/^Datatypes\:\:Motifs\:\:Motif_Hit/) {
					$hit = Serialization::Serializable::from_hash($hit);
					# debug ($hit->name . " : pval " . $hit->pvalue);
					if (defined ($pvals{$hit->{'name'}. "_".$hit->{accession}})) {
						push @{$pvals{$hit->{'name'}. "_".$hit->{accession}}->{hits}}, $hit;
					} else {
						$pvals{$hit->{'name'} . "_" . $hit->{accession}} = { 
							hits => [ $hit ] 
						};
					}
					++$hitsfound;
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
				}
			}

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
			
			while (my ($xk, $xv) = each (%pvals)) {
				my @pvals = map {$_->{pvalue}} @{$xv->{hits}};
				use Data::Dumper;
				if (scalar @pvals > 0 && defined $xv->{motifpresence}) {
					# debug($xv->{motifpresence}->collected_scores . " samples, " . Dumper(\@pvals));
					
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
						$testresult = $test->test($xv->{motifpresence}->collected_scores);
					} elsif ($self->overrep_test eq "poissonbinomial_R") {
						$test = Statistics::Overrepresentation_Test::PoissonBinomial->new (
							pvalues => \@pvals,
							threshold => $self->overrep_threshold
						);
						$testresult = $test->test($xv->{motifpresence}->collected_scores);
					} else {
						die "Unknown overrepresentation test " . $self->overrep_test;
					}
															
					## add the hits we accepted
					if (scalar @{$testresult->{indexes}} > 0) {
						if (!defined ($seq->{metadata})) {
							$seq->{metadata} = {};
						}
						if (!defined ($seq->{metadata}->{pssmtool})) {
							$seq->{metadata}->{pssmtool} = [  ];
						}

						$xv->{motifpresence}->sequence_name($seq->id);
						## if the PSSM Tool integrated test was used, this
						## has already happened above
						if ($self->overrep_test ne "binomial") {
							my $testpresence = $xv->{motifpresence}->cloneme();
							$testpresence->test($self->overrep_test);
							$testpresence->count(scalar @{$testresult->{indexes}});
							$testpresence->pvalue($testresult->{pvalue})
								unless $self->overrep_test =~ m/qvalue/;
							push @{$seq->{metadata}->{pssmtool}}, $testpresence;
						} else {
							push @{$seq->{metadata}->{pssmtool}}, $xv->{motifpresence};
						}
						++$presencefound;
						
						for my $ix (@{$testresult->{indexes}}) {
							my $hit = $xv->{hits}->[$ix];
							my $ug = Data::UUID->new();
							my $hit_strand = "positive";

							my $hit_start = $startpos + $factor * $hit->{five_prime_pos};
							my $hit_end =   $startpos + $factor * $hit->{three_prime_pos};

							## hit strands will be relative to the sequence strand
							if ($hit->{'strand'} eq 'positive') {
								$hit_strand = $seq->{strand};
							} else {
								$hit_strand = $seq->{strand} eq 'positive' ? 'negative' : 'positive';
							}

							my $gann = Sequences::Annotation->new(
								$ug->create_str(),
								"pssmhit",
								"PSSM Hit: " . $hit->{'name'} . " (" . $hit->{'accession'} . ") pval:" . $hit->{'pvalue'},
								$hit->{'name'},
								$hit_start,
								$hit_end,
								$hit_strand,
							);

							$gann->{hit} = $hit;
							$seq->add_annotation($gann);
						}						
					}
				}
			}
			
			push @sequences_to_return, $seq;
			debug ("Sequence " . $seq->id . ": $hitsfound hits and $presencefound overrepresented.");
			NEXTONE:
		}

		return Datatypes::Sequence::Set->new (sequences => \@sequences_to_return);
    }

}
