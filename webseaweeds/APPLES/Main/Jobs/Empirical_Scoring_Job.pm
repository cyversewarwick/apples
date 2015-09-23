#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Empirical_Scoring_Job extends Jobs::Job {

	with "Jobs::Roles::Dataitem_Selection_Job";

	use Runtime;

	use Datatypes::Moose_Types;

	use Scalar::Util qw(blessed);
	use JSON;

	use Jobs::Subtasks::Empirical_Scoring;
	use Serialization::Serializable_Array;
	use Datatypes::Sequence::Set;
    use Datatypes::Motifs::PSSM_Set;

	## which overrepresentation test to use
	has 'overrep_test' => (
		is => 'rw',
		isa => 'Overrep_Test',
		default => default_value( 'Overrep_Test' ),
		required => 1,
		documentation => "X1: Which overrepresentation filter to use.|".
		"Acceptable tests: none -- return all motifs that have a good enough p-value; binomial(_R) : binomial distribution overrepresentation filter. Use with short sequences. The _R version slightly slower and should only be used for testing/verification.; poissonbinomial: slightly more accurate than the binomial filter, use for short sequences only; qvalue: test all p-values also against a q-value threshold; qvalue_binomial: use optimal number of motifs from binomial filter, but only exclude motifs if their q-value is too low."
	);

	## threshold for accepting overrepresented motifs
    has 'overrep_threshold' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 0.02; },
        required => 1,
		documentation => "X2: Overrepresentation threshold for accepting matches.",
    );

	## p-value for accepting matches
    has 'pvalue' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 0.001; },
        required => 1,
		documentation => "1. Required p-value for accepting matches.",
    );

	## The score type
	has 'score_type' => (
        is => 'rw',
        isa => 'PSSM_Scoringfunction',
        default  => default_value('PSSM_Scoringfunction'),
        required => 1,
		documentation => "2. Which score type to use.|".
		"Acceptable are: bio or mult; if mult: you can specify the pseudocount type as well: mult none, mult linear, mult sqrt, or mult bifa",
	);

	## Use repeat masked sequences
	has 'repeatmasked' => (
		is => 'rw',
		isa => 'Bool',
		default => sub { return 0; },
		required => 1,
		documentation => '3. Use repeatmasked sequences.'
	);

	## this is where we write the output to, not cached, and only used
	## in postprocess
	has 'nc__output_dataitem' => (
		is            => "rw",
		isa           => "Str",
		default       => sub { return "PSSM analysis result"; },
		documentation => "O. The output dataitem to write to.",
	);

=head2 Overloaded validation method

=cut

	method validate () {
		if ( scalar @{ $self->selected_dataitems } < 1 ) {
			die "You need to select at least one data item.";
		}
		$self->SUPER::validate();
	}

=head2 Run analysis job

 Parameters:
 None

 Returns:
 A serializable array of Genomic_Sequence objects

=cut

	method _run () {
        my @sequences = $self->selected_objects_of_kind("Sequences::Genomic_Sequence");
        my @pssms = $self->selected_objects_of_kind("Datatypes::Motifs::PSSM");

		die "You must select at least one sequence"
			if scalar @sequences == 0;

		die "You must select at least one pssm set"
			if scalar @pssms == 0;

        my $j = Jobs::Subtasks::Empirical_Scoring->new (
            overrep_test => $self->overrep_test,
            overrep_threshold => $self->overrep_threshold,
            pvalue => $self->pvalue,
            sequences => Datatypes::Sequence::Set->new (
                sequences => \@sequences,
            ),
            pssms => Datatypes::Motifs::PSSM_Set->new (
                pssms => \@pssms
            ),
			score_type => $self->score_type,
			repeatmasked => $self->repeatmasked
        );

		return $j->run();
	}

=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		$self->put_dataitem(
			$result,
			$self->nc__output_dataitem()
		);
		return $result;
	}

=head2 Expiry time (disable cache)

(all caching is handled by the created BiFa subtasks)

=cut

	method expiry_time() {
		return 0;
	}

};
