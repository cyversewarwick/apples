#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Extract_Conserved_Blocks extends Jobs::Job {

	with "Jobs::Roles::Dataitem_Selection_Job";

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	use Jobs::Subtasks::Mark_Conservation;

	use Datatypes::Sequence::Set;
	use Datatypes::Results::Alignment_Plot;

	has 'pval_threshold' => (
		is => 'rw',
		isa => 'Num',
		default => 1e-3,
		documentation => "1. Pvalue threshold"
	);
	
	has 'masks' => (
		is => 'rw',
		isa => 'Str',
		default => "",
		documentation => "2. Masking values"
	);
	
	has 'nc__resultdataitem' => (
		is => 'rw',
		isa => 'Str',
		default => "Conserved Elements",
		documentation => "3. Dataitem to store conserved elements"
	);

=head2 Overloaded validation method

=cut

	method validate () {
		if ( scalar @{ $self->selected_dataitems } < 1 ) {
			die "You need to select at least one data item.";
		}
		
		my $hr = from_json ($self->masks);
		if (!Serialization::Serializable::is_hash ($hr)) {
			die "Invalid masks string.";
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
		
	}

=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		return $result;
	}
};
