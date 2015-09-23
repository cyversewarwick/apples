#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Sequence_Search extends Jobs::Job {

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	use Datatypes::Sequence::Search;
	use Datatypes::Sequence::Location;
	use Jobs::Subtasks::Sequence_Retrieval;

	## used by postprocess to specify the
	## place where the sequences go
	## this parameter is not cached since
	## it is only used in postprocess
	has 'nc__sequenceset' => (
							   is            => 'rw',
							   isa           => 'Str',
							   documentation => 'Sequence Set Name',
							   default       => sub { return "Sequences"; },
	);

	has 'search' => (
					  is            => 'rw',
					  isa           => 'Sequence_Search',
					  documentation => 'Sequence Query',
					  default       => sub { return "" },
	);

=head2 Run retrieval job

 Parameters: 
 None
 
 Returns: 
 A serializable array of Genomic_Sequence objects
 
=cut

	method _run () {
		## get the sequence, store it in the cache.

		my $s         = Datatypes::Sequence::Search->new;
		my $locations = $s->parse_search( $self->search );
		
		my $rj =
		  Jobs::Subtasks::Sequence_Retrieval->new(
		  							'output_dataset' => $self->output_dataset,
									'nc__sequenceset' => $self->nc__sequenceset,
									location => $locations,
									 );
		my $result = $rj->run;
		return $result;
	}
};
