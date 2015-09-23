#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Subtasks::Sequence_Retrieval extends Jobs::Job {

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;
	
	use Serialization::Serializable_Array;
	use Datatypes::Sequence::Set;
	require Datatypes::Sequence::Location;

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

	## a set of sequence locations
	has 'location' => (
		is            => 'rw',
		isa           => 'ArrayRef',
		documentation => "Sequence Locations",
		default       => sub { return [] },
	);

=head2 This method will be run by the cache after the result has 
been retrieved

 Parameters:
 $result : the computation result (cached or just computed)

 Returns:
 a postprocessed version of the result

=cut

	method postprocess ( Serialization::Serializable $result ) {
		my @arr = ();
		foreach my $gs ( @{ $result->data } ) {
			if ( ref($gs) eq "ARRAY" ) {
				push @arr, $_ foreach @$gs;
			} else {
				push @arr, $gs;
			}
		}

		$self->put_dataitem(
			Datatypes::Sequence::Set->new( sequences => \@arr ),
			$self->nc__sequenceset );

	  	$self->merge_dataitem_versions( $self->nc__sequenceset,
	  	 	\&Datatypes::Sequence::Set::merged_sets
		);

		return $result;
	}

=head2 Retrieve the sequences

 Parameters: 
 None
 
 Returns: 
 A serializable array of Genomic_Sequence objects
 
=cut

	method _run () {
		## get the sequence, store it in the cache.
		my @seqs = ();
		foreach my $seq ( @{ $self->location } ) {
			my $gs = $seq->get_sequence;
			if ( ref($gs) eq "ARRAY" ) {
				foreach my $ggs (@$gs) {
					push @seqs, $ggs;
				}
			} else {
				push @seqs, $gs;
			}
		}
		my $arr = Serialization::Serializable_Array->new( \@seqs );
		return $arr;
	}
};
