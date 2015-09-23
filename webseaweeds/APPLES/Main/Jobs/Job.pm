#!/usr/bin/perl -w

use MooseX::Declare;
use Data::Dumper;
=head1 Class Jobs::Job

=cut

class Jobs::Job extends Runtime::Computation {
	use JSON;

	use Runtime;
	use Serialization::Serializable;
	use Scalar::Util qw(blessed reftype);

	use Datatypes::Moose_Types;

	## Every job outputs to exactly one dataset
	has 'output_dataset' => (
		is            => "rw",
		isa           => "ResultDataset",
		default       => sub { return "Default" },
		documentation => "_ automatically filled in forms: the working dataset",
		required      => 1,
	);

=head2 This is a convenience function for deserializing a
Data Item from the result database

If the source points to a serializable class, we will
deserialize the data item.

Otherwise, we will return the item just as returned by ResultDB.

The created element of the data item will be used for
versioning: only the most recent item for this source
is returned.

 Parameters:
 $source : the source id of the data item
 $version : the version (undef for latest version)
            if the given version does not exist, the
            most recent version of the object is returned

 Returns:
 a single item or undef if this item does not exist yet

=cut

	method get_dataitem( Str $source, Str $version? ) {
		my $data =
		  jobservice->get_dataitem( $self->output_dataset, $source, $version );

		if ( !defined($data) ) {
			return undef;
		}
		if ( $data->{type} =~ m/json$/i ) {
			$data = $data->{data};
			if ( exists( $data->{SERIAL_VERSIONID} ) ) {
				$data = Serialization::Serializable::from_hash($data);
			}
		} else {
			$data = $data->{data};
		}

		return $data;
	}

=head2 Set a data item

Store an object as a data item (such that it can be retrieved as shown
above). Note that this function will store a new version. If multiple
jobs write the same data item, you might end up with many versions which
need to be merged at the end.

 Parameters:
 $item   : an object to store
 $source : (optional) a source identfier
           Default : if the item is blessed, its type will be
           used as default, otherwise we use ref ($item))
 $type   : (optional)

 Returns:
 Nothing

=cut

	method put_dataitem( Any $item, Str $source?, Str $type?) {
		if ( !defined($source) ) {
			$source = blessed($item) || ref($item);
		}

		if ( !defined($type) ) {
			$type = 'JSON';
		}

		if ( UNIVERSAL::isa( $item, 'Serialization::Serializable' ) ) {
			$item = Serialization::Serializable::to_hash($item);
		}

		jobservice->put_dataitem( $self->output_dataset, $source, $item,
			$type );
	}

=head2 Merge versions of a dataitem into one

 FIXME : we might consider parallel merging if this takes too
         long. This could be done by creating new jobs
         to merge in a tree.

 Parameters:
 $source : the source identifier
 $op     : a sub that takes two item values as its parameter, and
           returns one item.
 $type   : the type to merge into (defaults to JSON)
           note that not all data items processed must
           be of this type, this parameter only exists to make sure
           the output type is clearly defined.

           If data items of different types are stored with the
           same source ID, the merging operator $op must be
           able to distinguish between them, and merge them.

 Returns:
 Nothing

=cut

	method merge_dataitem_versions (Str $source, CodeRef $op, Str $type?) {
		my @versions =
		  @{jobservice->get_dataitem_versions( $self->output_dataset, $source )};
		my @versions_copy = @versions;

		if ( scalar @versions < 2 ) {
			return;
		}
		
		@versions = sort {$a <=> $b} @versions;

		my $data1 = $self->get_dataitem( $source, shift @versions );
		while ( scalar @versions > 0 ) {
			my $data2 = $self->get_dataitem( $source, shift @versions );

			$data1 = $op->( $data1, $data2 );
		}
		if ( !defined($type) ) {
			$type = "JSON";
		}
		$self->put_dataitem( $data1, $source, $type );
		jobservice->purge_dataitem_versions( $self->output_dataset, $source,
			\@versions_copy );
	}

=head2 Purge all but the last recent version of a data item

 Parameters:
 $source : the source identifier

 Returns:
 Nothing

=cut

	method purge_dataitem_versions ( Str $source ) {
		jobservice->purge_dataitem_versions( $self->output_dataset, $source );
	}

=head2 Override to specify the number of MPI cpus this job can use

 Returns: an Integer number of CPUs

=cut

	method ncpus () {
		return 1;
	}

=head2 Override to specify the number of threads this job can use

 Returns: an Integer number of threads

=cut

	method nthreads () {
		return 1;
	}

};
