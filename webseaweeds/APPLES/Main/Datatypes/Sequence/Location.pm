#!/usr/bin/perl

use MooseX::Declare;

=head1 Sequence Location Class

Stores sequence locations and databases.

=cut

class Datatypes::Sequence::Location extends Links::Node {

	with "Datatypes::Roles::Locatable";

	use Runtime;

	use Scalar::Util qw(blessed);
	use Datatypes::Moose_Types;

	require Sequences::Database::Relative_Location;

	## The ID of the sequence database to use for retrieving this sequence
	has 'db' => (
		is            => "rw",
		isa           => "Sequence_Database_ID",
		required      => 1,
		documentation => "Sequence database",
		default       => default_value("Sequence_Database_ID"),
		trigger       => sub {
			my $self = shift;
			$self->{unique_id} = $self->unique_id;
		},
	);

	## The relative location to use
	has 'location' => (
		is            => "rw",
		isa           => "Sequences::Database::Relative_Location",
		documentation => "Sequence location",
		default => default_value("Sequences::Database::Relative_Location"),
		trigger => sub {
			my $self = shift;
			$self->{unique_id} = $self->unique_id;
			$self->url ($self->get_meta_data->{url} || "");
		},
		required => 1
	);

=head2 Wrapper to get accession-style ids

 Parameters:
 $self : a self object

 Returns:
 a string with an ensembl/genbank accession

=cut

	sub accession {
		my $self = shift;
		return $self->location->identifier;
	}

=head2 Create a location from db id and accession (static)

 Parameters:
 $db : the db id
 $acc : the accession

 Returns:
 a new Datatypes::Sequence::Location

=cut

	sub new_simple_location {
		my ( $class, $db, $acc ) = @_;
		return
		  __PACKAGE__->new(
							db => $db,
							location =>
							  Sequences::Database::Relative_Location->new(
															 identifier => $acc,
							  )
		  );
	}

=head2 Get the sequence for this location

Try validating this sequence location. Throws an exception when things go wrong.

=cut

	method validate () {
		$self->SUPER::validate();
		$self->location->validate();
		
		if (lc ($self->db) ne "none") {
			get_sequence_database( $self->db )
			  ->validate_location( $self->location );
		}
	}

=head2 Unique ID override

=cut

	method unique_id () {
		return blessed($self) . " : $self->{db}" . $self->location->unique_id;
	}

=head2 Get the sequence for this location

 Parameters:
 $self : $self object

 Returns:
 a Genomic_Sequence corresponding to this location

=cut

	method get_sequence () {
		return get_sequence_database( $self->db )
		  ->get_sequence_by_location( $self->location );
	}

=head2 Get the metadata for this location

 Parameters:
 $self : $self object

 Returns:
 a metadata hash (see Sequences::Database::Metadata)
 corresponding to this location

 There might not be any metadata available, in which
 case we will return {}

=cut

	method get_meta_data () {
		
		if (lc ($self->db) ne "none") {
			my $sdb = get_sequence_database( $self->db );

			if ( UNIVERSAL::can( $sdb, 'get_meta_data' ) ) {
				return $sdb->get_meta_data( $self->accession );
			}
		}
		
		return {};
	}

=head2 Get URL via metadata

=cut

	## override this to change the color in derived classes
	method color () {
		return "#cfc";
	}

=head2 Create a dot label

=cut

	method dot_style () {
		return "label=\"" . $self->accession . "\" shape=\"box\"";
	}

}

1;
