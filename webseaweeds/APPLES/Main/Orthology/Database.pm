
=head1 Class Orthology::Database

A generic point-to-multipoint orthology database interface

=cut

package Orthology::Database;

use strict;

use Runtime;

use Scalar::Util qw(blessed);

use Serialization::Serializable;
use Serialization::Serializable_Array;

=head2 Given an ID, return orthologous ids in this database, restricting
       to a list of species

 Parameters:
 $self : an Orthology Database
 $id : a Datatypes::Concepts::Gene
 $species : an ARRAYREF for a list of species to match against

 Returns:
 a list of orthologous id records
=cut

sub fetch_orthologous_ids {
	my $self    = shift;
	my $id      = shift;
	my $species = shift;

	my $parameters = bless {
							 type    => 'ORTHOLOGY_INFO_' . blessed($self),
							 db      => $id->db,
							 acc      => $id->accession,
							 species => $species,
	  },
	  "Serialization::Serializable";

	my $result = cache->cache_get($parameters);

	if ( !defined($result) ) {
		$result = Serialization::Serializable_Array->new(
			@{ $self->_fetch_orthologous_ids( $id, $species ) }
		);
		cache->cache_put( $parameters, $result, 'never' );
	} else {
		debug ("Orthologues for " . $id->accession . " found in cache.");
	}

	return $result->data;
}

1;
