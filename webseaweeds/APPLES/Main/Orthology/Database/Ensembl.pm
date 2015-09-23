### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

=head1 Interface to Ensembl's Compara Orthology database

=cut

package Orthology::Database::Ensembl;

use strict;

require Sequences::Database::Sequence::Ensembl;
require Orthology::Database;
our @ISA = qw(Sequences::Database::Sequence::Ensembl Orthology::Database);

use Bio::EnsEMBL::Registry;
use List::Util qw(min max);
use Data::Dumper;

use Runtime;
use Carp;

use Datatypes::Concepts::Gene;

=head2 Constructor

 Parameters:
 $class : Orthology_Database_Ensembl
 $registry_location : a registry location

=cut

sub new {
	my $class = shift;
	my $registry_location = shift || "ensembl";

	my $self = Sequences::Database::Sequence::Ensembl->new($registry_location);

	bless $self, $class;
	return $self;
}

=head2 Given an ID, return orthologous ids in this database, restricting 
       to a list of species 

 Parameters:
 $self : an Orthology Database
 $id : a Datatypes::Concepts::Gene
 $list_of_species : an ARRAYREF for a list of species to match again
 
 Returns: 
 a list of orthologous id genes 
=cut

sub _fetch_orthologous_ids {
	my ( $self, $gene, $list_of_species ) = @_;

	my $id       = $gene->accession;
	my $registry = $self->registry;

	info("Finding orthologous IDs for $id");

	my $targetspecieshomologues = [];

	## find homologs in the list of species.
	## Get the compara member adaptor
	my $compara_adaptor_member =
	  $registry->get_adaptor( 'Multi', 'compara', 'Member' );

	## Get the compara homology adaptor
	my $compara_adaptor_homology =
	  $registry->get_adaptor( 'Multi', 'compara', 'homology' );

	my $qy_member =
	  $compara_adaptor_member->fetch_by_source_stable_id( "ENSEMBLGENE", $id );

	if ( !$qy_member ) {
		info(   "Skipping orthologue search for " 
			  . $id
			  . " because I could not find it in the database.\n" );
		next;    # exits the method before it crashes with 'undefined' error
	}

	info(   "Gene ENSEMBL ID: "
		  . $qy_member->stable_id
		  . ", gene common name: "
		  . $qy_member->display_label
		  . ". \n" );

	my $homologies =
	  $compara_adaptor_homology->fetch_all_by_Member_method_link_type(
														$qy_member,
														"ENSEMBL_ORTHOLOGUES" );
	foreach my $homology (@$homologies) {
		if ( $homology->description ne 'between_species_paralog' ) {
			foreach
			  my $member_attribute ( @{ $homology->get_all_Member_Attribute } )
			{
				my ( $member, $attribute ) = @{$member_attribute};
				next if ( $member->stable_id eq $qy_member->stable_id );
				
				eval {
					my ( $species, $object_type, $db_type ) =
					  $self->find_species_and_object_type( $member->stable_id );

					info( "Orthologous gene is found in $species: "
						  . $member->stable_id );

					my $species_match = 0;
					foreach my $sp ( @{$list_of_species} ) {
						$sp =~ s/\s/_/g;
						if ( $species =~ m/$sp/i ) {
							$species_match = 1;
							last;
						}
					}
					next unless $species_match;

					push @$targetspecieshomologues,
					  Datatypes::Concepts::Gene->new_gene(
							get_sequence_database_info( $self->{registry_location} )
							  ->{key},
							$member->stable_id
					  );
					info(   "Orthologous gene is added for $species: "
						  . $member->stable_id . " ( "
						  . $member->display_label
						  . " ). \n" );					
				};
				
				if ($@) {
					warn ( "Ignoring " . $member->stable_id . " because of an error: " . $@ );
				}
			}
		}
	}

	return $targetspecieshomologues;
}

1;
