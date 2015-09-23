#!/usr/bin/perl

use MooseX::Declare;

=head1 Gene_Set Class

A gene set is a collection of genes. There might be duplicates, which can be 
eliminated by calling combine_genes.

=cut

class Datatypes::Concepts::Gene_Set extends Serialization::Serializable {
	use Runtime;
	use Scalar::Util qw(blessed);
	
	use Serialization::Serializable;
	use Datatypes::Sequence::Search;
	
	require Datatypes::Concepts::Gene;
	
	## list of genes
	has 'genes' => (
		is => 'rw',
		isa => 'ArrayRef[Datatypes::Concepts::Gene]',
	);
	
=head2 Construct a new gene list

 Parameters:
 $class : Datatypes::Concepts::Gene_Set

 The remaining parameters can be any combination of
  Strings which get resolved using Datatypes::Sequence::Search
  Datatypes::Concepts::Gene which get added
  Arrayrefs

 Returns:
 a new Datatypes::Concepts::Gene_Set

=cut
	sub new_gene_set {
		my $class = shift;
		my @ids = @_;
		my %genes = ();
		
		while (0 < scalar @ids) {
			my $iid = shift @ids;
			next unless defined $iid;
			
			if ( ref ($iid) eq "ARRAY" ) {
				## expand subarrays
				push @ids, $_ foreach @$iid;
			} elsif (Serialization::Serializable::is_hash($iid) && blessed ($iid) eq "Datatypes::Concepts::Gene") {
				my $uid = $iid->accession;
				
				if (!defined $genes{$uid}) {
					 $genes{$uid} =  [ ];
				}
				
				push @{ $genes{$uid} }, $iid;
			} else {
				my $search = Datatypes::Sequence::Search->new;
				my $locs = $search->parse_search ( "$iid" );
				
				foreach my $loc (@{$locs}) {
					my $uid = $loc->accession;
					if (!defined $genes{$uid}) {
						 $genes{$uid} =  [ ];
					}

					push @{ $genes{$uid} }, Datatypes::Concepts::Gene->new (
						db => $loc->db,
						location => $loc->location,
					);					
				}
			}
		}
		
		my @geneset = ();
		
		while (my ($k, $v) = each (%genes)) {
			my $g0 = shift @$v;
			
			while (0 < scalar @$v) {
				my $g1 = shift @$v;
				$g0->merge($g1);
			}
			push @geneset, $g0;
		}
		
		return Datatypes::Concepts::Gene_Set->new (
			genes => \@geneset,
		);
	}
	
}
