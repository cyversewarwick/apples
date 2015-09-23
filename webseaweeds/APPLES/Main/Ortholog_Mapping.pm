### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Ortholog_Mapping class ###
# Ortholog_Mapping object holds info on ortholog assignments between two genomes

use MooseX::Declare;

class Ortholog_Mapping {
	
	has 'genome_a' => (is => 'rw', isa => 'Genome_Sequence_Database_Parameters'); 
	has 'genome_b' => (is => 'rw', isa => 'Genome_Sequence_Database_Parameters'); 
	has 'a_to_b' => (is => 'rw', isa => 'HashRef');
	has 'b_to_a' => (is => 'rw', isa => 'HashRef');

	method get_orthologous_ID_for_source_id (Str $query_id) {
	  my %a_to_b = %{$self->a_to_b}; 
	  my $ortholog = $a_to_b{$query_id};
	  return $ortholog;
	} # get_orthologous_ID_for_source_id #

} # Ortholog_Mapping #


class Syntenic_Orthology_Mapping extends Ortholog_Mapping {

# this class is not implemented yet, just illustrating an idea

} # Syntenic_Orthology_Mapping #
