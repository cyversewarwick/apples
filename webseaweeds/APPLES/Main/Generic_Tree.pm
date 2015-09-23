### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Generic_Tree Class ###

use MooseX::Declare;

class Generic_Tree {

	method tree_depth {
		my $depth = 0;
		foreach my $sub (@{$self->sub_trees}) {	
			my $int_depth = $sub->tree_depth;
			if ($int_depth > $depth) {
				$depth = $int_depth;
			}
		}
		$depth++;
		return $depth;
	} # tree_depth #

	method get_tree_species {
		my @tree_species;
		if ($self->root) {
			push (@tree_species, $self->root);
		}
		foreach my $sub (@{$self->sub_trees}) {
			my @inner_species = $sub->get_tree_species;
			push (@tree_species, @inner_species);
		}
		
		return @tree_species;
	} # get_tree_species #

} # Generic_Tree #

