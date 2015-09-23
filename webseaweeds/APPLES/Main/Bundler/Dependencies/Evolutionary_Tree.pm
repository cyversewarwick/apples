### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Evolutionary_Tree Class ###

use MooseX::Declare;

class Bundler::Dependencies::Evolutionary_Tree {
	use Bundler::Dependencies::APPLES_Datatypes qw(APPLESSpeciesName Boolean NonNegativeInt);
	use Bundler::Dependencies::General_Utilities;
	use Bundler::Dependencies::Character_Matrix;
	use Bundler::Dependencies::Character_Matrix_Maker;
	use Data::Dumper;
	use constant {FALSE => 0,
		      TRUE => 1};
	
	has 'sub_trees' => (is => 'ro', isa => 'ArrayRef[Evolutionary_Tree]', required => 1);
	has 'root' => (is => 'ro', isa => APPLESSpeciesName); 

	my $GU = Bundler::Dependencies::General_Utilities->new();

	method render_text(Boolean $abbreviate_species_names_by_numbers) {
	    my @tree_species = $self->get_tree_species(TRUE);
	    my %species_display_string;
	    my $count = 0;
	    my $unknown_species_label = '(no name)';
	    if ($abbreviate_species_names_by_numbers) {
		$GU->user_info(1,"List of species abbreviations:\n");
		$unknown_species_label = 'X';
	    }
	    foreach my $species (@tree_species) {
		my $string = $species;
		if ($abbreviate_species_names_by_numbers) {
		    $count++;
		    $string = $count;		
		    $GU->user_info(1,$string.": \t".$species."\n");
		}
		$species_display_string{$species} = $string;
	    }	
	    if ($abbreviate_species_names_by_numbers) {
		$GU->user_info(1,"\nTree:\n");
	    }
	    my $character_matrix = $self->private_render_text_subroutine(\%species_display_string,$unknown_species_label);
	    $GU->user_info(1,"\n");
	    $character_matrix->render_text();
	    $GU->user_info(1,"\n");
	} # render_text #

	method contains_species (APPLESSpeciesName $species) {
		if ($self->root) {
			if ($self->root eq $species) {
				return TRUE;
			}
		}
		my $result = FALSE;
		foreach my $sub (@{$self->sub_trees}) {
			if ($sub->contains_species($species)) {
				$result = TRUE;
			}
		}
		return $result; # returns Boolean 
	} # contains_species #

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
	
	method last_common_ancestor (APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB) {
		# given a pair of species, return last common ancester in tree
		
		foreach my $sub (@{$self->sub_trees}) {
			my $contains_A = $sub->contains_species($speciesA);
			my $contains_B = $sub->contains_species($speciesB);
			if ($contains_A && $contains_B) {
				return $sub->last_common_ancestor($speciesA,$speciesB);
			}
		}
		return $self;
	} # last_common_ancestor #
	
	method get_lineage_nodes (APPLESSpeciesName $species) {
		# given a (leaf) node, returns array of all nodes in lineage from root->leaf
		my @nodes;
		if (($self->contains_species($species) eq TRUE)) {
			push (@nodes, $self->root);			
		}
		foreach my $sub (@{$self->sub_trees}){
			my @inner_nodes = $sub->get_lineage_nodes($species);
			push (@nodes, @inner_nodes);
		}
		return @nodes;
	} # get_lineage_nodes #
	
	method check_all_nodes_unique {
		# given a tree, get all nodes, check for duplicates, throw exception if any found
		my @nodes = $self->get_tree_species(FALSE);
		my @unique_list;
		my %seen;
		foreach my $element( @nodes ){
			die 'duplicated node found!' if $seen{ $element }++ ;
			push (@unique_list, $element);
		}
		$GU->user_info( 1, "\nall nodes are unique\n" );
	} # check_all_nodes_unique #

	method is_leaf {
		# given a node in a tree (a subtree), determine if leaf or not 
		# i.e. does subtree's subtree have any subtrees, if not, it's a leaf
		# return boolean
		my $lower_tree = @{$self->sub_trees}[0];
		if ($lower_tree) { 
			return FALSE;
		}
		return TRUE;
	} # is_leaf #
	
	method get_species_levels (APPLESSpeciesName $species) {
	# given a species, create an array of arrays containing groups of species at each level in
	# tree, relative to that species
		my @result;
		my $boolean = TRUE;
		if ($self->root eq $species) {
			@{$result[0]} = $self->get_tree_species($boolean);
		}	
		else {
			@{$result[0]} = ();
			foreach my $sub (@{$self->sub_trees}) {
				my $contains_species = $sub->contains_species($species);
				if ($contains_species) {
					my @inner_species = $sub->get_species_levels($species);
					push (@result, @inner_species);
				}
				else {
					my @some_on_this_level = $sub->get_tree_species($boolean);
					push (@{$result[0]}, @some_on_this_level)
				}
			}
			if ($#result == 0) {
				die $species .' not in this tree';
			}
		}	
		return @result;
	} # get_species_levels #

	method get_subtree_lineage_species (APPLESSpeciesName $species) {
		# given a query species, get all species in that subtree
		my @subtree_species;
		
		if ($self->root) {
			if ($self->root eq $species) {
				push (@subtree_species, $species);
				return @subtree_species;
			}
		}
		else {
			foreach my $sub(@{$self->sub_trees}) {
				if ($sub->contains_species($species)) {
					# get all species in subtree
					@subtree_species = $sub->get_tree_species(TRUE);
				}
			}
		}
		return @subtree_species;
	} # get_subtree_lineage_species #
	
	method get_tree_species (Boolean $leaf) {
		my @tree_species;
		if ($self->root) {
			# if check for leaf requested
			if ($leaf) {
				my $bool = $self->is_leaf;
				if ($bool) {
					push (@tree_species, $self->root);
				}
			}
			else {
				push (@tree_species, $self->root);
			}
		}
		foreach my $sub (@{$self->sub_trees}) {
			my @inner_species = $sub->get_tree_species($leaf);
			push (@tree_species, @inner_species);
		}	
		return @tree_species;
	} # get_tree_species #
	
	method get_null_subtree_species (APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB) {
		
		my @null_subtree_species;
		if ($self->root) {
			if ($self->root eq $speciesA || $self->root eq $speciesB) {
				$GU->user_info( 1, "found one of the species is the root node\n" );
				return;
			}
		}
		foreach my $sub (@{$self->sub_trees}) {
			if ( !$sub->contains_species($speciesA) && !$sub->contains_species($speciesB) ) {
				my @inner = $sub->get_tree_species(TRUE);
				push (@null_subtree_species, @inner);
			}
		}
		
		return @null_subtree_species;
	} # get_null_subtree_species #

	method prune_tree (ArrayRef[APPLESSpeciesName] $chosen_species) {
	    # produces a tree that contains only species that are in $self and in $chosen_species, preserving
	    # the topology of $self

	    my @tree_species = $self->get_tree_species(TRUE);
	    my $overlap_exists = $GU->lists_overlap(\@{$chosen_species},\@tree_species);
	    if (!$overlap_exists) {
		die 'cannot produce pruned tree as it would not have any nodes.';
	    }
	    my $result;
	    my @empty_array;
	    if ($self->is_leaf()) {
		$result = Evolutionary_Tree->new(sub_trees => \@empty_array, root => 
						 $self->root());
	    }
	    else {
		my @pruned_trees;
		foreach my $sub (@{$self->sub_trees}) {
		    my @sub_tree_leaves = $sub->get_tree_species(TRUE);
		    my $sub_tree_overlaps = $GU->lists_overlap(\@sub_tree_leaves,\@{$chosen_species});
		    if ($sub_tree_overlaps) {
			my $one_pruned_sub_tree = $sub->prune_tree($chosen_species);
			push(@pruned_trees,$one_pruned_sub_tree);
		    }
		}
		if ($#pruned_trees<0) {
		    die 'Tree pruning not implemented for inner nodes/ancient species.';
		}
		if ($#pruned_trees==0) {
		    $result = $pruned_trees[0];
		}
		else {
		    $result = Evolutionary_Tree->new(sub_trees => \@pruned_trees);
		}
	    }
	    return $result;
	} # prune_tree #

	method private_render_text_subroutine(HashRef $hash_ref,Str $unknown_species_label) {
	    my $label = $unknown_species_label;
	    if ($self->root) {
		$label = ${$hash_ref}{$self->root};
	    }
	    my $number_of_sub_trees = @{$self->sub_trees};
	    my $character_matrix;
	    my $character_matrix_maker = Character_Matrix_Maker->new();
	    if ($number_of_sub_trees == 0) {
		my $width = length($label)+2;
		$character_matrix = $character_matrix_maker->make_character_matrix($width,1);
		$character_matrix->insert_line($label,0);
	    }
	    else {
		# get sub-matrices for sub-trees
		my @sub_matrices;
		foreach my $sub_tree (@{$self->sub_trees}) {
		    my $one_character_matrix = $sub_tree->private_render_text_subroutine($hash_ref,$unknown_species_label);
		    push(@sub_matrices,$one_character_matrix);
		}                   
		# work out width and height of joint matrix
		my $total_width = 0;	    
		my $maximum_height = 0;
		foreach my $sub_matrix (@sub_matrices) {
		    $total_width = $total_width+$sub_matrix->width;
		    if ($maximum_height<$sub_matrix->height) {
			$maximum_height = $sub_matrix->height;
		    }
		}
		my $height = $maximum_height+3;
		my @numbers = ($total_width,length($label)+2);
		my $width = $GU->maximum(\@numbers);
		# initialise character matrix
		$character_matrix = $character_matrix_maker->make_character_matrix($width,$height);
		# produce first line (label)
		$character_matrix->insert_line($label,0);
		# produce second (horizontal connectors) and third (vertical bars for each sub-tree) line 
		$character_matrix->insert_line('|',1);
		my $second_line = '';
		my $third_line = '';
		for (my $tree_number=0;$tree_number<$number_of_sub_trees;$tree_number++) {
		    my $one_sub_matrix = $sub_matrices[$tree_number];
		    my $width = $one_sub_matrix->width;
		    my $half = int($width/2);
		    my $remaining_length = $width-($half+1);
		    my $first_bit_spaces = $GU->repeat_character(' ',$half);
		    my $first_bit_underscores = $GU->repeat_character('_',$half);
		    my $last_bit_spaces = $GU->repeat_character(' ',$remaining_length);
		    my $last_bit_underscores = $GU->repeat_character('_',$remaining_length);
		    $third_line = $third_line.$first_bit_spaces.'|'.$last_bit_spaces;
		    if ($tree_number==0) {
			$second_line = $second_line.$first_bit_spaces.' '.$last_bit_underscores;
		    } else {
			if ($tree_number==$number_of_sub_trees-1) {
			    $second_line = $second_line.$first_bit_underscores.' '.$last_bit_spaces;
			} else {
			    $second_line = $second_line.$first_bit_underscores.'_'.$last_bit_underscores;
			}
		    }
		}
		$character_matrix->insert_line($second_line,1);
		$character_matrix->insert_line('|',1);
		$character_matrix->insert_line($third_line,2);
		# add sub-matrices into larger matrix
		my $next_position = 0;
		for (my $tree_number=0;$tree_number<$number_of_sub_trees;$tree_number++) {
		    my $one_sub_matrix = $sub_matrices[$tree_number];
		    $character_matrix->insert_matrix($one_sub_matrix,$next_position,3);
		    $next_position = $next_position + $one_sub_matrix->width;
		}
	    }
	    return $character_matrix;
	} # private_render_text_subroutine #

} # Evolutionary_Tree #
