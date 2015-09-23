### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Partial_Threshold_Matrix Class ###
# stores thresholds for some pairs of species covered in an evolutionary tree
# the evolutionary tree can be used to infer thresholds for other pairs of species
# thresholds are invariant under reversion of pairs

use MooseX::Declare;

class Partial_Threshold_Matrix {
    use APPLES_Datatypes qw(APPLESSpeciesName Boolean);
    use Evolutionary_Tree;
    use constant {FALSE => 0,
		  TRUE => 1};
    use Data::Dumper;
    use General_Utilities;

    has 'evolutionary_tree' => (is => 'ro', isa => 'Evolutionary_Tree', required => 1);
    has 'matrix' => (is => 'rw', isa => 'ArrayRef');
    has 'completeness_check' => (is => 'rw', isa => Boolean, required => 1, default => FALSE);
    has 'species_to_numbers' => (is => 'rw', isa => 'HashRef'); # private attribute

    my $GU = General_Utilities->new();

    method species_number(APPLESSpeciesName $species) {
	my $number = ${$self->species_to_numbers}{$species};
	return $number;
    } # species_number #
	
    method is_complete() {
	# check if thresholds for all pairs can be derived
	my @species = $self->evolutionary_tree->get_tree_species(TRUE);
	foreach my $species1 (@species) {
	    foreach my $species2 (@species) {
		unless ($species1 eq $species2) { 
		    my $present = $self->is_set($species1, $species2);
		    if (!$present) {
			$GU->user_info( 1, 'threshold for '.$species1.' vs '.$species2." NOT SET\n" );
			$self->completeness_check(FALSE);
			return FALSE;
		    }
		}
	    }
	}
	$self->completeness_check(TRUE);
	return TRUE;
	# sets completeness_check to TRUE only if all pair thresholds can be derived
    } # is_complete #
 	
    method consistency() {
	die 'method not implemented yet!';
 	# check that tree values set are consistent
    } # consistency #
 	
    method fill_species_to_numbers(ArrayRef $species_list_ref) {		
	my @species = sort { $a cmp $b } @{$species_list_ref};
	my %species_to_numbers = ();
	for (my $i=0; $i < @species; $i++) {
	    # assign a number to each species, in a hash
	    $GU->user_info( 3, $species[$i]."\t".$i."\n" );
	    $species_to_numbers{$species[$i]} = $i;
	    
	}
	$self->species_to_numbers(\%species_to_numbers); 
    } # fill_species_to_numbers #
 	
    method initialise_matrix() {
	my @species = $self->evolutionary_tree->get_tree_species(TRUE);	   
	# turn species to numbers
	my $species_hash_ref = $self->fill_species_to_numbers (\@species);	
	my @matrix;
	$self->matrix ( \@matrix );
    } # initialise_matrix #
 	
    method set_matrix_value_for_pair(APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB, Int $min, Int $max) {
	my $A_exists = $self->evolutionary_tree->contains_species($speciesA);
	my $B_exists = $self->evolutionary_tree->contains_species($speciesB);
	if ($A_exists && $B_exists) {
	    # sort $A and $B alphabetically
	    my @variables = ($speciesA, $speciesB);
	    my @sorted = sort { $a cmp $b } @variables;
	    $speciesA = $sorted[0];
	    $speciesB = $sorted[1];
	    # set or update thresholds for a pair of species
	    # look up key values
	    my $A = $self->species_number($speciesA);
	    my $B = $self->species_number($speciesB);
	    
	    # set matrix values
	    $self->matrix->[$A][$B][0] = $min;
	    $self->matrix->[$A][$B][1] = $max;
	}
    } # set_matrix_value_for_pair #
 	
    method is_set(APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB) {
	# checks if a pair threshold is set and returns Boolean
	# sort $A and $B alphabetically
	my @variables = ($speciesA, $speciesB);
	my @sorted = sort { $a cmp $b } @variables;
	
	$speciesA = $sorted[0];
	$speciesB = $sorted[1];

	my $A = $self->species_number($speciesA);
	my $B = $self->species_number($speciesB);

	if ( $self->matrix->[$A][$B][0] ) {
	    return TRUE;
	}
	return FALSE;
    } # is_set #
 	
 	method get_thresholds(APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB) {
	  #if ($speciesA eq 'vitis_vinifera') {
	   # $speciesA = 'grape';
	  #}
	  #if ($speciesB eq 'vitis_vinifera') {
	   # $speciesB = 'grape';
	  ##debug quick fix until find out apporporiate aliases
 		# sort $A and $B alphabetically
 		my @variables = ($speciesA, $speciesB);
		my @sorted = sort { $a cmp $b } @variables;
		$speciesA = $sorted[0];
		$speciesB = $sorted[1];
 		my $A = $self->species_number($speciesA);
 		my $B = $self->species_number($speciesB);
 		
 		my @result;
 		if ( $self->matrix->[$A][$B][0] && $self->matrix->[$A][$B][1]) {
 			push (@result, $self->matrix->[$A][$B][0], $self->matrix->[$A][$B][1]);
 			return @result;
 		}
 		return;
 	} # get_thresholds #
 	
 	method autofill() {
	    # complete as many unset thresholds in tree as possible, from existing values
	    my @unsorted_species = $self->evolutionary_tree->get_tree_species(TRUE);
	    my @species = sort { $a cmp $b } @unsorted_species; 
 		
	    foreach my $species1 (@species) {
		foreach my $species2 (@species) {
		    $GU->user_info( 3, $species1." vs ".$species2."\n" );
		    unless ($species1 eq $species2) {
			$self->infer_pair_threshold($species1, $species2);
		    }
		}
	    }
	    
	    $GU->user_info( 3, Dumper ($self->matrix) );
	    return;
 	} # autofill #
 	
 	method infer_pair_threshold(APPLESSpeciesName $speciesA, APPLESSpeciesName $speciesB) {
 		# if completeness_check eq TRUE {
 		
 		#}
 		
 		$GU->user_info( 3, Dumper ($self->matrix) );
 		
 		# sort $A and $B alphabetically
 		my @variables = ($speciesA, $speciesB);
		my @sorted = sort { $a cmp $b } @variables;
		$speciesA = $sorted[0];
		$speciesB = $sorted[1];
 		if ( $self->is_set($speciesA, $speciesB) ) {
 			$GU->user_info( 3, $speciesA."\t".$speciesB."\tthreshold already set\n" );
 			return;
 		}
 		else {
 			$GU->user_info( 3, $speciesA." vs ".$speciesB."\tvalue not set, trying to infer\n" );
	
			my $lca = $self->evolutionary_tree->last_common_ancestor($speciesA, $speciesB);
 			# retrieve thresholds for set X, Y and Z
 			
 			my @left_tree = $lca->get_subtree_lineage_species($speciesA);
 			$GU->user_info ( 3, Dumper (\@left_tree) );
 			my @right_tree = $lca->get_subtree_lineage_species($speciesB);
 			$GU->user_info( 3, Dumper (\@right_tree) );
 			my @other_trees = $lca->get_null_subtree_species($speciesA, $speciesB);
 			
 			my @setX = (); # left vs right thresholds

 			foreach my $leftspecies (@left_tree) {
 				foreach my $rightspecies (@right_tree) {
 					if ($self->is_set($leftspecies, $rightspecies)) {
 						my @threshold_to_add = $self->get_thresholds($leftspecies, $rightspecies);
 						push @setX, [ @threshold_to_add ];

 					}
 				}
 			}
 			# remove duplicates from @setX
 			my %seenX_values;
			my @setX_unique = grep { !($seenX_values{$_->[0]}{$_->[1]}++) } @setX;
			
			my @setY; # left vs other thresholds
 			foreach my $leftspecies (@left_tree) {
 				foreach my $otherspecies (@other_trees) {
 					if ($self->is_set($leftspecies, $otherspecies)) {
 						my @threshold_to_add = $self->get_thresholds($leftspecies, $otherspecies);
 						push @setY, [ @threshold_to_add ];

 					}
 				}
 			}
 			my %seenY_values;
			my @setY_unique = grep { !($seenY_values{$_->[0]}{$_->[1]}++) } @setY;
 			
 			my @setZ; # right vs other
 			foreach my $rightspecies (@right_tree) {
 				foreach my $otherspecies (@other_trees) {
 					if ($self->is_set($rightspecies, $otherspecies)) {
 						my @threshold_to_add = $self->get_thresholds($rightspecies, $otherspecies);
 						push @setZ, [ @threshold_to_add ];

 					}
 				}
 			}
 			my %seenZ_values;
			my @setZ_unique = grep { !($seenZ_values{$_->[0]}{$_->[1]}++) } @setZ;
 		
 			my $Xsetsize = @setX_unique;
 			
 			
 			# if X = 1: return result

			if ( $Xsetsize == 1 ) {
				$GU->user_info( 3,  $setX_unique[0][0]."\t". $setX_unique[0][1]."\n" );
			 	$self->set_matrix_value_for_pair($speciesA, $speciesB, $setX_unique[0][0], $setX_unique[0][1]);
			 	$GU->user_info( 3, $speciesA.":\t".$speciesB.":\t".$setX_unique[0][0]."\t".$setX_unique[0][1]."\t X=1, threshold pair set\n" );
			 	return;
			}
			# if X = 0
			if ($Xsetsize == 0) {
				$GU->user_info( 3, "X=0, no comparisons available, checking YZ set\n" );
				
				my @YZset = (@setY_unique, @setZ_unique);
				
				
				my %seenYZ_values;
				my @YZset_unique = grep { !($seenYZ_values{$_->[0]}{$_->[1]}++) } @YZset;
				
				my $YZsetsize = @YZset_unique;
				
				if ($YZsetsize == 0 ) {
					$GU->user_info( 3, "YZ set is also zero, cannot assign thresholds\n" );
					return;
				}
			# if YuZ = 1: return result
				if ($YZsetsize == 1) {
					$self->set_matrix_value_for_pair($speciesA, $speciesB, $YZset[0][0], $YZset[0][1]);
					$GU->user_info( 3, $speciesA.":\t".$speciesB.":\t".$YZset[0][0]."\t".$YZset[0][1]."\t YuZ = 1, threshold pair set\n" );
					return;
				}
			# if YuZ > 1: not resolvable
				if ($YZsetsize > 1) {
					$GU->user_info( 3, "YuZ >1, thresholds not resolvable\n" );
					return;
				}
			}
			# if X > 1
			if ( $Xsetsize > 1 ) {
				$GU->user_info( 3, "X > 1, trying to resolve\n" );
				# compute lists of species for each 'level' of tree, relative to each species
				my @relativeA = reverse ($self->evolutionary_tree->get_species_levels($speciesA));
				$GU->user_info( 3, Dumper (\@relativeA) );
				my @relativeB = reverse ($self->evolutionary_tree->get_species_levels($speciesB));					
				
				$GU->user_info( 3, "relative species lists". @relativeA . " ".@relativeB."\n" );

 				$GU->user_info( 3, Dumper (\@relativeB) );
				# initialise matrix (@AoA)
				my @tree_level_threshold_matrix = ();
				# initialise border matrix, and set all values FALSE
				my @border_matrix = ();
				foreach my $x (0..@relativeA-1) {
					foreach my $y (0..@relativeB-1) {
						$border_matrix[$x][$y] = FALSE;
					}
				}
				$GU->user_info( 3, Dumper (\@border_matrix) );
				# fill new matrix of thresholds
			
				# look up each pair of thresholds and store non-duplicates in @t_l_t_m
				for (my $k = 0; $k <= $#relativeA; $k++) {
					for (my $l = 0; $l <= $#relativeB; $l++) {
						my @cell_values;
						my @relA = @{$relativeA[$k]};
						my @relB = @{$relativeB[$l]};
						$GU->user_info( 3, "getting cell values:".$k." ".$l."\n" );
						foreach my $A (@relA) {
							foreach my $B (@relB) {
								$GU->user_info( 3, "look up ".$A." vs ".$B."\n" );
								if ( $self->is_set($A, $B) ) {
									my @threshold = $self->get_thresholds($A, $B);
									push (@cell_values, [ @threshold ]);
								}
							} 
						}
						# remove any duplicates
						my %seen;
						my @non_duplicate_cell_values = grep { !($seen{$_->[0]}{$_->[1]}++) } @cell_values;
						$GU->user_info( 3, Dumper (\@non_duplicate_cell_values) );
						# set values in tree matrix
						$tree_level_threshold_matrix[$k][$l] = [ @non_duplicate_cell_values ];
					}
				}	
				$GU->user_info( 3, Dumper (\@tree_level_threshold_matrix) );
				# set border matrix values to TRUE where values exist
				for (my $k = 0; $k <= $#tree_level_threshold_matrix; $k++) {
					$GU->user_info( 3, "k:".$k."\n" );
					for (my $l = 0; $l <= $#{$tree_level_threshold_matrix[$k]}; $l++) {
						$GU->user_info( 3, "l:".$l."\n" );
						
						if (@{$tree_level_threshold_matrix[$k][$l]} > 0) {
							$border_matrix[$k][$l] = TRUE;
						}
					}
				}
				# re-set, column-wise, values after a TRUE to FALSE
				
				for (my $k = 0; $k <= $#border_matrix; $k++) {
					my $reset = FALSE;
					for (my $l = 0; $l <= $#border_matrix; $l++) {
						if ($reset) {
							$border_matrix[$k][$l] = FALSE;
						}
						if ($border_matrix[$k][$l]) {
							# set values after to FALSE
							$reset = TRUE;
						}
						
					}
				}
				
				# same thing but row-wise
				
				for (my $k = 0; $k <= $#border_matrix; $k++) {
					my $reset = FALSE;
					for (my $l = 0; $l <= $#border_matrix; $l++) {
						if ($reset) {
							$border_matrix[$l][$k] = FALSE;
						}
						if ($border_matrix[$l][$k]) {
							# set values after to FALSE
							$reset = TRUE;
						}
						
					}
				}
				
				# get border values
				my @border_values;
				for (my $k = 0; $k <= $#border_matrix; $k++) {
					for (my $l = 0; $l <= $#border_matrix; $l++) {
						if ($border_matrix[$k][$l]) {
							push (@border_values, @{$tree_level_threshold_matrix[$k][$l]} );
						}
					}
				}
				# get unique border values
				my %seen;
				my @unique_border_values = grep { !($seen{$_->[0]}{$_->[1]}++) } @border_values;
				
				# if one value, assign
				if (@unique_border_values == 1) {
					$GU->user_info( 3, $speciesA."\t". $speciesB."\t" ).								$unique_border_values[0][0]."\t". $unique_border_values[0][1]. " unique, have set!\n";
					
					$self->set_matrix_value_for_pair($speciesA, $speciesB, 								$unique_border_values[0][0], $unique_border_values[0][1])
				}
				else {
					$GU->user_info( 3, "conflicting values, threshold cannot be inferred\n" );
				}
			}
 		}
 		return;
 	} # infer_pair_threshold #

} # Partial_Threshold_Matrix #
