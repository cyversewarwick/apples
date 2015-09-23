### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Partial_Threshold_Matrix_Maker Class ###
### need to fully parameterie the constructor so that user can set any pairwise threshold they like 

use MooseX::Declare;

class Partial_Threshold_Matrix_Maker {
	use Partial_Threshold_Matrix;

	use APPLES_Datatypes qw (EvolutionaryTreeKingdom);
	use constant {FALSE => 0,
		      TRUE  => 1};

	my $GU = General_Utilities->new();

	method make_partial_threshold_matrix(Evolutionary_Tree $tree, EvolutionaryTreeKingdom $kingdom) {
	  my $result;
	  # initialise a thresholds for some pairs of species
	  if ($kingdom eq 'plants') {
	    $result = $self->private_initialise_partial_threshold_matrix_plants ($tree);
	  }
	  elsif ($kingdom eq 'vertebrates') {
	    $result = $self->private_initialise_partial_threshold_matrix_vertebrates ($tree);
	  }
        elsif ($kingdom eq 'invertebrates') {
        $result = $self->private_initialise_partial_threshold_matrix_invertebrates ($tree);
        }
	  else {
	    die 'unknown kingdom, or cannot handle that kingdom setting yet!';
	  }
	  # run autofill-function through Job_Handler (for caching rather than parallelisation)
	  my $cache = TRUE;
	  my $cache_duration = 30;
	  my $job_handler = Job_Handler->new();
	  $job_handler->get_config_settings();
	  my $function = 'autofill';
	  my @parameters = ();
	  my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
						   wall_time_estimate => 172800,
						   recache => FALSE
	      );
	  my @cache_result = eval {
	      $job_handler->handle_APPLES_function($function, $result, \@parameters, $cache, $cache_duration, $job_parameters);
	  };
	  if ($@) {
	      my $exception_content = $@;
	      if (!UNIVERSAL::can($exception_content, 'isa')) {
		  die $exception_content;
	      }
	      else {
		  if ($exception_content->isa('Job_Information_Exception')) {
		      my $aggregate_exception = Aggregate_Exception->new();
		      $aggregate_exception->merge($exception_content);
		      $aggregate_exception->print_statistics_and_die();
		  }
		  else { 
		      die $exception_content;
		  }
	      }
	  }
	  $result = shift(@cache_result);
	  # check completeness (that is: thresholds for all pairs of species can be automatically inferred
          # from the given thresholds set below)
	  if (!$result->is_complete()) {
	      die 'partial threshold matrix is incomplete for this tree - provide more thresholds';
	  }
	  return $result;
	} # make_partial_threshold_matrix #

	method private_initialise_partial_threshold_matrix_plants(Evolutionary_Tree $tree) {	  
	    my $PTM = Partial_Threshold_Matrix->new(evolutionary_tree => $tree);	    
	    $PTM->initialise_matrix();
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'maize', 60, 70);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'tomato', 60, 70);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'cotton', 60, 70);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'grape', 62, 82);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('brassica', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('papaya', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('cotton', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('medicago', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('soybean', 'poplar', 62, 78);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'brassica', 70, 80);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'papaya', 62, 75);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'medicago', 60, 70);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'rice', 60, 70);
	    $PTM->set_matrix_value_for_pair('arabidopsis', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('brassica', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('papaya', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('cotton', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('medicago', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('soybean', 'castor', 60, 70);
	    $PTM->set_matrix_value_for_pair('glaberrima', 'indica', 60, 70);
	    $PTM->set_matrix_value_for_pair('glaberrima', 'maize', 60, 70);
	    $PTM->set_matrix_value_for_pair('medicago', 'soybean', 60, 70);
	    $PTM->set_matrix_value_for_pair('maize', 'sorghum', 60, 70);
	    $PTM->set_matrix_value_for_pair('castor', 'poplar', 60, 70);
	    return $PTM;
	} # private_initialise_partial_threshold_matrix_plants #

	method private_initialise_partial_threshold_matrix_vertebrates(Evolutionary_Tree $tree) {
	  my $PTM = Partial_Threshold_Matrix->new(evolutionary_tree => $tree);  
	  $PTM->initialise_matrix();
	  
	  $PTM->set_matrix_value_for_pair('human', 'mouse', 83, 99);
	  $PTM->set_matrix_value_for_pair('mouse', 'rat', 90, 110);
	  $PTM->set_matrix_value_for_pair('mouse', 'chimpanzee', 83, 99);
	  $PTM->set_matrix_value_for_pair('mouse', 'dog', 83, 99);
	  $PTM->set_matrix_value_for_pair('mouse', 'cow', 83, 99);
	  $PTM->set_matrix_value_for_pair('mouse', 'opossum', 64, 80);
	  $PTM->set_matrix_value_for_pair('mouse', 'chicken', 63, 77);
	  $PTM->set_matrix_value_for_pair('mouse', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('mouse', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('mouse', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('mouse', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('human', 'chimpanzee', 99, 120);
	  $PTM->set_matrix_value_for_pair('human', 'dog', 83, 99);
	  $PTM->set_matrix_value_for_pair('human', 'cow', 83, 99);
	  $PTM->set_matrix_value_for_pair('human', 'opossum', 64, 80);
	  $PTM->set_matrix_value_for_pair('human', 'chicken', 63, 77);
	  $PTM->set_matrix_value_for_pair('human', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('human', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('human', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('human', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('dog', 'cow', 83, 99);
	  $PTM->set_matrix_value_for_pair('dog', 'opossum', 64, 80);
	  $PTM->set_matrix_value_for_pair('dog', 'chicken', 63, 77);
	  $PTM->set_matrix_value_for_pair('dog', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('dog', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('dog', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('dog', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('cow', 'opossum', 64, 80);
	  $PTM->set_matrix_value_for_pair('cow', 'chicken', 63, 77);
	  $PTM->set_matrix_value_for_pair('cow', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('cow', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('cow', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('cow', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('opossum', 'chicken', 63, 77);
	  $PTM->set_matrix_value_for_pair('opossum', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('opossum', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('opossum', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('opossum', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('chicken', 'xenopus', 63, 75);
	  $PTM->set_matrix_value_for_pair('chicken', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('chicken', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('chicken', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('xenopus', 'fugu', 61, 72);
	  $PTM->set_matrix_value_for_pair('xenopus', 'tetraodon', 61, 72);
	  $PTM->set_matrix_value_for_pair('xenopus', 'zebrafish', 61, 72);
	  $PTM->set_matrix_value_for_pair('fugu', 'tetraodon', 83, 99);
	  $PTM->set_matrix_value_for_pair('fugu', 'zebrafish', 70, 90);
          # add further thresholds as appropriate	  
	
	  return $PTM;
	} # private_initialise_partial_threshold_matrix_vertebrates #
    
    method private_initialise_partial_threshold_matrix_invertebrates(Evolutionary_Tree $tree) {
        my $PTM = Partial_Threshold_Matrix->new(evolutionary_tree => $tree);
        $PTM->initialise_matrix();
        
        $PTM->set_matrix_value_for_pair('nasonia vitripennis', 'apis mellifera', 83, 99);
    }

} # Partial_Threshold_Matrix_Maker #
