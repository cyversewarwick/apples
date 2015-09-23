### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReMo_Set class ###
use MooseX::Declare;
use Serialization::Serializable;

class ReMo_Set extends Serialization::Serializable {
  use ReMo;
  use General_Utilities;
  use Genome_DB_Utilities;
  use APPLES_Datatypes qw (APPLESSpeciesName PositiveNum);

  use Moose;
  has 'remo_set' => (is => 'rw', isa => 'ArrayRef[ReMo]', required => 1); 
  
  my $GU = General_Utilities->new();
  my $GDBU = Genome_DB_Utilities->new();

  # Thresholds the ReRe based on a specific beliefe score. Looks at the top ReMo and checks that the belief score is 
  # greater or equal to the threshold.
  # 
  # Input:     $belief_threshold - PositiveNum sepcifying the belief score threshold.
  # 
  # Output:    TRUE if ReMo set passes the threshold, FALSE otherwise.	
  method threshold(PositiveNum $belief_threshold) {
  	my $accept = 0;
  	my $remo = @{$self->remo_set}[0];
  	
  	if($remo->belief_score >= $belief_threshold) {
  		$accept = 1;
  	}
  	
  	return $accept;
  }

  method render_text {
    print "\nremo_set:\n";
    foreach my $remo ( @{$self->remo_set} ) {
      $remo->render_text();
    }
  } # render_text # 
    
    method render_tabular_for_one_species (APPLESSpeciesName $species, Str $core_reg_loc_id, Int $core_reg_loc_position) {
      foreach my $remo ( @{$self->remo_set} ) {
	$remo->render_tabular_for_one_species($species, $core_reg_loc_id, $core_reg_loc_position);
      }
    } # render_tabular_for_one_species #

  method number_of_remos() {
    my $result = @{$self->remo_set};
    return $result;
  } # number_of_remos #

  method filter_remo_set_by_db_IDs (ArrayRef[Str] $db_IDs) {
    # produces a new ReMo_Set containing sequences from only the databases specified by user
   
    my @remos;
    foreach my $remo (@{$self->remo_set}) {
      my $remo_DB_ID = $GDBU->get_ID_for_database($remo->genome_db);
      my @remo_DB_ID_as_list = ();
      push(@remo_DB_ID_as_list,$remo_DB_ID);
      my $lists_overlap = $GU->lists_overlap($db_IDs,\@remo_DB_ID_as_list);
      if ($lists_overlap) {
	push(@remos,$remo);
      }
    }
    my $result = ReMo_Set->new(remo_set => \@remos);
    return $result;
  } # filter_remo_set_by_db_IDs #

} # ReMo_Set #
