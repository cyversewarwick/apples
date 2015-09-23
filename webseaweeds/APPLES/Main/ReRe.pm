### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReRe class ###
# A ReRe is a group of ReMo_Sets for one specified regulated locus ("regulatory region")
use MooseX::Declare;
use Serialization::Serializable;

class ReRe extends Serialization::Serializable {
  use Bio::SeqIO;
  use General_Utilities;
  use Parameters;
  use ReMo_Reg_Loc_Assignment;
  use Reg_Loc;
  use ReMo_Set;
  use Data::Dumper;
  use APPLES_Datatypes qw (APPLESSpeciesName PositiveNum Boolean);

  use Moose;
  has 'core_reg_loc' => (is => 'rw', isa => 'Reg_Loc', required => 1);
  has 'remo_sets' => (is => 'rw',
		      isa => 'ArrayRef[ReMo_Set]',
		      required => 1);
  has 'ns__homologous_reg_loc_list_ref' => (is => 'rw', isa => 'ArrayRef[Reg_Loc]'); # just an idea, not yet implemented
  has 'ns__remo_reg_loc_assignment_list_ref' => (is => 'rw', isa => 'ArrayRef[ReMo_Reg_Loc_Assignment]');  # just an idea, not yet implemented
  
  my $GU = General_Utilities->new();


  # Simple method for creating a FASTA file from this ReRe. The sequence is written to file specified.
  #
  # Input:		$filename - Str filename where the fasta data is stored (".fasta" is appended automatically).
  #             $details - Boolean to specify if details should be written too.
  #
  # Returns:	1 - if all was success.
  method create_fasta_file(Str $filename, Boolean $details) {
  	
  	# for each of the ReMo's get the sequence.
  	foreach my $ReMo_Set ( @{$self->remo_sets}){
	    $GU->user_info( 3, $ReMo_Set."\n" );
	    my @array_of_remos = @{$ReMo_Set->remo_set};
	    
	    foreach my $remo (@array_of_remos){
	      
	      my $gi_sequence = $remo->get_sequence();
	      my $interval = $remo->five_prime_pos . ':' . $remo->three_prime_pos;
	      $GU->user_info( 3, $remo->get_repeatmasked_sequence()."\n" ); 
	      
	      
	      my $desc;
	      if($details) {
	      	$desc = "dbname=".$remo->genome_db->dbname." conservation=".$remo->conservation." belief_score=".$remo->belief_score;
	      } else {
	      	$desc = "";
	      }
	      my $seq = Bio::Seq->new( -seq => $gi_sequence,
			     -id => $remo->label.'('.$interval.')',
			     -description => $desc
			     );
  	
  		  # N.B.: '>>' means append to the file.
  		  my $output_stream = Bio::SeqIO->new(-file => '>>'.$filename.'.fasta', -format => 'Fasta');
    	  $output_stream->write_seq($seq);
	      
	    }#foreach
	  }#foreach
    
    return 1;
  }

  # Thresholds the ReRe based on a specific beliefe score. This iterates over all remo_sets and throws out the once that
  # are not accepted using ReMo->threshold method.
  # Input:     $belief_threshold - PositiveNum sepcifying the belief score threshold.	
  # Output:    TRUE if the ReRe is accepted, FALSE otherwise.
  method threshold(PositiveNum $belief_threshold) {
  	my @thresholded_remos;
  	
  	my $remo_sets_length = @{$self->remo_sets};
  	my $index;
  	for($index = 0; $index < $remo_sets_length; $index++ ) {
  		my $ReMo_Set = @{$self->remo_sets}[$index];
  		
  		my $accepted = $ReMo_Set->threshold($belief_threshold);
  		
  		if($accepted) {
  			push(@thresholded_remos, $ReMo_Set);
  		}
  	
  	} #foreach ReMo_Set
  	
  	$self->{remo_sets} = \@thresholded_remos;
  	
  	$remo_sets_length = @{$self->remo_sets};
  	
  	if($remo_sets_length > 0) {
  		return 1;
  	} else {
  		return 0;
  	}
  }

  method number_of_remo_sets() {
    my $result = @{$self->remo_sets};
    return $result;
  } # number_of_remo_sets #

  method give_rere_species() {
    
    die 'method not implemented yet!';
    
    # Input all ReMos that make up the ReRe
    # Check what the species assignment for each ReMo is
    # Returns the species if they are all the same,
    # Else returns an array/error/dies if they do not all match

  } # give_rere_species #

    method get_rere_sequence() { # adds N-spacers in between discontinuous genomic intervals to form one sequence
      
      die 'method not implemented yet!';
      
    } # get_rere_sequence #
      
  method render_text() {
      print "Core Reg_Loc:\t".$self->{core_reg_loc}->gene_ID." ***\n";
      foreach my $remo_set ( @{$self->{remo_sets}} ) {
	$remo_set->render_text();
      }
  } # render_text #
	
  method render_tabular_for_one_species (APPLESSpeciesName $species) {
    my $core_reg_loc_id = $self->{core_reg_loc}->gene_ID;
    my $core_reg_loc_position = $self->{core_reg_loc}->position;
    if (scalar(@{$self->{remo_sets}}) == 0) {
      print $core_reg_loc_id."\tno remos found\n";
    }
    foreach my $remo_set ( @{$self->{remo_sets}} ) {
      $remo_set->render_tabular_for_one_species($species, $core_reg_loc_id, $core_reg_loc_position);
    }
  } # render_tabular_for_one_species # 

  method make_gi_set_from_rere (Boolean $turn_to_positive_strand) {
      # function takes all the remos in the ReRe (may be from multiple species), and creates a
      # Genomic_Interval_Set containing non-overlapping intervals (overlaps are merged)
      
      # first make an array of all the intervals (which may be overlapping)
      my @remos;
      foreach my $remo_set ( @{$self->{remo_sets}} ) {
	  foreach my $remo ( @{$remo_set->remo_set} ) {
	      push (@remos, $remo);
	  }
      }
      my $gi_set_maker = Genomic_Interval_Set_Maker->new();
      my $gi_set = $gi_set_maker->make_gi_set_from_gi_array(\@remos);
      my $result = $gi_set->merge($turn_to_positive_strand);
      return $result;
  } # make_gi_set_from_rere #

  method filter_rere_by_db_IDs (ArrayRef[Str] $db_IDs) {
    # produces a new rere containing sequences from only the databases specified by user
   
    my @remo_sets;
    foreach my $remo_set (@{$self->remo_sets}) {
      my $filtered_remo_set = $remo_set->filter_remo_set_by_db_IDs($db_IDs);
      if ($filtered_remo_set->number_of_remos() > 0) {
	push(@remo_sets,$filtered_remo_set);
      }
    }
    my $result = ReRe->new(core_reg_loc => $self->core_reg_loc,
			   remo_sets => \@remo_sets);
    return $result;
  } # filter_rere_by_db_IDs # 

} # ReRe #
