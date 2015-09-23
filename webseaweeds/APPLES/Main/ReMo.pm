### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReMo class ###
use MooseX::Declare;

class ReMo extends Genomic_Interval {
  use Parameters;
  use APPLES_Datatypes qw (APPLESSpeciesName);

  use Moose;
  has 'belief_score' => (is => 'ro', isa => 'Num');
  has 'repeat_ratio' => (is => 'ro', isa => 'Num');
  has 'conservation' => (is => 'ro', isa => 'Num');
  
  method get_bisis {
    
    die "method not implemented yet!\n";
    
    my $ReMo_BiSi_List;
    return $ReMo_BiSi_List;
    
  } # get_bisis #
    
    method render_text {
      
      print $self->{genome_db}->alias." ".$self->coord_sys_name." ".$self->region." ".$self->five_prime_pos."-".$self->three_prime_pos." (".$self->strand." strand)";
      if (defined $self->label) {
	print " ".$self->label;
      }
      
      if (defined $self->belief_score) {
	printf " \[belief score/repeat ratio/conservation: %.2f %s %s\n", $self->belief_score, $self->repeat_ratio, $self->conservation."\]";
      }
      print $self->get_working_sequence()."\n";

    } # render_text #
    
      method render_tabular_for_one_species (APPLESSpeciesName $species, Str $core_reg_loc_id, Int $core_reg_loc_position) {
	if ($self->{genome_db}->alias eq $species) {
	  print $core_reg_loc_id."\t";
	  if (defined $self->belief_score) {
	    my $rounded_belief_score = sprintf("%.2f", $self->belief_score);
	    my $rounded_repeat_ratio = sprintf("%.2f", $self->repeat_ratio);
	    print $rounded_belief_score."\t". $rounded_repeat_ratio."\t". $self->conservation."\t";
	  }
	  else {
	    print "-\t-\t-\t";
	  }
	  my $seq = $self->get_working_sequence();
	  my $seq_length = length($seq);
	  print $seq_length."\t";
	  
	  my $distance = $core_reg_loc_position - $self->three_prime_pos;
	  if ($self->strand eq 'negative') {
	    $distance = $self->three_prime_pos - $core_reg_loc_position;
	  }
	  print $distance."\n";
	}
      } # render_tabular_for_one_species #

} # ReMo #
