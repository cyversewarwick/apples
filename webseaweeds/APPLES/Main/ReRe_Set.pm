### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReRe_Set class ###
use MooseX::Declare;
use Serializable;

class ReRe_Set extends Serializable {
  use Parameters;
  use ReRe;	
  use General_Utilities;
  use Data::Dumper;
  use APPLES_Datatypes qw (APPLESSpeciesName PositiveInt PositiveNum Boolean);
  
  has 'rere_set_members' => (is => 'rw', isa => 'ArrayRef[ReRe]'); #ArrayRef[ReRe]
  
  my $GU = General_Utilities->new();
  
  # This method creates a single FASTA file for the whole ReRe Set by appending to the specified file name.
  # Input:     $output_fasta_filename - Str specifying the output file name.
  #            $details - Boolean to specify if details should be written too.
  #
  # Return: Returns 1 if the operation was successful.              
  method create_FASTA_from_rere_set(Str $output_fasta_filename, Boolean $details){
    $GU->user_info( 1, "Creating FASTA file for ReRe set\n" );
    my $agi_gene_id;
    my $gi_sequence;
    foreach my $rere (@{$self->rere_set_members}) {     
      $rere->create_fasta_file($output_fasta_filename, $details);
      
    }#foreach
    
    return 1;
    
  } # create_FASTA_from_rere_set #
  
  # Thresholds the ReRe set based on a specific beliefe score.
  #Input:     $belief_threshold - PositiveNum sepcifying the belief score threshold.	
  method threshold(PositiveNum $belief_threshold) {
  	my @thresholded_reres;
  	foreach my $rere ( @{$self->rere_set_members} ) {
  		my $accepted = $rere->threshold($belief_threshold);
	
  		if($accepted) {
  			push(@thresholded_reres, $rere);
  		}
  	}
  	
  	 $self->{rere_set_members} = \@thresholded_reres;
  	
  }
  
  method render_text() {
    my $counter = 1;
    print "\nrendering ReRe_Set as text\n";
    foreach my $rere ( @{$self->rere_set_members} ) {
      print "*** ReRe  ".$counter.", ";
      $rere->render_text();
      print "*** end of ReRe ".$counter." ***\n\n";
      $counter ++;
    }
  } # render_text #

    method render_simple_count(Int $conservation) {
      
      foreach my $rere ( @{$self->rere_set_members} ) {
	my $rere_count = 0;
	my $remoset_count = 0;
	my $remo_count = 0;
	my $remos_above_threshold = 0;
	$rere_count++;
	foreach my $remo_set ( @{$rere->{remo_sets}} ) {
	  $remoset_count++;
	  foreach my $remo ( @{$remo_set->remo_set} ) {
	    $remo_count++;
	    if ($remo->conservation >= $conservation) {
	      $remos_above_threshold++;
	    }
	  }
	}
	my $id = $rere->core_reg_loc->gene_ID;
	print "$id, $remoset_count, $remo_count, $remos_above_threshold\n";
      }
    #  print "rere count in this set = $rere_count\nremo set count (conserved regions in arabidopsis and at least one other species) = $remoset_count\nremo count (total number of conserved regions in all species) = $remo_count\nremos greater than or equal to conservation score($conservation) = $remos_above_threshold\n";
    } # render_simple_count #

  method render_conservation_scores {
    my @scores;
    foreach my $rere ( @{$self->rere_set_members} ) {
      foreach my $remo_set ( @{$rere->{remo_sets}} ) {
	  foreach my $remo ( @{$remo_set->remo_set} ) {
	    push (@scores, $remo->conservation);
	  }
	}
      }
    foreach my $score (@scores) {
      print $score."\n";
    }
  } # render_conservation_scores #
    
  method make_gi_set_from_rere_set(Boolean $turn_to_positive_strand) {
      # function takes all the remos in the ReRe (may be from multiple species), and creates a
      # Genomic_Interval_Set containing non-overlapping intervals (overlaps are merged)
    
    my $progress_count = 0;
    my $list_length = @{$self->rere_set_members};
    $GU->user_info(2, $list_length." reres to process:\n" );
      my @gi_sets;
      foreach my $rere (@{$self->rere_set_members}) {
	  my $one_gi_set = $rere->make_gi_set_from_rere($turn_to_positive_strand);
	  push(@gi_sets,$one_gi_set);
	  $progress_count++;
	  if (($progress_count % 10) == 0) {
	    $GU->user_info(2,$progress_count."/".$list_length." reres processed.\n");
	  }
      }
      my $result;
      if ($#gi_sets>=0) {
	  my $main_gi_set = $gi_sets[0];
	  
	  $main_gi_set->add_gi_sets(\@gi_sets); # this will duplicate genomic intervals of the first set
	  
	  $result = $main_gi_set->merge($turn_to_positive_strand);
      }
      else {
	  my @emptyarray = ();
	  my $gi_set_maker = Genomic_Interval_Set_Maker->new();
	  $result = $gi_set_maker->make_gi_set_from_gi_array(\@emptyarray);
      }
      return $result;
  } # make_gi_set_from_rere_set #

  method filter_rere_set_by_db_IDs (ArrayRef[Str] $db_IDs) {
    # produces a new rere_set containing sequences from only the databases specified by user
   
    my @reres;
    foreach my $rere (@{$self->rere_set_members}) {
      my $filtered_rere = $rere->filter_rere_by_db_IDs($db_IDs);
      if ($filtered_rere->number_of_remo_sets() > 0) {
	push(@reres,$filtered_rere);
      }
    }
    my $result = ReRe_Set->new(rere_set_members => \@reres);
    return $result;
  } # filter_rere_set_by_db_IDs #

  method first_n(PositiveInt $n) {
      # returns a ReRe_Set containing only the first $n ReRes

      my @all_reres = @{$self->rere_set_members};
      my $last = $n-1;
      my $length = @all_reres;
      if ($n>$length) {
	  $last = $length-1;
      }
      my @subset = ();
      if ($last>=0) {
	  @subset = @all_reres[0..$last];
      }
      my $new_rere_set = ReRe_Set->new(rere_set_members => \@subset);
      return $new_rere_set;
  } # first_n #

  method all_gene_ids() {
      # returns gene-IDs in same order as in $self->rere_set_members

      my @result;
      my @full_array = @{$self->rere_set_members};
      my $length = @full_array;
      for (my $i=0;$i<$length;$i++) {
	  my $member = $full_array[$i];
	  if (defined $member->core_reg_loc->gene_ID) {
	      my $one_gene_id = $member->core_reg_loc->gene_ID;
	      push(@result,$one_gene_id);
	  }
	  else {
	      die 'gene-ID for at least one ReRe not defined.';
	  }
      }
      return @result;
  } # all_gene_ids # 

  method render_tabular_for_one_species (APPLESSpeciesName $species) {
    # note - making this a dbname rather than APPLESSpeciesName would be more generic

    my $counter = 1;
    print "\nrendering ReRe_Set information for ".$species." in tabular format\n";
    foreach my $rere ( @{$self->rere_set_members} ) {
      #print "*** ReRe ".$counter.", ";
      $rere->render_tabular_for_one_species($species);
      #print "*** end of ReRe ".$counter." ***\n\n";
      $counter ++;
    }
  } # render_tabular_for_one_species #

} # ReRe_Set #
  
