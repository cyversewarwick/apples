### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReRe_Set_Maker class ###
use MooseX::Declare;

class ReRe_Set_Maker {
  use Parameters;
  use Reg_Loc;
  use ReRe_Maker;
  use Reg_Loc_Maker;
  use Data::Dumper;
  use Exception;
  use APPLES_Datatypes qw (Boolean);
  use constant {FALSE => 0,
		TRUE	=> 1};
  use General_Utilities;
  
  my $GU = General_Utilities->new();
  
  method rere_set_maker(ReRe_Set_Constructor_Parameters $parameters, ArrayRef[Reg_Loc] $reg_loc_list_ref, Boolean $each_rere_through_job_handler) {
      my @rere_set_members = ();
      my $rere_maker = ReRe_Maker->new();
      my $exception_info = Aggregate_Exception->new();
      my $stats_required = FALSE;
      my $progress_count = 0;
      my $number_of_reg_locs = @{$reg_loc_list_ref};
      foreach my $reg_loc (@{$reg_loc_list_ref}) {
	  $progress_count++;
	  $GU->user_info(2,$progress_count."/".$number_of_reg_locs.": ".$reg_loc->gene_ID."\n");
	  my $rere;
	  if ($each_rere_through_job_handler) {
	      eval {
		  $rere = $rere_maker->rere_maker_through_job_handler($parameters->rere_constructor_parameters, $reg_loc);
	      };
	      if ($@) {
		  my $exception_content = $@;
		  my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$exception_info);
		  if ($outcome[0]) {
		      $stats_required = TRUE;
		  }
		  $exception_info = $outcome[1];
	      }
	      else {
		  push (@rere_set_members, $rere);
	      }
	  }
	  else {	
	      $rere = $rere_maker->rere_maker($parameters->rere_constructor_parameters, $reg_loc);
	      push (@rere_set_members, $rere);
	  }	 
      }
      if ($stats_required eq TRUE) {
	  $exception_info->print_statistics_and_die;
      }
      my $rere_set = ReRe_Set->new(rere_set_members => \@rere_set_members);
      return $rere_set;
  } # rere_set_maker #

  method rere_set_maker_through_job_handler (ReRe_Set_Constructor_Parameters $parameters, ArrayRef[Reg_Loc] $reg_loc_list_ref, Boolean $each_rere_through_job_handler) {
      my $object = ReRe_Set_Maker->new();
      my $function = 'rere_set_maker';
      my @parameters = ($parameters,$reg_loc_list_ref,$each_rere_through_job_handler);
      my $high_memory = TRUE;
      my @result = $GU->standard_call_to_job_handler($object,$function,\@parameters,$high_memory,FALSE);
      return $result[0];
  } # rere_set_maker_through_job_handler #
  
  method make_rere_set_for_all_reg_locs_in_a_genome(ReRe_Set_Constructor_Parameters $parameters,
						    Genome_Sequence_Database_Parameters $genome_db) {
      my $reg_loc_maker = Reg_Loc_Maker->new();
      my @genome_wide_reg_locs = $reg_loc_maker->make_all_reg_locs_for_a_genome($genome_db);
      my $rere_set = $self->rere_set_maker($parameters,\@genome_wide_reg_locs,FALSE);
      return $rere_set;
  } # make_rere_set_for_all_reg_locs_in_a_genome #

  method make_rere_set_for_all_reg_locs_in_a_genome_through_job_handler(ReRe_Set_Constructor_Parameters $parameters,
									Genome_Sequence_Database_Parameters $genome_db) {
      my $rere_set_maker = ReRe_Set_Maker->new();
      my $cache = TRUE;
      my $cache_duration = 360;
      my $job_handler = Job_Handler->new();
      $job_handler->get_config_settings();
      my $function = 'make_rere_set_for_all_reg_locs_in_a_genome';
      push (my @parameters, $parameters,$genome_db);
      my $job_parameters = Job_Parameters->new(memory_requirement_high => TRUE,
					       wall_time_estimate => 172800
	  );
      my $aggregate_exception = Aggregate_Exception->new();
      my @cache_result = eval {
	  $job_handler->handle_APPLES_function($function, $rere_set_maker, \@parameters, $cache, $cache_duration, $job_parameters);
      };
      if ($@) {
	  my $exception_content = $@;
	  my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
	    # no parallelism here, though
	  $aggregate_exception = $outcome[1];
	  if ($outcome[0]) {
	      $aggregate_exception->print_statistics_and_die;
	  } else {
	      die 'did not expect this could happen at design time of this code.';
	  }
      }
      my $rere_set = $cache_result[1];
      return $rere_set;
  } # make_rere_set_for_all_reg_locs_in_a_genome_through_job_handler #

} # ReRe_Set_Maker #
  
