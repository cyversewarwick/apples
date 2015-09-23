### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Motif_Finder class: methods to predict over-represented motifs ###
use MooseX::Declare;

class Motif_Finder {
	use Parameters;
	use Data::Dumper;
	use Bio::SeqIO;
	use General_Utilities;
	use APPLES_Datatypes qw (Boolean);
	use File::Temp qw (tempdir);
	use constant {FALSE => 0,
				  TRUE => 1};

	my $GU = General_Utilities->new();
	
	method run_MEME_motif_finder( ArrayRef[Genomic_Interval_Set] $genomic_interval_sets, MEME_Parameters $meme_parameters) {
		
	
		
		my @all_motif_finding_results;
		
		my $aggregate_exception = Aggregate_Exception->new();
		my $stats_required = FALSE;
		my $meme_wrapper =  MEME_wrapper->new();		
		foreach my $genomic_interval_set ( @$genomic_interval_sets ){
					
			# CALL MEME THROUGH THE JOB HANDLER:
			my $cache = TRUE; 
			my $cache_duration = 10; # number of days
			my $job_handler = Job_Handler->new(); 
			$job_handler->get_config_settings();
			my $function = 'run'; 
			push (my @parameters, $genomic_interval_set, $meme_parameters); 
			my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
														wall_time_estimate => 172800
													);
			
			# CALL THE FUNCTION VIA THE JOB_HANDLER - WRAP IN AN EVAL STATEMENT TO CATCH ERRORS PROPERLY (thrown by Job_Handler as Job_Information_Exception objects) 

			my @cache_result = eval {
				$job_handler->handle_APPLES_function($function, $meme_wrapper, \@parameters, $cache, $cache_duration, $job_parameters);
			};			
			if ($@) {
				my $exception_content = $@;
				my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
				if ($outcome[0]) {
					$stats_required = TRUE;
				}
				$aggregate_exception = $outcome[1];
			}
			else {
				my $motif_finding_result_object = $cache_result[1]; # the remaining array is the result of the function you called
				push( @all_motif_finding_results, $motif_finding_result_object );
			}
			
		}

		if ($stats_required) {
		  die $aggregate_exception;
		}
	 	return \@all_motif_finding_results;
	} # run_MEME_motif_finder #
	
}	
