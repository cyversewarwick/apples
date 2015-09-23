### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### MEME_wrapper Class ###
# wrapper for MEME binary

use MooseX::Declare;

class MEME_wrapper {
  use Genomic_Interval_Set;
  use General_Utilities;
  use Motif_Finding_Result_Maker;
  use File::Path qw (rmtree);
  use constant {FALSE => 0,TRUE => 1};

  my $GU = General_Utilities->new();

  method run (Genomic_Interval_Set $genomic_interval_set, MEME_Parameters $meme_parameters) {
      
      # Get unique temporary directory
      my $tempdir = $GU->get_temp_random_directory(FALSE);
      my $running_dir = $tempdir;
      # Set this temp dir as the running/output dir for this MEME job
      $meme_parameters->{-oc} = $running_dir;
      
      # Get MEME binary location
      my $job_handler = Job_Handler->new(); # Instantiate a Job_Handler object 
      $job_handler->get_config_settings();
      my $run_command = $job_handler->job_handler_settings()->{meme_binary};
      
      my $fasta_filename = $running_dir."meme_input.fasta";
      my $meme_input_filepath = $genomic_interval_set->create_FASTA_from_genomic_interval_set($fasta_filename);
      
      $run_command .= ' '.$meme_input_filepath;
      
      my @meme_arguments = keys %{$meme_parameters};
      
      foreach my $arg (@meme_arguments) {
	  # If argument is -revcomp or -pal or -dna then because they have no required parameters do not print the value associated with that key
	  $GU->user_info ( 3, "ARG= $arg\n" );
	  if( ($arg eq '-revcomp') or  ($arg eq '-pal') or ($arg eq '-dna') ) {
	      if ($$meme_parameters{$arg}) {
		  $run_command .= " $arg";	
	      }
	  }
	  else{  
	      $run_command .= " $arg $$meme_parameters{$arg}";	 
	  	}			
      }
   
      $GU->user_info ( 1, "Running MEME now\n" );
      $GU->user_info( 1, $run_command."\n" ); 
      #system(`$run_command`) == 0 || die "System error!";
      system(`$run_command`); # Sascha: took out check of return code as there were cases where the return code was unequal zero
                              #         even though MEME ran successfully
      
      $GU->user_info( 1, "Returning meme output file address\n" );
      my $meme_output_filepath = $$meme_parameters{-oc}.'/meme.txt';# Needs to change!
     	
      
      my $motif_finding_result_maker = Motif_Finding_Result_Maker->new();
      my $motif_finding_result = $motif_finding_result_maker->make_motif_finding_result_from_meme_text_output($meme_output_filepath);  
      
     # rmtree($tempdir);
     system("rm -rf ".$tempdir);

      return $motif_finding_result;
  } # run #
  
} # MEME_wrapper #
