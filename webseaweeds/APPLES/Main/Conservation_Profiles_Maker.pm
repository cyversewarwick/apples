### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Conservation_Profiles_Maker Class ###
# Constructor class to create window pair alignment conservation profiles (Conservation_Profiles objects)

use MooseX::Declare;

class Conservation_Profiles_Maker {
  use Conservation_Profiles;
  use Star_Bundler;
  use Data::Dumper;
  use APPLES_Datatypes qw (Boolean);
  
  use constant {FALSE => 0, TRUE => 1};
  use APPLES_Datatypes qw (Boolean);

  my $GU = General_Utilities->new();

  method make_conservation_profiles (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, ArrayRef $reg_loc_ref) {
    # makes conservation profiles for one sequence against a set of other sequences (using Star_Bundler)

    my @conservation_profiles;
    my $exception_info = Aggregate_Exception->new();
    my $stats_required = FALSE;

    foreach my $reg_loc (@{$reg_loc_ref}) {
      my $conservation_profile;
      eval {
	$conservation_profile = $self->private_make_conservation_profile($parameters, $reg_loc);
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
	push (@conservation_profiles, $conservation_profile);
      }
    }
    if ($stats_required) {
      $exception_info->print_statistics_and_die;
    }
    return @conservation_profiles;
  } # make_conservation_profiles #

  method make_conservation_profiles_all_pairs_of_genomic_intervals (Window_Pair_Algorithm_Parameters $window_pair_parameters, Genomic_Interval_Set $gi_set) {
      my @genomic_intervals = @{$gi_set->genomic_interval_set};
      my $number_of_intervals = @genomic_intervals;
      my %mapping = ();
      my $key;
      my $aggregate_exception = Aggregate_Exception->new();
      my $stats_required = FALSE;
      my $window_length = $self->private_get_window_length($window_pair_parameters);
      my @step_widths = $self->private_get_step_widths($window_pair_parameters);
      for (my $i=0;$i<$number_of_intervals;$i++) {
	  my $start_point = $i+1;
	  for (my $j=$start_point;$j<$number_of_intervals;$j++) {
	      my @profiles;
	      eval {
		  @profiles = $self->private_get_conservation_profile_for_two_genomic_intervals($window_pair_parameters,$genomic_intervals[$i],$genomic_intervals[$j]);
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
		  my @profile_1 = @{$profiles[0]};
		  my @profile_2 = @{$profiles[1]};
		  my $conservation_profile_pair = Conservation_Profile_Pair->new(
		      profile1 => \@profile_1,
		      profile2 => \@profile_2,
		      window_length => $window_length,
		      step_width1 => $step_widths[0],
		      step_width2 => $step_widths[1],
		      );
		  $key = $i."_".$j; 
		  $mapping{$key} = $conservation_profile_pair; # insert Conservation_Profile_Pair into hash, where key is made of indices referring to the GI array
	      }
	  }
      }
      if ($stats_required) {
	  $aggregate_exception->print_statistics_and_die;
      }
      my $conservation_profiles_object = Conservation_Profiles->new(
	  genomic_intervals => \@genomic_intervals,
	  profile_pairs => \%mapping # HashRef
	  );
      return $conservation_profiles_object;
  } # make_conservation_profiles_all_pairs_of_genomic_intervals #

  method private_get_conservation_profile_for_two_genomic_intervals (Window_Pair_Algorithm_Parameters $window_pair_parameters, Genomic_Interval $first_gi, Genomic_Interval $second_gi) {
      my $star_bundler = Star_Bundler->new();
      my @profiles;
      # create and fill in $parameters
      # (these parameters are not actually needed, the get_conservation_profiles-method of Star_Bundler is asking too much, so we make it up here)
      my $kingdom = 'plants';
      my $evolutionary_tree_maker = Evolutionary_Tree_Maker->new();
      my $tree = $evolutionary_tree_maker->make_evolutionary_tree($kingdom);
      my $partial_threshold_matrix_maker = Partial_Threshold_Matrix_Maker->new();
      my $partial_threshold_matrix = $partial_threshold_matrix_maker->make_partial_threshold_matrix($tree, $kingdom);
      my $star_bundler_params = Star_Bundler_Parameters->new(partial_threshold_matrix => $partial_threshold_matrix);
      my $sequence_params = Sequence_Parameters->new(region => 'upstream');
      my @homolog_databases_for_phylogenetic_remos;
      my $parameters = ReMo_Set_Phylogenetic_Constructor_Parameters->new(sequence_databases_to_use_for_homologs => \@homolog_databases_for_phylogenetic_remos,
									 window_pair_algorithm_parameters => $window_pair_parameters,
									 star_bundler_parameters => $star_bundler_params,
									 sequence_parameters => $sequence_params);      
      @profiles = $star_bundler->get_conservation_profiles($parameters, $first_gi, $second_gi);
      return @profiles;
  } # private_get_conservation_profile_for_two_genomic_intervals #

  method private_make_conservation_profile (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, Reg_Loc $reg_loc) {
    
    my @all_gis;
    my @conservation_profile_pairs;
	
    my $gdbu = Genome_DB_Utilities->new();
    # get genomic intervals (method of genome_db_utilities) and pass to Star_Bundler to obtain profiles
    my $centre_gi = $gdbu->get_genomic_sequence($reg_loc, $parameters->sequence_parameters);	
    my $gi_set_maker = Genomic_Interval_Set_Maker->new();
    my $phylo_gi_set = $gi_set_maker->make_gi_set_for_phylogenetic_remos($parameters, $reg_loc);
    
    push (@all_gis, $centre_gi, @{$phylo_gi_set->genomic_interval_set}); # make an array of all the genomic intervals

    my %mapping = ();
    my $key;
    my $aggregate_exception = Aggregate_Exception->new();
    my $stats_required = FALSE;
    my $window_length = $self->private_get_window_length($parameters->window_pair_algorithm_parameters);
    my @step_widths = $self->private_get_step_widths($parameters->window_pair_algorithm_parameters);
    for ( my $j = 1; $j < scalar(@all_gis); $j++ ) {
	my @profiles;
	eval {
	    @profiles = Star_Bundler->new()->get_conservation_profiles($parameters, $centre_gi, $all_gis[$j]);
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
	    my @profile_1 = @{$profiles[0]};
	    my @profile_2 = @{$profiles[1]};
	    my $conservation_profile_pair = Conservation_Profile_Pair->new(
		profile1 => \@profile_1,
		profile2 => \@profile_2,
		window_length => $window_length,
		step_width1 => $step_widths[0],
		step_width2 => $step_widths[1],
		);
	    $key = "0"."_".$j; # centre_gi is in position 0, the comparitor is indexed by $j
	    $mapping{$key} = $conservation_profile_pair; # insert Conservation_Profile_Pair into hash, where key is made of indices referring to the GI array
	}
    }
    if ($stats_required) {
	die $aggregate_exception;
    }
    my $conservation_profiles_object = Conservation_Profiles->new(
	genomic_intervals => \@all_gis,
	profile_pairs => \%mapping # HashRef
	);
    return $conservation_profiles_object;
  } # private_make_conservation_profile #
  
  method private_get_window_length (Window_Pair_Algorithm_Parameters $wpap) {
      my $result;
      if ($wpap->isa('Ott_Algorithm_Parameters')) {
	  $result = $wpap->windowlength;
      } elsif ($wpap->isa('Seaweed_Algorithm_Parameters')) {
	  $result = $wpap->windowlength;
      } else {
	  die 'could not establish window length for this type of window pair algorithm.';
      }
      return $result;
  } # private_get_window_length #

  method private_get_step_widths (Window_Pair_Algorithm_Parameters $wpap) {
      my @result;
      if ($wpap->isa('Ott_Algorithm_Parameters')) {
	  push(@result,$wpap->stepwidth1);
	  push(@result,$wpap->stepwidth2);
      } elsif ($wpap->isa('Seaweed_Algorithm_Parameters')) {
	  push(@result,$wpap->stepwidth);
	  push(@result,$wpap->stepwidth);
      } else {
	  die 'could not establish step width for this type of window pair algorithm.';
      }
      return @result;
  } # private_get_step_widths #

} # Conservation_Profiles_Maker #
