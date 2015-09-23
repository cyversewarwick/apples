### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Conservation_Profiles Class ###
# class to contain window pair alignment conservation profiles

use MooseX::Declare;

class Conservation_Profile_Pair {
  # technical class to contain a pair of conservation profiles

    use APPLES_Datatypes qw (PositiveInt);

    has 'profile1' => (is => 'rw', isa => 'ArrayRef[Num]', required => 1); # profile1
    has 'profile2' => (is => 'rw', isa => 'ArrayRef[Num]', required => 1); # profile2
    has 'window_length' => (is => 'ro', isa => PositiveInt, required => 1);
    has 'step_width1' => (is => 'ro', isa => PositiveInt, required => 1);
    has 'step_width2' => (is => 'ro', isa => PositiveInt, required => 1);

} # Conservation_Profile_Pair #

class Conservation_Profiles {
  use Chart::Gnuplot;
  use Data::Dumper;
  use Genome_DB_Utilities;
  use APPLES_Datatypes qw (IntPercentage PositiveInt);
  use File::Temp qw (tempdir);
  use File::Path qw (rmtree);
  use constant {FALSE => 0, TRUE => 1};

  has 'genomic_intervals' => (is => 'rw', isa => 'ArrayRef[Genomic_Interval]', required => 1); # arrayref of genomic intervals
  has 'profile_pairs' => (is => 'rw', isa => 'HashRef', required => 1); # holds a hash of pairs of conservation profiles, the keys being made from index pairs of the genomic_intervals array

  my $GU = General_Utilities->new();
  my $gdbu = Genome_DB_Utilities->new();
  
  method render_profiles_to_graph (Int $gi_index, Str $directory_name, IntPercentage $max_repeat_percentage) {
      # renders the conservation profile(s) available, for the GI index requested (first index 0), as a graph
      # all windows which exceed the maximum repeat percentage will be shown as scoring 0
      
      $GU->user_info(1, "rendering graph...\n");
      # get the data     
      my @datasets;     
      my $key;
      my $profilepair_exists;
      my $profile;
      my $window_length = undef;
      my $next_window_length;
      my $step_width = undef;
      my $next_step_width;
      for (my $j = 0; $j < scalar(@{$self->genomic_intervals}); $j++) {
	  if ($gi_index != $j) {
	      # do profiles exist for this pair?
	      $profilepair_exists = $self->private_profile_exists($gi_index, $j);
	      # get profiles for this key, if it exists
	      if ($profilepair_exists) {
		  $profile = $self->private_obtain_one_profile($gi_index, $j);
		  # check window length
		  $next_window_length = $self->private_obtain_window_length($gi_index, $j);
		  if (defined $window_length) {
		      if ($window_length != $next_window_length) {
			  die 'cannot plot profiles for different window lengths in one plot.';
		      }
		  }
		  $window_length = $next_window_length;
		  # check step width
		  $next_step_width = $self->private_obtain_step_width($gi_index, $j);
		  if (defined $step_width) {
		      if ($step_width != $next_step_width) {
			  die 'plotting of profiles for varying step widths not implemented yet.';
		      }
		  }
		  $step_width = $next_step_width;
		  # create data set
		  my $dataset = $self->private_create_dataset($profile, $gi_index, $j, $max_repeat_percentage, $window_length, $step_width);
		  push (@datasets, $dataset);
	      }
	  }
      }
      # create chart object
      if (defined $window_length) {
	  my $query_name = ${$self->genomic_intervals}[$gi_index]->label;
	  my $outputname = $query_name."_conservation_profiles.png";
	  my $title = "Conservation Profile for ".$query_name." and Orthologues";
	  my $chart = Chart::Gnuplot->new(
	      output => $directory_name.'/'.$outputname,
	      title  => $title,
	      xlabel => "Position Along Sequence (bp)",
	      ylabel => "Alignment Score",
	      yrange => [0, $window_length],
	      );
	  # Plot the data set onto the chart      
	  $chart->plot2d(@datasets);
	  $GU->user_info(1, "done\n");
      } else {
	  $GU->user_info(1, "plot not produced - there was no data.\n");
      }
      return;
  } # render_profiles_to_graph #

  method private_get_hash_key (Int $i, Int $j) {
    my $key = $i."_".$j;
    if ($j < $i) {
      $key = $j."_".$i;
    }
    return $key;
  } # private_get_hash_key #

  method private_profile_exists(Int $i, Int $j) {
    # returns boolean stating existence of a profile pair for a pair of GI indices

    my $key = $self->private_get_hash_key($i, $j);
      
    if ($i == $j) {
      # won't have profiles of an interval against itself
      return FALSE;
    }
    my %hash = %{$self->profile_pairs};
    if (exists $hash{$key}) {
      return TRUE;
    }
    else {
      return FALSE;
    }
  } # private_profile_exists #

  method private_obtain_one_profile (Int $i, Int $j) {
      # obtains profile for $i relative to $j

      # get key
      my $key = $self->private_get_hash_key($i, $j);
      # get profile
      my $profile_pair = ${$self->profile_pairs}{$key};
      my @result;
      if ($i<$j) {
	  @result = @{$profile_pair->profile1};
      } else {
	  @result = @{$profile_pair->profile2};
      }
      return \@result;
  } # private_obtain_one_profile #

  method private_obtain_step_width (Int $i, Int $j) {
      # gets the step width used for $i when compared with $j

      my $key = $self->private_get_hash_key($i, $j);
      my $profile_pair = ${$self->profile_pairs}{$key};
      my $result;
      if ($i<$j) {
	  $result = $profile_pair->step_width1;
      } else {
	  $result = $profile_pair->step_width2;
      }
      return $result;
  } # private_obtain_step_width #

  method private_obtain_window_length (Int $i, Int $j) {
      my $key = $self->private_get_hash_key($i, $j);
      my $profile_pair = ${$self->profile_pairs}{$key};
      my $result = $profile_pair->window_length;
      return $result;
  } # private_obtain_window_length #

  method private_create_dataset (ArrayRef[Num] $profile,Int $main_index, Int $j,IntPercentage $max_repeat_percentage,PositiveInt $window_length,PositiveInt $step_width) {
      my $interval = @{$self->genomic_intervals}[$j];
      # work out label
      my $label;
      if (defined $interval->label) {  
	  $label = $interval->label;
      } else {
	  $label = 'Sequence '.$j;
      }
      my $species = $gdbu->get_display_name_for_database($interval->genome_db);
      my $label_and_species = $label." ".$species;
      # take out repetitive windows
      my $main_interval = @{$self->genomic_intervals}[$main_index];
      my $plot_profile = $self->private_set_repetitive_windows_to_zero($profile, $main_interval, $max_repeat_percentage, $window_length, $step_width);
      # create data set
      my $profile_dataset = Chart::Gnuplot::DataSet->new(
	  ydata => $plot_profile,
	  title => $label_and_species,
	  style => "points",
	  ); 
      return $profile_dataset;
  } # private_create_dataset #
  
  method private_set_repetitive_windows_to_zero (ArrayRef[Num] $profile,Genomic_Interval $gi,IntPercentage $max_repeat_percentage,PositiveInt $window_length,PositiveInt $step_width) {

      my @result;
      my $length = @{$profile};
      my $gi_length = $gi->get_length();
      if ($gi_length >= 10000000) {
	  die 'must adapt C-code to deal with a sequence of this length';
      }
      my $sequence = uc($gi->get_sequence());
      my $masked_sequence = uc($gi->get_repeatmasked_sequence());      
      my $tempdir = $GU->get_temp_random_directory(FALSE);
      $GU->write_items_to_file($tempdir.'max_repeat_percentage',FALSE,[$max_repeat_percentage]);
      $GU->write_items_to_file($tempdir.'window_length',FALSE,[$window_length]);
      $GU->write_items_to_file($tempdir.'step_width',FALSE,[$step_width]);
      $GU->write_items_to_file($tempdir.'sequence',FALSE,[$sequence]);
      $GU->write_items_to_file($tempdir.'masked_sequence',FALSE,[$masked_sequence]);
      $GU->write_items_to_file($tempdir.'profile_length',FALSE,[$length]);
      $GU->write_items_to_file($tempdir.'input_profile',FALSE,$profile);
      my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
      my %APPLES_config = $APPLES_conf->getall();
      my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
      my $path_to_executable = $APPLES_C_binaries."/set_repetitive_windows_to_zero";
      system($path_to_executable.' '.$tempdir.'max_repeat_percentage '.$tempdir.'window_length '.$tempdir.'step_width '.$tempdir.'sequence '.$tempdir.'masked_sequence '.$tempdir.'profile_length '.$tempdir.'input_profile '.$tempdir.'profile') == 0 || die "System error!";
      open RESULT_FILE, $tempdir.'profile';
      my $line = <RESULT_FILE>;
      chomp($line);
      my @result_profile = split(/:/, $line);
      close RESULT_FILE;
      rmtree($tempdir);
      my $result_length = @result_profile;
      if ($length != $result_length) {
	  die 'An error must have occurred: the profile has changed its length.';
      }
      return \@result_profile;
  } # private_set_repetitive_windows_to_zero #

} # Conservation_Profiles #

