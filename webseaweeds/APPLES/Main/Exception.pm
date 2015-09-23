### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Exception Class ###
use MooseX::Declare;

class Exception {
  has 'message' => (is => 'rw', isa => 'Str');
} # Exception #

class Job_Information_Exception extends Exception {
  has 'job_type' => (is => 'rw', isa => 'Str', required => 1);
  has 'CPU_time_estimate' => (is => 'ro', isa => 'Num');
  has 'job_info' => (is => 'ro', isa => 'HashRef');
} # Job_Information_Exception #

class Aggregate_Exception extends Exception {
  use Data::Dumper;
  use APPLES_Datatypes qw (NonNegativeNum NonNegativeInt);
  use General_Utilities;

  use Moose;
  has 'messages' =>  (is => 'rw', isa => 'ArrayRef[Str]');
  has 'total_CPU_time' => (is => 'rw', isa => 'HashRef[NonNegativeNum]');
  has 'counters' => (is => 'rw', isa => 'HashRef[NonNegativeInt]');
  has 'total_job_info' => (is => 'rw', isa => 'ArrayRef[Any]'); # strictly should be ArrayRef[HashRef]

  my $GU = General_Utilities->new();

  method merge (Job_Information_Exception $job_info_exception) {
    # message
    if (!defined $self->messages) {
      my @empty_array=();
      $self->messages(\@empty_array); # so that the dereferencing of an undefined arrayref below does not cause a crash
    }
    my @messages = @{$self->messages};
    push(@messages, $job_info_exception->message);
    $self->messages(\@messages);# syntax was $self->messages =\@messages;
    # CPU-time
    if (!defined $self->total_CPU_time) {
      my %empty_hash = ();
      $self->total_CPU_time(\%empty_hash);
      # so that the dereferencing of an undefined hashref below does not crash
    }
    if (defined $job_info_exception->CPU_time_estimate) {
      my $current_total = ${$self->total_CPU_time}{$job_info_exception->job_type};
      if (!defined $current_total) {
	$current_total = 0;
      }
      my $new_total = $current_total+$job_info_exception->CPU_time_estimate;
      ${$self->total_CPU_time}{$job_info_exception->job_type} = $new_total;
      
    }
    # count
    if (!defined $self->counters) {
      my %empty_hash = (); 
      $self->counters(\%empty_hash); # so that the dereferencing of an undefined hashref below does not cause a crash
    }
    my $current_count = ${$self->counters}{$job_info_exception->job_type};
    if (!defined $current_count) {
      $current_count = 0;
    }
    my $newcount = $current_count+1;
    ${$self->counters}{$job_info_exception->job_type} = $newcount;

    # job info
    if (!defined $self->total_job_info) {
      my @empty_array = ();
      $self->total_job_info(\@empty_array); # so that dereferencing of an undefined arrayref below does not cause a crash
    }
    my @total_job_info = @{$self->total_job_info};
    push (@total_job_info, $job_info_exception->job_info);
    $self->total_job_info(\@total_job_info);
  } # merge #
	
  method merge_aggregate (Aggregate_Exception $exception) {

    # message
    if (!defined $self->messages) {
      my @empty_array = ();
      $self->messages(\@empty_array); # so that the dereferencing of an undefined arrayref below does not cause a crash
    }
    my @onearray = @{$self->messages};
    my @otherarray = @{$exception->messages};
    my @joint_array = (@onearray,@otherarray);
    $self->messages(\@joint_array);

    # CPU-time
    if (!defined $self->total_CPU_time) {
      my %empty_hash = ();
      $self->total_CPU_time(\%empty_hash); # so that the dereferencing of an undefined hashref below does not cause a crash
    }
    my %one_hash = %{$exception->total_CPU_time};
    my %other_hash = %{$self->total_CPU_time};
    my @one_set_of_keys = keys %one_hash;
    my @other_set_of_keys = keys %other_hash;
    my @all_keys = (@one_set_of_keys,@other_set_of_keys);
    my %joint_hash;
    foreach my $key (@all_keys) {
      $GU->user_info(3,"key: ".$key."\n");
      my $one = 0;
      my $other = 0;
      if (exists($one_hash{$key})) {
	$one = $one_hash{$key};
      }
      if (exists($other_hash{$key})) {
	$other = $other_hash{$key};
      }
      #$joint_hash{$key} = $one_hash{$key}+$other_hash{$key};
      $joint_hash{$key} = $one + $other;
    }
    $self->total_CPU_time(\%joint_hash);
    # count
    if (!defined $self->counters) {
      my %empty_hash = ();
      $self->counters(\%empty_hash); # so that the dereferencing of an undefined hashref below does not cause a crash
    }
    my %one_count_hash = %{$exception->counters};
    my %other_count_hash = %{$self->counters};
    my @one_set_count_keys = keys %one_count_hash;
    my @other_set_count_keys = keys %other_count_hash;
    my @all_count_keys = (@one_set_count_keys, @other_set_count_keys);
    my %joint_count_hash;
    foreach my $count_key (@all_count_keys) {
      my $one_count = 0;
      my $other_count = 0;
      if (exists($one_count_hash{$count_key})) {
	$one_count = $one_count_hash{$count_key};
      }
      if (exists($other_count_hash{$count_key})) {
	$other_count = $other_count_hash{$count_key};
      }
      $joint_count_hash{$count_key} = $one_count + $other_count;
    }
    $self->counters(\%joint_count_hash);
    # job info
    if (!defined $self->total_job_info) {
      my @empty_array = ();
      $self->total_job_info(\@empty_array); # so that the dereferencing of an undefined arrayref below does not cause a crash
    }
    my @one_info_array = @{$self->total_job_info};
    my @other_info_array = @{$exception->total_job_info};
    my @joint_info_array = (@one_info_array,@other_info_array);
    $self->total_job_info(\@joint_info_array);
  } # merge_aggregate #

  method print_statistics_and_die () {

    my @count_keys = keys %{$self->counters};
    my @CPU_time_keys = keys %{$self->total_CPU_time};
    my @joint_keys = (@count_keys,@CPU_time_keys);

    my @non_dup = $GU->remove_duplicates_from_list(\@joint_keys);

    my $global_mode = $main::APPLES_running_mode;
    $GU->user_info(1, "\n\n****** Job Statistics: ******\n");
    $GU->user_info(1, "** (running mode: ".$global_mode.")\n\n");
    foreach my $key (@non_dup) {
      my $count = ${$self->counters}{$key};
      my $CPU_time = ${$self->total_CPU_time}{$key};
      my $jobstring = 'jobs';
      if ($count == 1) {
	  $jobstring = 'job';
      }
      my $line = $key.":\t".$count." ".$jobstring.",\ttotal of ".$CPU_time." seconds estimated.\n";
      $GU->user_info(1,$line);
    }
    $GU->user_info(1, "** End of Job Statistics: ***\n\n\n");
    die 'Cannot proceed further.';
  } # print_statistics_and_die #

  method print_job_info_to_file (Str $filename) {
  	
  	
  	my $APPLES_DAT = $ENV{'APPLES_DAT'};
    my $APPLES_conf = new Config::General($APPLES_DAT);
    my %APPLES_config = $APPLES_conf->getall();
    my $queue = $APPLES_config{queue};

    $GU->user_info(1, "printing job info to file\n");
    # prints job info (MD5sum and queue directory name for each submitted job) to a user-specified file
    #print "\n\n\n\nprinting job info to file $filename \n\n\n";
    my @array_of_hashrefs = @{$self->total_job_info};
   # print Dumper (\@array_of_hashrefs);
    foreach my $entry (@array_of_hashrefs) {
      my %hash = %{$entry};
      my $MD5 = $hash{MD5sum};
      my $tmpdir = $hash{tempdir};
      my $method = $hash{method};
      open (MYFILE, ">>$filename") or die "Can't open $filename: $!";
      print MYFILE $queue."/".$tmpdir."\t".$MD5."\t".$method."\n";
      close (MYFILE); 
      $GU->user_info(1, "Written job information to file $filename\n");
    }
  } # print_job_info_to_file
} #  Aggregate_Exception #
  
class No_Neighbouring_Gene_Exception extends Exception {

} # No_Neighbouring_Gene_Exception #

class Web_Service_Parameters_Inconsistent extends Exception {

} # Web_Service_Parameters_Inconsistent #
