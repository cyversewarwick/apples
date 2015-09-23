### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Window_Pair_Alignment Class ###
# class to wrap calls to Window Pair Algorithm binaries

use MooseX::Declare;

class Window_Pair_Alignment {
  use File::Temp qw (tempdir);
  use File::Spec;
  use Cwd;
  use Genomic_Interval;
  use ReMo_Data;
  use Data::Dumper;
  use General_Utilities;

  my $GU = General_Utilities->new();

  method run(Genomic_Interval $first_sequence, Genomic_Interval $second_sequence, Str $tempdir, Str $binary_path, Window_Pair_Algorithm_Parameters $parameters) {
  
    my $stepwidth1;
    my $stepwidth2;
    
    my $command = $binary_path; # full path and binary filename is specified in Job_Handler_Config.dat

    if ($parameters->isa('Ott_Algorithm_Parameters')) {
      $stepwidth1 = $parameters->stepwidth1;
      $stepwidth2 = $parameters->stepwidth2;
    }
    elsif ($parameters->isa('Seaweed_Algorithm_Parameters')) {
      $stepwidth1 = $parameters->stepwidth;
      $stepwidth2 = $stepwidth1;
    }
    else {
      die 'unknown algorithm parameter type';
    }
    my $windowlength = $parameters->windowlength;
    my $method_threshold = $parameters->cutoff_for_uninteresting_alignments;
    
    my $file1 = $tempdir."/".'tempfile1'; # need to save first sequence to this file and put in stepwidth1 stepwidth2 windowlength and cutoff values in first line
    open FILE1, '>', $file1 or die $!;
    print FILE1 $stepwidth1,"\t",$stepwidth2,"\t",$windowlength,"\t",$method_threshold,"\t\n",$first_sequence->get_sequence(),"\n"; # newline is essential to get correct profile length from Ott algorithm!
    close FILE1;
    
    my $file2 = $tempdir."/"."tempfile2";
    open FILE2, '>', $file2 or die $!;
    print FILE2 $second_sequence->get_sequence(),"\n"; # newline is essential to get correct profile length from Ott algorithm!
    close FILE2;
    my $resultdir = $tempdir."/";
    my $profile1 = "./profile1";
    my $profile2 = "./profile2";
    my $result = "./result.txt";

    my $commandline = $command. " $file1 $file2 $result $profile1 $profile2"."\n";

    $GU->user_info (3, "\nCommand to run:\n$commandline\n");
    return $commandline;
  } # run #

  method get_result(Str $tempdir, Window_Pair_Algorithm_Parameters $parameters) { 
    my @alignment_result;
    my $file;
    # first get the windowpair_result
    $file = "/result.txt";
    my $data = $self->private_get_raw_data($tempdir, $file);
    my @windowpair_result = $self->private_make_windowpair_result ($data, $parameters);
    # now get the first profile
    $file = "/profile1"; # IS NAME SAME FOR EACH BINARY? (seaweed=profilea ?)
    $data = $self->private_get_raw_data($tempdir, $file);
    my @profile1 = $self->private_make_profile_result($data);
    # now get the second profile
    $file = "/profile2"; # IS NAME SAME FOR EACH BINARY? (seaweed=profileb ?)
    $data = $self->private_get_raw_data($tempdir, $file);
    my @profile2 = $self->private_make_profile_result($data);

    # push all 3 results into the alignment_result array
    @alignment_result = (\@windowpair_result, \@profile1, \@profile2);
    
    # return the alignment_result array
    return @alignment_result;#return result (arrayrefs to windowpairs+profile1+profile2)
  } # get_result #

  method private_get_raw_data(Str $tempdir, Str $file) {  
    my $resultsfile = $tempdir . $file;
    my $fh;
    $GU->wait_for_file_existence ($resultsfile); # will wait up to 20 seconds for file to exist
    open $fh,$resultsfile or die "Can't open $resultsfile";
    my $data = do {local $/; <$fh>};	# slurp in the scores
    close $fh;
    return $data; # string
  } # private_get_raw_data #
  
  method private_make_windowpair_result (Any $data, Window_Pair_Algorithm_Parameters $parameters) { # $data - type?
    
    my @regions = split "\n",$data;
    shift @regions; # get rid of first line - just indicates the cutoff limit
   
    # get the windowpair results:
    my @windowpairs;
    foreach my $line (@regions) {
      my @values = split(/\t/, $line);
      
      if ($values[2] >= $parameters->cutoff_for_uninteresting_alignments) {
	my $windowpair = WindowPair->new(offset1=>$values[0], offset2=>$values[1], score=>$values[2]);#join "," , ($start, $start + $self->window-1, $score);
	push @windowpairs,$windowpair;
      }
    }
    
    return @windowpairs;

  } # private_make_windowpair_result #
  
  method private_make_profile_result(Any $data) { # datatype = Str?
    my @profile = split "\n",$data;
    # NB last element of Seaweed profile is a zero, not a windowpair score! But profile length is correct. Bug reported for WindowAlignment_BSP_darwin_default_release binary.
    return @profile;
  } # private_make_profile_result #
  
} # Window_Pair_Alignment #

class Ott_Window_Pair_Alignment extends Window_Pair_Alignment {

} # Ott_Window_Pair_Alignment #

class Seaweed_Window_Pair_Alignment extends Window_Pair_Alignment { 

} # Seaweed_Window_Pair_Alignment #






 


