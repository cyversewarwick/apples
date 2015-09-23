### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### General_Utilities Class ###
use MooseX::Declare;

class Bundler::Dependencies::General_Utilities {
	use feature qw(state);
	use Bundler::Dependencies::APPLES_Datatypes qw(Boolean);
	use Bundler::Dependencies::Running_Perseverance;
	use Config::General;
	use Digest::MD5 qw(md5_hex);
	use constant {FALSE => 0,
		      TRUE  => 1};
	use Data::Dumper;
	use File::Spec;
	use File::Temp qw (tempdir);	
	use Bundler::Dependencies::Exception;
	use XML::Simple qw/XMLin/;

	method remove_duplicates_from_list (ArrayRef $list_ref) {
	  # takes an array, removes any duplicate elements, returns array of only the unique elements
	  my @unique_list;
	  my %seen;
	  foreach my $element( @$list_ref ){
	    next if $seen{ $element }++;
	    push (@unique_list, $element);
	  }
	  return @unique_list;
	} # remove_duplicates_from_list #

	method list_has_redundancy (ArrayRef $list_ref) {
	    my $length = @{$list_ref};
	    my @non_redundant = $self->remove_duplicates_from_list($list_ref);
	    my $non_redundant_length = @non_redundant;
	    my $result = TRUE;
	    if ($length == $non_redundant_length) {
		$result = FALSE;
	    }
	    return $result;
	} # list_has_redundancy #
	
	method lists_overlap (ArrayRef $listeref1, ArrayRef $listref2) {
	  # take 2 arrays, determines if any of the list elements overlap, returns TRUE or FALSE
	  my $result = FALSE;
	  foreach my $element1 (@{$listeref1}) {
	    foreach my $element2 (@{$listref2}) {
	      if ($element1 eq $element2 ) {
		$result = TRUE;
	      }
	    } 
	  }
	  return $result;
	} # lists_overlap #

	method list_intersection(ArrayRef $reference_list, ArrayRef $subject_list) {
	    # returns all elements of $subject_list that are also in $reference_list
	    my %hash;
	    foreach my $element (@{$reference_list}) {
		$hash{$element} = 13;
	    }
	    my @result;
	    foreach my $element (@{$subject_list}) {
		if (defined $hash{$element}) {
		    push(@result,$element);
		}
	    }
	    return @result;
	} # list_intersection #

	method subset_array_by_booleans (ArrayRef $array_ref, ArrayRef $booleans) {
	    # both arrays must be same length

	    my @result;
	    my $length = @{$array_ref};
	    my $other_length = @{$booleans};
	    if ($length != $other_length) {
		die 'arrays must be of same length! ('.$length.' vs '.$other_length.')';
	    }
	    for (my $i=0;$i<$length;$i++) {
		if (${$booleans}[$i]) {
		    push(@result,${$array_ref}[$i]);
		}
	    }
	    return @result;
	} # subset_array_by_booleans #

	method user_info (Int $level, Str $text) {
	  # to do: make input parameter a type rather than any integer
	  state $global_level = $main::global_print_level;
	  
	  if (!defined $main::global_print_level || $global_level < 0 || $global_level > 3) {
	    die 'You must set $global_print_level to a value between 0 and 3';
	  }
	  elsif ( $global_level == 0 ) {
	    # no print statements
	    return;
	  }
	  elsif ( $global_level >= $level ) {
	    # print the contents of $text
	    print $text;
	  }
	  return;
	} # user_info #
	
	method timestamp () {
	    # do not change, Scheduler makes use of this particular formatting

	  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time); 
          $year=$year+1900;
	  $mon++;
          return $year . sprintf('%02u%02u_%02u%02u%02u', $mon, $mday, $hour, $min, $sec); 
	} # timestamp #

	method load_includes (Str $APPLES_dat_path) {
	  $self->user_info(2,"LOAD INCLUDES\n");
	  # loads ALL required includes from APPLES.dat
	  # sets APPLES_DAT environmental variable
	  
	  $Storable::canonical = 1; # -> "Storable will store hashes with the elements sorted by their key"
	  # turn relative path into absolute path
	  $APPLES_dat_path = File::Spec->rel2abs($APPLES_dat_path);
	  $APPLES_dat_path =~ tr|\\|/|;

	  my $APPLES_conf = new Config::General($APPLES_dat_path);
	  # read APPLES config file
	  my %APPLES_config = $APPLES_conf->getall();

	  # reset PERL5LIB environment variable
	  $ENV{PERL5LIB}='/common/perl-5.10.0/bin/';
	# set BLAST environmental variables
	$ENV{BLASTDIR} = '/common/bin/'; # should be taken from conf file
	  use lib '/common/perl-5.10.0/bin/'; # testing override of PERL5LIB - this line probably not necessary
	  no lib '/System/Library/Perl/5.8.6/darwin-thread-multi-2level'; # hack to get rid of wsbc cluster environmental variable
	  no lib '/System/Library/Perl/5.8.6'; # hack to get rid of wsbc cluster environmental variable
	  no lib '/System/Library/Perl/Extras/5.8.6/darwin-thread-multi-2level'; # hack to get rid of wsbc cluster environmental variable
	  my $bio_perl = $APPLES_config{bio_perl};
	  my $ensembl = $APPLES_config{ensembl};
	  my $ensembl_compara = $APPLES_config{ensembl_compara};
	  my $apples_main = $APPLES_config{APPLES_main};
	  my $input_dir = $APPLES_config{input_dir};
	  my $non_system_temp_dir = $APPLES_config{non_system_temp_dir};
	  my $data_dumps = $APPLES_config{data_dumps};
	  my $blast_db = $APPLES_config{blast_db};
	  my $path_to_perl = $APPLES_config{path_to_perl};
	  my $ensembl_registry_conf = $APPLES_config{ensembl_registry_conf};
	  push(@INC, $ensembl);
	  push(@INC, $ensembl_compara);
	  push(@INC, $bio_perl);
	  push(@INC, $apples_main);
	  push(@INC, $apples_main."/DataInterface/");
	  push(@INC, $apples_main."/DataInterface/XML");
	  push(@INC, $input_dir);
	  push(@INC, $non_system_temp_dir);
	  push(@INC, $data_dumps);
	  push(@INC, $blast_db);
	  $ENV{'APPLES_DAT'} = $APPLES_dat_path;
	  $ENV{'input_dir'} = $input_dir;
	  $ENV{'non_system_temp_dir'} = $non_system_temp_dir;
	  $ENV{'data_dumps'} = $data_dumps;
	  $ENV{'blast_db'} = $blast_db;
	  $ENV{'path_to_perl'} = $path_to_perl;
	  $ENV{'ensembl_registry_conf'} = $ensembl_registry_conf;
	} # load_includes
	
	method md5sum (Str $filepath) {
	  my $md5_generator = Digest::MD5->new();
	  open(FILE, $filepath) or die "Can't open '$filepath': $!";
	  binmode(FILE);
	  my $result = $md5_generator->addfile(*FILE)->hexdigest or die "Cannot make md5 of $filepath.\n"; # obtain MD5 of file at $filepath
	  close(FILE);
	  return $result;
	} # md5sum #

	method wait_for_file_existence (Str $filename) {	  
	  my $counter = 0;
	  my $limit   = 4;
	  while ( ( !-e $filename )
		 && ( $counter < $limit ) ) {
	      sleep 5;
	      $counter++;
	  }
	  if ( $counter > 0 ) {
	    my $time = 5 * $counter;
	    my $onemessage = "Had to wait for file " . $filename . ": " . $time . " seconds.";
	    $self->user_info(2, $onemessage);  
	  }
	  if ( $counter >= $limit ) {
	   $self->user_info(2, "$filename does not exist");  
	    die 'file $filename does not exist';
	  }
	  else {
	    return 1; # file exists
	  }
	} # wait_for_file_existence #

	method wait_until_one_can_write_to_file (Str $filename) {
	    $self->wait_for_file_existence($filename);
	    my $done = FALSE;
	    while (!$done) {
		if (open(FILEHANDLE,'>>'.$filename)) {
		    $done = TRUE;
		};
		sleep 1;
	    }
	    close FILEHANDLE;
	} # wait_until_one_can_write_to_file #
	
	method get_temp_random_directory(Boolean $use_perseverance) {
	    # creates a temporary directory within the APPLES non-system temporary directory, perseverance is optional

	    my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
	    my %APPLES_config = $APPLES_conf->getall();
	    my $non_system_temp_dir = $APPLES_config{non_system_temp_dir};
	    my $tempdir;
	    if ($use_perseverance) {
		$tempdir = $self->get_temp_dir_and_use_perseverance($non_system_temp_dir);
	    } else {
		$tempdir = tempdir (DIR => $non_system_temp_dir); # prefixes a temp directory with this location
	    }
	    chmod (0777, $tempdir);
	    $tempdir = $tempdir."/";		
	    return $tempdir;
	} # get_temp_random_directory #

	method get_temp_dir_and_use_perseverance (Str $basedir) {
	    # creates a temporary directory within given base directory
	    # can be used as a more resilient replacement for: tempdir ( DIR => $basedir )

	    my $result;
	    while (!$main::global_perseverance->stop_trying) {
		eval {
		    $result = tempdir ( DIR => $basedir );
		};
		my $error = $@;
		$main::global_perseverance->decide_on_rerun('TempDirs', TRUE, $error);
	    }
	    $main::global_perseverance->stop_trying(FALSE);
	    chmod (0777, $result);
	    return $result;
	} # get_temp_dir_and_use_perseverance #

	method maximum(ArrayRef $list_ref) {
	    my @numbers = @$list_ref;

	    my $result = $numbers[0];
	    foreach my $number (@numbers) {
		if ( $number > $result ) {
		    $result = $number;
		}
	    }
	    return $result;
	} # maximum #

	method minimum(ArrayRef $list_ref) {
	    my @numbers = @$list_ref;

	    my $result = $numbers[0];
	    foreach my $number (@numbers) {
		if ( $number < $result ) {
		    $result = $number;
		}
	    }
	    return $result;
	} # minimum #

	method standard_exception_handling_for_parallelisation(Any $exception_content, Aggregate_Exception $aggregate_exception) {
	    my @result;
	    my $stats_required = FALSE;
	    
	    if (!UNIVERSAL::can($exception_content, 'isa')) {
			die $exception_content; # throw error string
	    }
	    else {
			if ($exception_content->isa('Job_Information_Exception')) {
		    	$stats_required = TRUE;
		    	$aggregate_exception->merge($exception_content);
			}
			elsif ($exception_content->isa('Aggregate_Exception')) {
		    	$aggregate_exception->merge_aggregate($exception_content);
		    	$stats_required = TRUE;
			}
			elsif ($exception_content->isa('Exception')) {
		    	$self->user_info(3,"catching exception:\n".$exception_content->message."\n");
		    	die $exception_content;
			}
			else {
		    	$self->user_info(3,"caught non-APPLES exception\n");
		    	die $exception_content; # throw error object
			}
	    }
	    $result[0] = $stats_required;
	    $result[1] = $aggregate_exception;
	    return @result;
	} # standard_exception_handling_for_parallelisation #
	
	method repeat_character(Str $character, Int $n) {
	    my $result = '';
	    for (my $i=0;$i<$n;$i++) {
		$result = $result.$character;
	    }
	    return $result;
	} # repeat_character #

	method repeat_array_element(Any $element, Int $number) {
	    # returns an array containing $number copies of $element

	    my @result = ();
	    if ($number>0) {
		for (my $i=0;$i<$number;$i++) {
		    push(@result,$element);
		}
	    }
	    return @result;
	} # repeat_array_element #

	method standard_call_to_job_handler(Any $object, Str $function, ArrayRef $parameters, Boolean $high_memory, Boolean $cache_in_memory) {
	  # this call to the Job_Handler will not return the (possibly altered) $object, but only the results directly returned by $function
          # (in most cases the object is not needed)

	    my $cache = TRUE;
	    my $cache_duration = 180;
	    my $job_handler = Job_Handler->new();
	    $job_handler->get_config_settings();
	    my $job_parameters = Job_Parameters->new(memory_requirement_high => $high_memory,
						     wall_time_estimate => 172800,
						     cache_in_memory => $cache_in_memory
		);
	    my $aggregate_exception = Aggregate_Exception->new();
	    my @cache_result = eval {
		$job_handler->handle_APPLES_function($function, $object, $parameters, $cache, $cache_duration, $job_parameters);
	    };
	    if ($@) {
		my $exception_content = $@;
		my @outcome = $self->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
		  # no parallelism here, though
		$aggregate_exception = $outcome[1];
		if ($outcome[0]) {
		    $aggregate_exception->print_statistics_and_die;
		} else {
		    die 'did not expect this could happen at design time of this code.';
		}
	    }
	    shift @cache_result; # throwing away the object
	    return @cache_result;
	} # standard_call_to_job_handler #

	method standard_call_to_job_handler_without_exception_handling(Any $object, Str $function, ArrayRef $parameters, Boolean $high_memory, Boolean $cache_in_memory) {
	  # this call to the Job_Handler will not return the (possibly altered) $object, but only the results directly returned by $function
          # (in most cases the object is not needed)

	    my $cache = TRUE;
	    my $cache_duration = 180;
	    my $job_handler = Job_Handler->new();
	    $job_handler->get_config_settings();
	    my $job_parameters = Job_Parameters->new(memory_requirement_high => $high_memory,
						     wall_time_estimate => 172800,
						     cache_in_memory => $cache_in_memory
		);
	   
	    my @cache_result = $job_handler->handle_APPLES_function($function, $object, $parameters, $cache, $cache_duration, $job_parameters);
	    shift @cache_result; # throwing away the object
	    return @cache_result;
	} #  standard_call_to_job_handler_without_exception_handling #

	method append_string_to_a_file(Str $filename, Str $outputstring) {
	    open FILETOAPPENDTO, '>>'.$filename;
	    print FILETOAPPENDTO $outputstring;
	    close FILETOAPPENDTO;
	} # append_string_to_a_file #

	method count_Ns(Str $sequence) {
	    my $upper_case = uc($sequence);
	    my $result = ($upper_case =~ tr/N//);
	    return $result;
	} # count_Ns #

	method write_items_to_file(Str $full_file_name, Boolean $append, ArrayRef $items) {
	    # $items must contain elements that can be printed by the print-function

	    my $prefix = '>';
	    if ($append) {
		$prefix = '>>';
	    }
	    open TARGET_FILE, $prefix.$full_file_name or die "Cannot write to file\n";
	    foreach my $item (@{$items}) {
		print TARGET_FILE $item."\n";
	    }
	    close TARGET_FILE;
	} # write_items_to_file #

	method read_in_web_service_parameters(Str $parameters) {
	 	# parse web service parameters from XML input
		my $ref = XMLin($parameters, NoAttr => 1);
		#PEB added NoAttr as I have added attributes to xml to help web interface
		#not relevant to tool and would complicate parsing if present
	    my %result;
	    ## loop through command line parameters
	    while ( my ($name, $value) = each (%{$ref}) ) {
			# in XML, labeling parameters with PARAMx_... is not necessary
			# we just ignore the param bit here. 
			
			$name =~ s/^PARAM[0-9]+\_//;
				
			# translate booleans
		    if ($value =~ m/true/i) {
				$value = 1;
		    }
		    if ($value =~ m/false/i) {
				$value = 0;
		    }
		
			# renaming to capital letters. this could be 
			# removed if web services wrote OUTPUTDIR directly
			if($name eq 'output_tmp_dir') {
				$name = 'OUTPUTDIR';
			}

			if(ref($value) eq 'ARRAY') {
				push @{$result{$name}}, $_ foreach (@$value);				
			} else {
				push @{$result{$name}}, $value;
			}
	    }
		
	    ## check and output directory was specified
	    my $output_dir;
	    if (!defined $result{'OUTPUTDIR'}) {
		    my $web_service_exception = Web_Service_Parameters_Inconsistent->new();
			die $web_service_exception;
	    } else {
			my @directories = @{$result{'OUTPUTDIR'}};
			if ($#directories < 0) {
			    my $web_service_exception = Web_Service_Parameters_Inconsistent->new();
			    die $web_service_exception;
			}
			$output_dir = $directories[0];
	    }
	
	    ## write parameters output file
	    my $output_filename = $output_dir.'/parameters.txt';
	    open PARAMETERS_OUTPUT, '>'.$output_filename or die "System error!";
	    print PARAMETERS_OUTPUT "Web service parameters received:\n";
	    print PARAMETERS_OUTPUT Dumper($parameters);
	    print PARAMETERS_OUTPUT "\n\nParameter data structure after parsing command line:\n";
	    print PARAMETERS_OUTPUT Dumper(\%result);
	    print PARAMETERS_OUTPUT "\n";
	    close PARAMETERS_OUTPUT;
	    ## return result
	    return \%result;
	} # read_in_web_service_parameters #

} # General_Utilities #
