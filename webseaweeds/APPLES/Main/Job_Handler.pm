### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Job_Handler Class ###
# The Job_Handler provides an interface between APPLES code and the results databases
# A 'job' is defined as any computation that can be run in parallel - e.g. C binaries or 
# APPLES functions wrapped as Perl scripts
# The Job_Handler is the ONLY class that should use the $running_mode global variable 
use MooseX::Declare;

class Job_Handler {
  use APPLES_Datatypes qw(Boolean AlignmentAlgorithm);
  use Data::Dumper;
  use Exception;
  use General_Utilities;
  use Cache;
  use constant {FALSE => 0,
		  TRUE	=> 1};
  use File::Temp qw (tempdir);
  use File::Path qw (rmtree);
  use File::Basename;
  use Window_Pair_Alignment;
  use MEME_wrapper;
  use Scheduler;
  use Storable qw (nstore retrieve);
  use Digest::MD5 qw(md5_hex);# md5_hex($data) returns the digest in hexadecimal form. The length of the returned string is 32 and only contains characters from this set: '0'..'9' and 'a'..'f'.
  use Cwd;

  use Moose;
  has 'APPLES_settings' => (is => 'rw', isa => 'HashRef');
  has 'job_handler_settings' => (is => 'rw', isa => 'HashRef'); # software settings (binaries etc)
  has 'cache_settings' => (is => 'rw', isa => 'HashRef');

  my $GU = General_Utilities->new();
  my $md5_generator = Digest::MD5->new();

  method get_config_settings () {
    # we should add some 'or die' statements in here, in case the *.dat files are not found
    $GU->user_info (3, "getting job handler config settings\n");
    use Config::General;
    my $APPLES_DAT = $ENV{'APPLES_DAT'};
    my $APPLES_conf = new Config::General($APPLES_DAT);
    # read APPLES config file
    my %APPLES_config = $APPLES_conf->getall();
    $self->APPLES_settings(\%APPLES_config);
    # determine job handler config file from APPLES config file
    my $job_handler_config_file = $APPLES_config{job_handler_config};
    # read job handler config file
    my $job_handler_conf = new Config::General($job_handler_config_file);
    # set config parameters, e.g. location of binaries
    my %job_handler_parameters = $job_handler_conf->getall();
    $self->job_handler_settings(\%job_handler_parameters);
    # read cache config file
    my $cache_config_file = $APPLES_config{cache_config};
    my $cache_conf = new Config::General($cache_config_file);
    my %cache_settings = $cache_conf->getall();
    $self->cache_settings(\%cache_settings);
    # confirm correct setting of global variable
    my @global_mode = ($main::APPLES_running_mode);
    my @allowed_modes = ('statistics','statistics_and_retrieval','preparation','normal');
    my $test_result = $GU->lists_overlap(\@global_mode,\@allowed_modes);
    if (!$test_result) {
	die 'You must have mistyped the global running mode.';
    }
  } # get_config_settings #
  
  method handle_alignment_job (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, Genomic_Interval $first_gi, Genomic_Interval $second_gi, AlignmentAlgorithm $alignment_algorithm, Boolean $cache_result, Int $cache_duration) {  
    # handler for Ott and Seaweed algorithms
    $GU->user_info(3, "handling $alignment_algorithm job\n");
    $self->get_config_settings();
    my $global_mode = $main::APPLES_running_mode;
    
    my %keyfile_contents = ();
    my $relevant_parameters = $parameters->window_pair_algorithm_parameters; # should seaweed windowlength NOT be used as part of keyfile?
    $keyfile_contents{'parameters'} = $relevant_parameters;
    $keyfile_contents{'first_gi'} = $first_gi->get_working_sequence();#$first_gi; altered because of caching consistency
    $keyfile_contents{'second_gi'} = $second_gi->get_working_sequence();#$second_gi; altered because of caching consistency

    $GU->user_info( 3, "global running mode is $global_mode\n");

    my $window_pair_algorithm_parameters = $parameters->window_pair_algorithm_parameters;
    my $inputdir = $ENV{'input_dir'};
    my $tempdir = tempdir ( DIR => $inputdir ) or die "Cannot create job folder .\n";
    chmod (0777, $tempdir);
    my $method;
    my $binary;
    if ($alignment_algorithm eq 'ott') {
      $method = 'APPLES_Ott_Job';
      $binary = ${$self->job_handler_settings}{ott_binary};
    }
    elsif ($alignment_algorithm eq 'seaweed') {
      $method = 'APPLES_Seaweed_Job';
      $binary = ${$self->job_handler_settings}{seaweed_binary};
      $self->private_check_seaweeds_job_size ($first_gi, $second_gi);
    }
    my $keyfilename = 'keyfile.txt';
    
    $Storable::canonical = 1; # -> "Storable will store hashes with the elements sorted by their key"
    nstore (\%keyfile_contents, "$tempdir/$keyfilename") or die "Cannot store keyfile contents to file $tempdir/$keyfilename\n"; # serialise RELEVANT parameters to parameters file
    my $filepath = "$tempdir/$keyfilename";
    my $keyfileMD5 = $GU->md5sum($filepath);
    my $username = ${$self->APPLES_settings}{username};
    
    # check cache for result
    my $result = $self->private_check_cache_for_result($keyfileMD5, $method, $username);
    
    my $job_result;
    
    my $CPU_estimate = $self->private_cpu_estimate_alignment($first_gi->get_working_sequence(), $second_gi->get_working_sequence(), $alignment_algorithm);
    if ($result eq 'NULL') { # result doesn't exist in cache
	if ($global_mode eq 'statistics') {
	    rmtree($tempdir);
	    my $exception = Job_Information_Exception->new(message => "$alignment_algorithm job needs computing", job_type => $method, CPU_time_estimate => $CPU_estimate);
	    die $exception;
	}
	elsif ($global_mode eq 'statistics_and_retrieval') {
	    rmtree($tempdir);
	    my $exception = Job_Information_Exception->new(message => "$alignment_algorithm job needs computing", job_type => $method, CPU_time_estimate => $CPU_estimate);
	    die $exception;
	}
	my $job_mode;
	if ($global_mode eq 'preparation') {
	    $job_mode = 'queue';
	}
	elsif ($global_mode eq 'normal') {
	    $job_mode = 'session';
	}
	else {
	    die 'unknown global mode';
	}
	# V2:check queue(s) - the job might already have been submitted (if so, don't submit it again)
	# V2:check cache again (result might have got placed in cache while you were looking in queue(s))
	
	# result not available in either case (V1:queued jobs notwithstanding), so will firstly submit job to Scheduler, 'session' mode will wait until it is available, and then return result, 'queue' mode will exit after job submission.
	my $wrapper = Window_Pair_Alignment->new();
	
	my $command_line = $wrapper->run($first_gi, $second_gi, $tempdir, $binary, $window_pair_algorithm_parameters); # returns commandline for Scheduler
	# prepare files etc for Scheduler input
	# Str $command_line, Str $method, Str $parameters_file, Str $calling_object_file, Str $result_file, Str $user_id, JobMode $job_mode, Boolean $cache_results, Int $cache_duration
	$GU->user_info(3,"command line to be run: ".$command_line."\n");
	my $result_dir = $tempdir; 
	
	# create Job_Parameters object
	my $wall_time_estimate = $CPU_estimate; 
	my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
						 wall_time_estimate => $wall_time_estimate);
	# create + call Scheduler
	my $cache = Cache->new(cache_settings => $self->cache_settings);
	my $scheduler = Scheduler->new(cache=>$cache);
	$scheduler->get_config_settings();
	
	$job_result = $scheduler->queue_job($command_line, $method, $keyfileMD5, $result_dir, $username, $job_mode, $cache_result, $cache_duration, $job_parameters); #returns result - as a (temp) directory, so get result from here later. V2: add cpu estimates to 'batching'?
	$GU->user_info(3,$job_result."\n");
	my $job_info_record = {
			   MD5sum  => $keyfileMD5,
			   tempdir => $job_result,
			   method => $method
			  };
	if ($job_mode eq 'queue') {
	    $GU->user_info (3,  "result not currently available, appropriate job has been submitted\n");
	    # result not available, so will submit job to Scheduler, and exit via a job_information_exception, providing the job_number
	    my $exception = Job_Information_Exception->new(message=>"result not currently available, $alignment_algorithm job has been submitted",
							   job_type => $method, CPU_time_estimate => $CPU_estimate, job_info => $job_info_record); 
	    die $exception;
	}
	elsif ($job_mode eq 'session') {
		# in session mode, the result should have been run within Scheduler.pm 
		# and now be available in the in the cache database. 
		# hence, we check again and have the code below pick things up. 
		# (peterkrusche)
		$result = $self->private_check_cache_for_result($keyfileMD5, $method, $username);
	}	
    }
    # we check again since in session mode, this variable might have been updated 
    if($result) { # result exists!
	# $result is the result unzipped into a temporary directory
	if ($global_mode eq 'statistics') {
	    rmtree($tempdir);
	    rmtree($result);
	    my $exception = Job_Information_Exception->new(message => "$alignment_algorithm result available",job_type => $method, CPU_time_estimate => $CPU_estimate);
	    die $exception;
	}
	my $wrapper = Window_Pair_Alignment->new();
	my @windowpair_result = $wrapper->get_result($result, $window_pair_algorithm_parameters); # needs to return all 3
	rmtree($tempdir);
	rmtree($result);

	## alignment result sanity check. Ignore result, delete from cache + die with an error if the profiles of the result are not of the expected length
	
	my $result_is_ok;
	
	$result_is_ok = $self->private_alignment_result_sanity_check($first_gi, $second_gi, \@windowpair_result, $parameters->window_pair_algorithm_parameters);
	
	if ($result_is_ok) {
	    return @windowpair_result; # array of all 3 results (windowpairs, profile1, profile2).
	}
	else {
	    # delete result from cache # comment from laura - this hasn't occurred yet, so i don't know if this part of the code works as expected!
	    $self->private_delete_results_in_cache($keyfileMD5, $method, $username);
	    $GU->user_info(1,"Alignment result failed sanity check and has been discarded (profile length in result is not the expected length), please recompute");
	    die 'FAILED_PROFILE_LENGTH_CHECK';
	}
	## end sanity check
    } else { # if for some reason we have failed to retrieve from cache twice
             # we end up here, where we can only fail.
	    $GU->user_info(1, "Result could not be computed, most likely, we are in session mode and the result was not entered into the cache database. ");
    	die "FAILED_RESULT_NOT_COMPUTED";
    }      
  } # handle_alignment_job #
  
  method handle_APPLES_function (Str $method, Any $object, ArrayRef $parameters_ref, Boolean $cache_result, Int $cache_duration, Job_Parameters $job_parameters) {
    
    my $path_to_apples_dat = File::Spec->rel2abs($ENV{'APPLES_DAT'});
    $path_to_apples_dat =~ tr|\\|/|;
    
    my $username = ${$self->APPLES_settings}{username};

    ### NEEDS TO GET ALL NECESSARY ENVIRONMENTAL VARIABLES FROM USER + HARD-WIRE INTO THE SCRIPT 

    # to handle APPLES functions as perl scripts
    # inputs are method name, the calling object, and parameters (as an ArrayRef)
    # freeze calling object to disk, and thaw inside the perl script
    # freeze the parameters to disk, and thaw inside the perl script
    # use eval in perl script to catch exceptions
    # two output files: one for exception (freeze $@), two for results (freeze $result, freeze $object)
    # cache, scheduler, cache (proceed as before)
    
    my $inputdir = $ENV{'input_dir'};
    my $tempdir = $GU->get_temp_dir_and_use_perseverance($inputdir);
    $tempdir =~ tr|\\|/|;
    my $objectfile = "object.dump";
    nstore (\$object, "$tempdir/$objectfile") or die "Cannot store object to $tempdir/$objectfile\n"; # store calling object to file
    my $parametersfile = "parameters.dump";
    my @parameters = @{$parameters_ref};
    nstore (\@parameters, "$tempdir/$parametersfile") or die "Cannot store parameters to $tempdir/$parametersfile"; # store dereferenced parameters array to file
    my $keyfilename = 'keyfile.txt';
    my @inputs;
    push (@inputs, $object, @{$parameters_ref}, $method);
    $Storable::canonical = 1; # -> "Storable will store hashes with the elements sorted by their key"
    nstore (\@inputs, "$tempdir/$keyfilename") or die "Cannot store inputs to $tempdir/$keyfilename\n"; # serialise object, parameters and method ($inputs) to keyfile
    my $filepath = "$tempdir/$keyfilename";
    my $keyfileMD5 = $GU->md5sum($filepath);
    $GU->user_info(3,"MD5SUM: $keyfileMD5 4 $filepath\n");
    my $globalperseverancefile = "perseverance.dump";
    nstore (\$main::global_perseverance, "$tempdir/$globalperseverancefile") or die "Cannot store global perseverance object to $tempdir/$globalperseverancefile"; # store global perseverance settings to file
    

    # create a perl script using generic/exhaustive 'use ...;' statements 
    my $filename = "perlscript.pl"; 
    my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
    my %APPLES_config = $APPLES_conf->getall();
    my $apples_main = $APPLES_config{APPLES_main};
    
    my $function_actually_takes_parameters = TRUE;
    if ($#{$parameters_ref}<0) {
	$function_actually_takes_parameters = FALSE;
    }

    open (FILE, ">$tempdir/$filename") or die "Cannot create $tempdir/$filename.\n";
    
    print FILE "BEGIN {\n
     our (\$APPLES_running_mode, \$global_print_level, \$global_perseverance);
     \$global_print_level = \"3\";
      use lib \"$apples_main\";  
      use lib '/common/perl-5.10.0/bin/';
      use General_Utilities;
      my \$GU = General_Utilities->new();
      \$GU->load_includes(\"$path_to_apples_dat\");
    }\n\n";
 	
    print FILE "use APPLES_library;\n"; # APPLES_library will contain use statements for all APPLES packages
    print FILE "use APPLES_dependencies;\n"; # APPLES_dependencies will contain use statements for all Perl module dependencies in APPLES
	print FILE "use Storable qw(nstore retrieve);\n";
	print FILE "use Data::Dumper;\n";

    print FILE "our (\$APPLES_running_mode, \$global_perseverance, \$memory_cache);\n";
    print FILE "\$memory_cache = Memory_Cache->new();\n";
    print FILE "\$APPLES_running_mode = \"preparation\";";
    print FILE "\nno warnings qw(uninitialized);\n";
    print FILE "my \$tempdir = \'$tempdir\';\n";

    # thaw 'parent' global perseverance object file into the script
    print FILE "my \$parent_global_perseverance_settings = \${ retrieve (\"$tempdir/$globalperseverancefile\") } or die \"Cannot retrieve object from disk\";\n"; 
    print FILE "\$global_perseverance = \$parent_global_perseverance_settings;\n";

    # thaw object and parameters files into the script
    #print FILE "no lib \'\/System\/Library\/Perl\/Extras\/5.8.6\';\n";
    #print FILE "INC derived:\n\@INC;\n";
    print FILE "my \$object = \${ retrieve (\"$tempdir/object.dump\") } or die \"Cannot retrieve object from disk\";\n"; 
    print FILE "print \"OBJ\".Dumper(\$object);\n";
    if ($function_actually_takes_parameters) {
	print FILE "my \@parameters = \@{ retrieve (\"$tempdir/parameters.dump\") } or die \"Cannot retrieve parameters from disk\";\n"; 
	print FILE "my \@result;\n\neval{\n\t\@result = \$object->$method(\@parameters);\n};\n"; # eval statement to catch exceptions
    } else {
	print FILE "my \@result;\n\neval{\n\t\@result = \$object->$method();\n};\n"; # eval statement to catch exceptions
    }
    print FILE "print \"RES\".Dumper(\@result);\n";
    # write any $@ errors to error file 
    print FILE "if (\$\@)\{\n";  				
    print FILE "\tmy \$errorlog = \"$tempdir/errorlog\";\n";
    print FILE "\tnstore(\\\$\@, \$errorlog);\n";
    # and make separate error file (errordumperlog) of Dumper $@
    print FILE "\tmy \$errordumperlog = \"$tempdir/errordumperlog\";\n";
    print FILE "\topen (ERRORLOG, \">\$errordumperlog\") or die \"Cannot open \$errordumperlog\";\n";
    print FILE "\tprint ERRORLOG Dumper (\$\@);\n";
    print FILE "\tclose (ERRORLOG);\n";
    print FILE "\}\n"; # end if ($@)
    # write result to result file (this will be a directory to get result from)
    print FILE "else \{\n";   
    print FILE "\tmy \$resultfile = \"resultfile\";\n";
    print FILE "\tnstore (\\\@result, \$resultfile) or die \"Cannot store result to file\";\n";	
    # serialise calling object (possibly now modified) to a file on disk
    print FILE  "\tmy \$objectfile = \"objectfile\";\n";
    print FILE "\tnstore (\\\$object, \$objectfile) or die \"Cannot store object to file\";\n"; 
    print FILE "\}\n";

    close (FILE);

    chmod (0777, "$tempdir/$filename"); 
    my $path_to_perl = $ENV{'path_to_perl'}; # wsbc: "/common/perl-5.10.0/bin/perl"
    my $command_line = $path_to_perl." ".$tempdir."\/".$filename."\n"; # NB must be perl 5.10 for MooseX etc
	  #$GU->user_info(3,"->".$command_line."<-\n");
    my $result_dir = $tempdir; # i.e. /cluster/data/webservices/apples/input/some_temp_dir
    
    ### delete previous result from cache if required ###
    my $result;
    my $perl_result;
    my $global_mode = $main::APPLES_running_mode;
    my $job_mode;
    if ($job_parameters->recache) {
        $self->private_delete_results_in_cache($keyfileMD5, $method, $username);
    }
    ### check memory cache for result ###
    if ($job_parameters->cache_in_memory) {
	if (!$job_parameters->recache) {
	    $result = $main::memory_cache->retrieve($keyfileMD5, $method);
	    if (defined $result) {
		$GU->user_info(3,"Retrieved a result from memory cache (method ".$method.").\n");
		rmtree($tempdir);
		return @{$result};
	    }
	}
    }
    ### get result from cache database if needed ###
    $result = $self->private_check_cache_for_result($keyfileMD5, $method, $username);
    if ($result ne 'NULL') { # result exists, $result is a directory
	# thaw exception, result and object files. Deal with exception file first: throw exception if errors found. Else, return array of object and result.
	### NOTE - there should be no need to error handle at this point - all error handling is done BEFORE caching a result
	if ($global_mode eq 'statistics') {
	    rmtree($result);
	    rmtree($tempdir);
	    my $exception = Job_Information_Exception->new(message => "APPLES function result is available",
							   job_type => $method, CPU_time_estimate => $job_parameters->wall_time_estimate);
	    die $exception; # in statistics mode, throw an error saying the result exists
	}
	# in all other cases, return the result
	$GU->user_info (3, "result was found in cache database (method ".$method.", MD5-key: ".$keyfileMD5.").\n");	
	my $processed_object = ${retrieve ("$result/objectfile")};
	my @result = @{retrieve ("$result/resultfile")};
	push (my @final_result, $processed_object, @result);	
	File::Path::rmtree ($result); # Remove the temporary directory now the result has been extracted
	File::Path::rmtree ($tempdir);
	### insert result into memory cache ###
	if ($job_parameters->cache_in_memory) {
	    $main::memory_cache->store($keyfileMD5, $method, \@final_result);
	}
	return @final_result;
    }
    ### result was neither in memory cache nor in cache database, have to deal with running a job ###
    else { # result not available, so submit script to Scheduler
	# deal with different modes: statistics, statistics_and_retrieval, preparation, normal
	if ( ($global_mode eq 'statistics') || ($global_mode eq 'statistics_and_retrieval') ) {
	    $GU->user_info(3, "This (perl) job has not been previously requested: the result is unavailable\n");
	    rmtree($tempdir);
	    my $exception = Job_Information_Exception->new(message => "This APPLES function has not been previously requested with the required parameters, result unavailable",
							   job_type => $method, CPU_time_estimate => $job_parameters->wall_time_estimate);
	    die $exception;
	}
	else { # if preparation or normal, submit job in 'queue' mode
	    if ($global_mode eq 'preparation') {
		$job_mode = 'queue';
	    }
	    elsif ($global_mode eq 'normal') {
		$job_mode = 'session';
	    }
	    my $cache = Cache->new(cache_settings => $self->cache_settings);
	    my $scheduler = Scheduler->new(cache=>$cache);
	    $scheduler->get_config_settings();
	    # submit job (perl script) to Scheduler
	    $perl_result = $scheduler->queue_job($command_line, $method, $keyfileMD5, $result_dir, $username, $job_mode, $cache_result, $cache_duration, $job_parameters); # session=zipped result filepath; queue=job directory name
	    #$GU->user_info(3,"perl_result: ".$perl_result."\n");
	    $GU->user_info(2, "submitted perl script to scheduler\n");
	    #$GU->user_info(3,"our global mode is: ".$global_mode."\n");
	    
	    if ($job_mode eq 'queue') { 
	      my $job_info_record = {
			   MD5sum  => $keyfileMD5,
			   tempdir => $perl_result,
			   method => $method
			  };
		$GU->user_info (3, "result not currently available, appropriate job has been submitted\n");
		# result not available, so will submit job to Scheduler, and exit via a job_information_exception, providing the job_number
		my $exception = Job_Information_Exception->new(message=>"result of APPLES function not currently available, the job has been submitted",
							       job_type => $method, CPU_time_estimate => $job_parameters->wall_time_estimate, job_info => $job_info_record );
		die $exception;
	    }
	    elsif ($job_mode eq 'session') { # the result file(path) is returned - get result from here
		$GU->user_info(3, "now we're fetching the result in session mode\n");
		my @return_result = $self->private_get_APPLES_result($perl_result, $result_dir);
		rmtree($tempdir);
		return @return_result; # first element is the object, the rest is the result array
	    }
	    else {
		die 'Unknown mode.';
	    }
	}
    }
  } # handle_APPLES_function #
    
  method handle_MEME_job (Genomic_Interval_Set $genomic_interval_set, MEME_Parameters $meme_parameters, Boolean $cache_result, Int $cache_duration) {
    # need to implement via cache/scheduler route
    die 'method not implemented via Job Handler, Scehduler and Cache yet';

    # generate Job_Parameters object

    $self->get_config_settings();
    my $global_mode = $main::APPLES_running_mode;
    # check cache for result
    my @result = Cache->new(cache_settings => ${$self->cache_settings})->check_cache();
    if ($result[0] eq 'NULL') { # job needs running
       if ($global_mode eq 'statistics') {
	my $exception = Job_Information_Exception->new(message => 'meme job'); # would need to get MD5sum and queue directory of job to fill job_info attribute
	die $exception;
      }
       elsif ($global_mode eq 'preparation') {
	 die 'not configured for this option yet';
       }
       elsif ($global_mode eq 'normal') {
	 # bypassing scheduler for now...
	 my $binary = ${$self->job_handler_settings}{meme_binary};
	 my $result = MEME_wrapper->new()->run($genomic_interval_set, $meme_parameters, $binary);
       }
    }
    else {
      return @result;
    }
  } # handle_MEME_job #

  method private_check_seaweeds_job_size (Genomic_Interval $first_gi, Genomic_Interval $second_gi) {
  	# size limit removed, new seaweed binary can run jobs of any size
  	# (peterkrusche)
    return; # do nothing, the job is small enough to 
  } # private_check_seaweeds_job_size #

  method private_cpu_estimate_alignment (Str $seq1, Str $seq2, Str $algorithm) {
    if ($algorithm eq 'ott') {
      my $total = (length ($seq1)/1000) * (length ($seq2)/1000);
      my $CPU = 7200*($total/10000); #seconds
      return $CPU;
    }
    elsif ($algorithm eq 'seaweed') {
      my $total = (length ($seq1)/1000) * (length ($seq2)/1000);
      my $CPU = 720*($total/10000); #seconds
      return $CPU;
    }
    else {
      die 'unknown algorithm type\n';
    }
  } # private_cpu_estimate_alignment #

  method private_cpu_estimate_meme () {
    die 'method not implemented yet\n';
    # CPU time is quadratic with input sequence length
  } # private_cpu_estimate_meme #
 
  method private_check_cache_for_result (Str $keyfileMD5, Str $method, Str $username) {
    $GU->user_info(3,"checking cache for result (".$method.")\n");
    my $cache = Cache->new(cache_settings => $self->cache_settings);
    my %cache_settings = %{$self->cache_settings};
    my @caches;
    for my $key ( keys %cache_settings ) {
      push (@caches, $key);
    }
    my %db_handles = $cache->connect_all_caches(\@caches); # returns db handle
    
    my $result = $cache->retrieve_result(\%db_handles, $keyfileMD5, $method, $username); 
    return $result; # this is --- the temp dir?
  } # private_check_cache_for_result #

  method private_delete_results_in_cache (Str $keyfileMD5, Str $method, Str $username) {
    my $cache = Cache->new(cache_settings => $self->cache_settings);
    my %cache_settings = %{$self->cache_settings};
    my @caches;
    for my $key ( keys %cache_settings ) {
      push (@caches, $key);
    }
    my %db_handles = $cache->connect_all_caches(\@caches); # returns db handle
    
    $cache->delete_results(\%db_handles, $keyfileMD5, $method, $username); 
  } # private_delete_results_in_cache #

  method private_get_APPLES_result (Str $zip_file, Str $resultpath) { 
    # unzip + parse file
    my $cwd = getcwd(); 
    chdir $resultpath;
    $GU->user_info(3,"Trying to unzip result file:".$zip_file." at this path: ".$resultpath."\n");
    if ($^O eq 'MSWin32') {
	my $intermediateFile = substr($zip_file,0,-3);
	system ("gunzip $zip_file") == 0 || die "System error!";
	system ("tar -xf $intermediateFile") == 0 || die "System error!";
	unlink ($intermediateFile);
    }
    else {
      	system ("tar -zxf $zip_file") == 0 || die "System error!";
        unlink "$zip_file"; # delete zipped file
    }
    my $processed_object = ${retrieve ("$resultpath/objectfile")};
    my @result = @{retrieve ("$resultpath/resultfile")};
    push (my @APPLES_result, $processed_object, @result);
    chdir $cwd;
    return @APPLES_result;  
  } # private_get_APPLES_result #

  method private_alignment_result_sanity_check(Genomic_Interval $first_gi, Genomic_Interval $second_gi, ArrayRef $windowpair_result, Window_Pair_Algorithm_Parameters $parameters) {
    my $check;
    my $stepwidth1;
    my $stepwidth2;
    if ($parameters->isa('Ott_Algorithm_Parameters')) {
      $stepwidth1 = $parameters->stepwidth1;
      $stepwidth2 = $parameters->stepwidth2;
    }
    if ($parameters->isa('Seaweed_Algorithm_Parameters')) {
      $stepwidth1 = $parameters->stepwidth;
      $stepwidth2 = $parameters->stepwidth;
    }
    my $windowlength = $parameters->windowlength;
    my $profile1_check = $self->private_check_profile_length ($first_gi, $stepwidth1, $windowlength, 1, $windowpair_result);
    my $profile2_check = $self->private_check_profile_length ($second_gi, $stepwidth2, $windowlength, 2, $windowpair_result);
    if ($profile1_check && $profile2_check) {
      $check = TRUE; # everything is ok
    }
    return $check;
  } # private_alignment_result_sanity_check #

  method private_check_profile_length (Genomic_Interval $gi, Int $stepwidth, Int $windowlength, Int $profile_int, ArrayRef $windowpair_result) {
    my $sequence = $gi->get_working_sequence();
    my $seqlength = length($sequence);
    my $shouldbeprofilelength = int( ($seqlength) / $stepwidth );
    $shouldbeprofilelength -= int( $windowlength / $stepwidth ) - 1;

    while ( ( $shouldbeprofilelength > 0 ) && ( ( ( $shouldbeprofilelength - 1 ) * $stepwidth ) + $windowlength > $seqlength ) )  {
      $shouldbeprofilelength--;
    }
    if ( $shouldbeprofilelength < 0 ) {
      $shouldbeprofilelength = 0;
    }

    my @deref_result = @{$windowpair_result}; # dereference windowpair result
    my @profile_result = @{$deref_result[$profile_int]}; # select + deref actual profile
    my $isprofilelength = @profile_result; # get length of actual profile
    # check length of profile 
    my $check;
    $GU->user_info(3, "Length of profile is " . $isprofilelength . ", and should be " . $shouldbeprofilelength ."\n");
    if ( $isprofilelength == $shouldbeprofilelength ) {
      $check = TRUE; # this profile is ok
    }
    return $check;
  } # private_check_profile_length #

} # Job_Handler #
  
