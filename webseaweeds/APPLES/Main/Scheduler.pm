### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Scheduler Class ###
#
# Purpose:  Implements a scheduler which can make use of a number of slots in a qsub cluster to run jobs, or run them locally if not run on a grid.  Results from jobs are optionally cached.
#
# Jay Moore
# Version 0.9  27/10/2009
#
# Known bugs:
#
# None.
#
#
# Changes:
#
# None.
#


# If a job's job_mode is 'session' it is executed immediately, either through qsub or directly depending on $grid, and results are returned.
#ÊIf a job's job_mode is 'queued', it is passed to a series of listeners for execution, either through qsub or directly depending on $grid, and a job no is returned.
# A scheduler object relates to a directory $queue containing queued jobs
# A number of listeners are assigned to the queue, one of each listener_type.
# Each listener_type controls one step in the process of executing a job, cacheing its results, and tidying up
# Each listener uses a listener-specific semaphore in the queue to take control of the relevant step (could parallelise this by adding a listener id to listener semaphores)
# A listener scans the queue for relevant jobs, takes control of a job using a semaphore, carries out its task then releases the job on to the next listener_type
# A single listener can control several steps if it is called with $listener_types = 'request|running|complete' for example 
# In addition to listener-specific semaphores, the queue directory contains one subdirectory per job, and associated job-specific semaphores
# Job subdirectories are created with a temporary name and an associated .request semaphore by the queue_job method
# A 'request' listener watches for .request semaphores, calls qsub if on a grid else executes the command, and creates a .running semaphore
# A 'running' listener watches for completion of a running job if $grid, and creates a .complete or .error semaphore
# A 'complete' listener watches for complete jobs and caches the results if needed and creates a .trash semaphore
# An 'error' listener watches for failed jobs and logs the error and creates a .trash semaphore
# A 'trash' listener removes cached and failed jobs from the queue, and any associated semaphores
# Listeners respond to a file stopfile.listener_types in the queue by quitting

# each job in the queue has: calling object, called method, set of parameters, [submitter], [start datetime], [finish datetime], status
# status can take one of these values:
# request - job is awaiting start by listener
# running - job has been started by listener
# complete - job has finished and results have not been cached
# error - job has finished and failed
# trash - job is ready for deletion

# parameters of queued jobs are written to a parameter file in the job subdirectory, which is added to with further status updates as the job progresses

use MooseX::Declare;

my $debug = 1;

class Scheduler {
  use constant {FALSE => 0,
		TRUE  => 1};	
  use APPLES_Datatypes qw (Status JobMode ResultType Boolean ListenerType MethodType CommandLine);
  use General_Utilities;
  use Cache;
  use Parameters;
  use Cwd;
  use File::Temp qw/ tempfile tempdir /;
  use Config::General;
  use File::Basename;
  use File::Path qw (mkpath rmtree);
  use Data::Dumper;
  use File::Copy;

  has 'settings' => (is => 'rw', isa => 'HashRef'); # software settings (paths etc), must be filled by calling get_config_settings after instantiation
  has 'cache' => (is => 'ro', isa => 'Cache', required => 1);

  my $GU = General_Utilities->new;
  #Cache->new(); # needs cache_settings, obtained from Job_Handler

  method get_config_settings() {
    #$GU->user_info(3,"getting config settings!\n");
    
    my $APPLES_DAT = $ENV{'APPLES_DAT'};
    my $APPLES_conf = new Config::General($APPLES_DAT);

    #my $APPLES_conf = new Config::General("APPLES.dat") or die 'cannot open APPLES.dat config file';
    # read APPLES config file
    my %APPLES_config = $APPLES_conf->getall();
    # determine application from APPLES config file
    #my $username = $APPLES_config{username};
    #my $application = $APPLES_config{application}; # Name of application (e.g. 'APPLES') read from config (APPLES.dat)
    #my $queue = $APPLES_config{queue}; # Directory path to queue read from config
    #my $queue_id = $APPLES_config{queue_id}; # Identifier for queue read from config
    #my $cluster_dir = $APPLES_config{cluster_dir}; # Directory path to qsub command read from config
    #my $grid = $APPLES_config{grid}; # Are we running on a grid or not? read from config
    #my $job_parameters_file = $APPLES_config{job_parameters_file}; # Name of file to save job parameters in ('job_parameters.txt')
    
    #my $job_check_delay = $APPLES_config{job_check_delay}; #seconds - brought in from config
    #my $grid_type = $APPLES_config{grid_type};# config parameter for grid type: PBS (CSC francesca), SGE (wsbc)

    #my $cache_name = $APPLES_config{scheduler_cache_name}; # cache database for Scheduler, i.e. to store results

    $self->settings(\%APPLES_config);

    # this file contains:
    # calling object MD5 (on initiation)
    # calling_method (on initiation)
    # parameters MD5 (on initiation)
    # results filename (on initiation)
    # user ID (on initiation)
    # queued timestamp (on initiation)
    # cache results (y/n) (on initiation)
    # cache duration (days?) (on initiation)
    # cluster job no (on running)
    # started timestamp (on running)
    # ended timestamp (on completion)	
  } # get_config_settings #	

  method queue_job(Str $command_line, Str $method, Str $keyfile, Str $result_file, Str $user_id,
	JobMode $job_mode, Boolean $cache_results, Int $cache_duration, Job_Parameters $job_parameters) {

    # (LB) - there will now be only ONE keyfile, which will be parameters_object or just $parameters, depending on if it's a binary or a Perl script that is being run.
    # Add a job to the queue, or execute it immediately
    # $command_line, $method, $parameters_file, $calling_object_file, $result_file are parameters needed to run the job
    # $job_mode may be 'session' or 'queue':
    # Execute a 'session' job immediately, wait for results, and return them; pass a 'queue' job to a listener for execution and return $job
    # On completion, cache results if $cache_results, and add $cache_duration to the cached object

	# Tidy up and validate incoming parameters
	### $command_line must exist (and not do anything dangerous?)
	### $parameters_file must exist
	### $calling_object_file must exist
	### $method?
	### $cache_results t/f
	### if $cache_results then $cache_duration must be > 0, or use a default?

	### may need to add an outputdirpath to $command_line

      # first make sure data for job parameters file is well-defined
          chomp $method;
	  if (($keyfile eq '')||($method eq '')||($result_file eq '')||($user_id eq '')) {
	      die 'One parameter was an empty string - not allowed.';
	  }

    my $cache = $self->cache; #debugging
    my $grid_type = ${ $self->settings }{grid_type};
    $GU->user_info( 3, "we're running on grid type: " . $grid_type . "\n" );

    my $queue = ${ $self->settings }{queue};
    
    # create a directory for the job to run in
    my $job = tempdir( DIR => $queue ) or die "Cannot create job folder in $queue.\n";
    chmod( 0777, "$job" );
    if ($grid_type eq 'PBS') {
      chmod ( 0777, "$job" );
    }
    ### JOB PRIORITISATION NOTE:
    ### Jobs are placed in the job queue.
    ### The job queue is sorted, ascibetically: ! 0 9 A Z _ a z ~
    ### The jobs are then run in that order (i.e. not by submission timestamp!), which is essentially random.
    ### However, therefore you can use these letters (! or 0, z or ~) to prefix tempdirs, to alter the prioritisation of your jobs accordingly.
    ### (Listeners need to be stopped and started to reset the loop to the beginning).

    ### example temporary fix - prepend zzz to tempdir if user_id eq laurabaxter (should then leave these jobs until last to be picked up)

    if ( $user_id eq 'laurabaxter' ) {
      $GU->user_info( 3, "prepending zzz to job name" );
      my $dir_name = basename($job);
      
      #$GU->user_info(3,$dir_name."\n");
      my $renamed_dir = "zzz" . $dir_name;
      
      #$GU->user_info(3,$renamed_dir."\n");
      # make the directory named above
      my $renamed_dir_fullpath = $queue . '/' . $renamed_dir;
      mkpath($renamed_dir_fullpath);
      
      # remove original temp directory
      rmtree($job);
      
      # assign new directory path to $job
      $job = $renamed_dir_fullpath;
      chmod( 0777, "$job" );    # allows group to read/write
      if ($grid_type eq 'PBS') {
	chmod( 0777, "$job" );
      }

    }
    ### end temporary fix
    
    $GU->user_info( 3, "job->" . Dumper($job) );
    
    # we remove path from $job
    my ( $jobbase, $jobpath, $jobtype ) = fileparse($job);
    $job = "$jobbase";
    $GU->user_info( 3, "job now->" . Dumper($job) );
    my $job_parameters_file = ${ $self->settings }{job_parameters_file};
    
    # write the parameters, user_id and timestamp to a parameters file in the job directory
    open( JOBPARAMS, ">$queue/$job/$job_parameters_file" ) or die "Cannot create $queue / $job / $job_parameters_file.\n";
    print JOBPARAMS "$keyfile\n$method\n$result_file\n$user_id\n" . $GU->timestamp . "\n$cache_results\n$cache_duration";
    close(JOBPARAMS);
    chmod( 0777, "$queue/$job/$job_parameters_file" );
    if ($grid_type eq 'PBS') {
	chmod( 0777, "$queue/$job/$job_parameters_file" );
      }
    my $executable =  ${ $self->settings }{application}.'.queue.'.${ $self->settings }{queue_id}; # executable gets assigned this name so it can be identified in qstat
    # write the command to a shell script in the job directory
    if ( $^O eq 'MSWin32' ) { 
      $executable = $executable . '.bat';
    }
    $GU->user_info( 3, "Writing shell script to the job directory!\n" );
    open( JOBCOMMAND, ">$queue/$job/$executable" ) or die "Cannot create executable $queue/$job/$executable.\n";
    $GU->user_info( 3, "Command line: " . $command_line );
    # make header line necessary for IBM cluster
    my $header_line = $self->private_add_job_information_line_to_shell_script($job_parameters);
    # add lines necessary for IBM cluster
    if ( $grid_type eq 'PBS' ) {
      my $change_dir_line = "cd " . "$queue/$job" . "\n";
      print JOBCOMMAND $header_line;
      print JOBCOMMAND $change_dir_line; # or cd $PBS_O_WORKDIR ?
    }
    if ($grid_type eq 'SGE') {
      my $bash_line = "#!/bin/bash";
      print JOBCOMMAND $bash_line."\n";
      my $set_shell_line = "#\$ \-S /bin/bash";
      print JOBCOMMAND $set_shell_line."\n";
      my $export_lib = "export DYLD_LIBRARY_PATH=/cluster/laurabaxter/binaries/Seaweeds/bin/alignmentplots_769"; # this path is currently hard-wired, but it should probably be a config
      print JOBCOMMAND $export_lib."\n";
      print JOBCOMMAND 'echo $DYLD_LIBRARY_PATH'."\n"; # for debugging
    }

    print JOBCOMMAND 'pwd' . "\n"; # for debugging
    print JOBCOMMAND 'echo $PERL5LIB' . "\n"; # for debugging
    print JOBCOMMAND $command_line . "\n";
    close(JOBCOMMAND);
    chmod( 0775, "$queue/$job/$executable" ); # alter permissions
    if ($grid_type eq 'PBS') {
      chmod ( 0777, "$queue/$job/$executable" );
    }
    my $grid = ${ $self->settings }{grid};
    $GU->user_info( 3, "!!! MARKER 1 !!!\n" );
    if ( $job_mode eq 'queue' ) {    # This means queue the job for execution by a listener
      # create a request semaphore for the request listener
      open( REQUESTLOCK, ">$queue/$job.request" ) or die "Cannot create $queue/$job.request.\n";
      $GU->user_info( 3, $queue . "/" . $job . "request" . "\n" );
      close(REQUESTLOCK); # request is now in the queue ready for the listener
      chmod( 0660, "$queue/$job.request" );
      if ($grid_type eq 'PBS') {
	chmod ( 0666, "$queue/$job.request" );
      }
      $GU->user_info( 3, "$job.request\n" );
      $GU->user_info( 3, "Job:\n" );
      $GU->user_info( 3, Dumper($job) );
      return $job;
      
      # returns $job if it managed to add the command file to the queue
    }
    elsif ( $job_mode eq 'session' ) {    # This means execute the job immediately, wait, and return the result
      # create a session semaphore to lock the job and allow group permissions
      open( SESSIONLOCK, ">$queue/$job.session" )  or die "Cannot create $queue/$job.session.\n";
      chmod( 0660, "$queue/$job.session" );
      if ($grid_type eq 'PBS') {
	chmod ( 0666, "$queue/$job.session" );
      }
      my ( $queued_timestamp, $started_timestamp, $ended_timestamp, $results );
      if ($grid) {    # we are on a grid so pass the job to listener for execution
	# allow read write permissions for group on target result directory
	chmod( 0770, $result_file );
	if ($grid_type eq 'PBS') {
	  chmod ( 0777, $result_file );
	}
	# create a request semaphore for the request listener, and allow group permissions
	open( REQUESTLOCK, ">$queue/$job.request" ) or die "Cannot create $queue/$job.request.\n";
	chmod( 0660, "$queue/$job.request" );
	if ($grid_type eq 'PBS') {
	  chmod ( 0666, "$queue/$job.request" );
	}
	close(REQUESTLOCK); # request is now in the queue ready for the listener

	### put in a timeout?
	while (( !( -e "$queue/$job.error" ) )  && ( !( -e "$queue/$job.trash" ) ) ) { # loop until job is complete or in error
	  sleep ${ $self->settings }{job_check_delay}; # sleep (seconds) to save CPU while waiting for listeners to execute job
	}
	
	my $job_status;
	
	if ( -e "$queue/$job.complete" ) {
	  $job_status = 'complete';
	  die "Job $queue $job is in complete state for some reason...\n";
	}
	elsif ( -e "$queue/$job.error" ) {
	  $job_status = 'error';
	  die "Job $queue $job is in error state for some reason...\n";
	}
	elsif ( -e "$queue/$job.trash" ) {
	  $job_status = 'complete';
	}
	else {
	  die "Cannot determine job status of $queue $job\n";
	}
	my %jobs;
	my %other_cluster_jobs;
	my $running_jobs = $self->private_get_cluster_job_status( \%jobs,
								  \%other_cluster_jobs, $job_status, TRUE ); # (Int $running_jobs not needed)
	$queued_timestamp  = $jobs{$job}{queued_timestamp};
	$started_timestamp = $jobs{$job}{started_timestamp};
	$ended_timestamp   = $jobs{$job}{ended_timestamp};
      }
      
      else {     # not on a grid
	$queued_timestamp = $GU->timestamp;
	my $cwd = getcwd();
	chdir "$queue/$job" or die "Cannot chdir $queue/$job.\n";
	$started_timestamp = $GU->timestamp;
	chmod( 0755, $executable );
	if ($grid_type eq 'PBS') {
	  chmod ( 0777, $executable );
	}
	if ( $^O eq 'MSWin32' ) {
	  system("$executable > $executable.log 2> $executable.err") == 0
	    or die "Cannot execute $queue/$job/$executable.\n :$?";
	}
	else {
	  system("./$executable > $executable.log 2> $executable.err") == 0
	    or die "Cannot execute $queue/$job/$executable.\n :$?";
	}
	$ended_timestamp = $GU->timestamp;
	chdir $cwd or die "Cannot chdir $cwd.\n";
	
	if ($cache_results) {
	  
	  my $cache = $self->cache;
	  
	  my %all_cache_settings = %{ $cache->cache_settings };
	  my $cache_name = ${ $self->settings }{scheduler_cache_name};
	  my %cache_settings = %{ $all_cache_settings{$cache_name} };
	  
	  my $tablename = $cache_settings{db_cache_table};
	  
	  # zip contents of result dir and store in cache
	  my $cwd = getcwd();
	  chdir "$queue/$job" or die "Cannot change into $queue/$job\n";
	  
	  #$GU->user_info(3,"now in $queue/$job\n");
	  #$GU->user_info(3,"job: $job\n");

	  # remove stdout and stderr output files
	  unlink("$executable.log");
	  unlink("$executable.err");

	  if ( $^O eq 'MSWin32' ) {
	    system("tar -cf ../$job.results *") == 0 || die "System error!";
	    system("gzip ../$job.results") == 0 || die "System error!";
	  }
	  else {
	    system("tar -zcf ../$job.results.gz *") == 0 || die "System error!";
	  }
	  chdir $cwd or die "Cannot change directory to $cwd";
	  
	  my $db = $cache->connect($cache_name);
	  my $result = "$queue/$job.results.gz";
	  $GU->user_info( 3, "$db, $keyfile, $method, $user_id, $queued_timestamp, $started_timestamp, $ended_timestamp, $cache_duration, $result, $tablename\n");
	  $cache->store_result($db, $cache_name, $keyfile, $method, $user_id, $queued_timestamp, $started_timestamp, $ended_timestamp, $cache_duration, $result, $tablename);
	  $cache->disconnect($db, $cache_name);

	  # ought to check cache, to see that result has been stored ok
	}
	File::Path::rmtree("$queue/$job");
      }    # end not on grid
      #chmod (0775, $result_file); # grant permission to write to this folder
      #move ("$queue/$job.results.gz", $result_file) or die "Copy failed for $queue $job.results.gz to $result_file: $!"; # copy the result to $result_file location - this is provided by user. Currently failing to write.
      close(SESSIONLOCK);    # session is now released
      unlink "$queue/$job.session";
      
      ###ÊCheck for existence of $queue/$job.error - if exists report error
      ###ÊCheck content of $queue/$job/$executable.err - if not empty report error
      
      #rmtree("$queue/$job") or die "Cannot remove $queue/$job.\n";
      
      ### remove all semaphores? are there any?
      if ($grid) {
	$GU->user_info( 3,"returning $result_file/$job.results.gz to user\n" );
	return "$result_file/$job.results.gz"; # this is the zipped file location
      }
      else {
	$GU->user_info( 3, "returning $queue/$job.results.gz to user\n" );
	return "$queue/$job.results.gz";  # this is the zipped file location
      }
    }    # end elsif job mode eq 'session'
    else {
      die "Don't understand job mode $job_mode.\n";
    }
  }    # queue_job #
  
  method cancel_job (Str $job, Str $user_id) {
    
    my $success;
    ### $job = tempfile;
    ### kill qsub job if there is one
    ### remove all semaphores, waiting for lock if needed
    ### remove job
    
    # if job is in request state:
    #  get lock
    #  delete job folder
    #  release lock
    #  $success = 1;
    # elsif ($state eq 'running')
    #  get lock
    #  cancel qsub job
    #  change state to 'error'? 'cancelled'?
    #  $success = 1;
    die 'Method not implemented yet!';
    ###Êdo something with user id?
    return $success;
    # returns $job if it managed to add the command file to the queue
  } # cancel_job #

  method retrieve_queue (Str $user_id) {

    # Get details of the jobs in a queue
	  
    my %jobs;
    my %other_cluster_jobs;
    # get list of each type of jobs in queue
    my %queue_jobs;
    my @listener_types = qw(request running complete error trash);
    foreach my $listener_type (@listener_types) {
      my $running_jobs = $self->private_get_cluster_job_status(\%jobs, \%other_cluster_jobs, $listener_type, TRUE); # (Int $running_jobs not needed)
      $queue_jobs{$listener_type} = %jobs;
    }
    return %queue_jobs;
  } # retrieve_queue #

  method listen (Str $listener_types, Int $max_allowed_jobs, Int $no_available_nodes, Int $min_free_nodes, Int $listen_loop_delay, Str $username) {
    my $cache = $self->cache;

    # Start up a listener and have it do stuff
    my %jobs;
    my %other_cluster_jobs;
    my $queue = ${$self->settings}{queue};
    my $grid = ${$self->settings}{grid};
    my $job_parameters_file = ${$self->settings}{job_parameters_file};
    my $executable = ${$self->settings}{application}.'.queue.'.${$self->settings}{queue_id}; # executable gets assigned this name so it can be identified in qstat
    my $cluster_dir = ${$self->settings}{cluster_dir};
    # sort out which listener type(s) the listener will respond to
    my @listener_types = split (/_/, $listener_types);
    my $listening = 1; # anything goes wrong on startup and this gets set to 0, listener doesn't start 
    my $grid_type = ${$self->settings}{grid_type};
    # try to get a lock on the queue with a semaphore for each relevant listener type

    # neeed to open filehandles before we get into listener_type loop (and not close them), to prevent multiple listeners of same type from being started up
    my @queuelocks;
 
    foreach my $listener_type (@listener_types) {
      if ( ($listener_type eq 'request') || ($listener_type eq 'running') || ($listener_type eq 'complete') || ($listener_type eq 'error') || ($listener_type eq 'trash') ) { # $listener_type is a valid listener_type
	open (QUEUELOCK, ">$queue/$listener_type.lock"); # if we can lock one of the relevant listener types, do start
	if (flock QUEUELOCK,2) {
	  push (@queuelocks, *QUEUELOCK);
	  chmod (0664, "$queue/$listener_type.lock");#chmod
	  if ($grid_type eq 'PBS') {
	    chmod ( 0666, "$queue/$listener_type.lock" );
	  }
	}

	else {
	  $listening = 0;
	}
      } 
      else { # $listener_type is not a valid listener_type
	$listening = 0;
	die "$listener_type is not a valid listener_type.\n";
      }
    }
	  
    if ($listening) {		
      # declare a database handle in case we need one
      my $db;

      # on restart of listener, tidy up the queue:
	    
      # remove previous stop request if there is one
      unlink "$queue/$listener_types.stop";

      # tidy up previous jobs current listener should be dealing with
      foreach my $listener_type (@listener_types) {
	if ($listener_type eq 'complete') {
	  # start up a connection to the cache
	  my $cache_name = ${ $self->settings }{scheduler_cache_name};
	  $db = $cache->connect($cache_name);
	}
	# get list of jobs in my queue, and other jobs on cluster
	$self->private_get_cluster_job_status(\%jobs, \%other_cluster_jobs, $listener_type, TRUE);
	foreach my $job (keys %jobs) {
	  if ($listener_type eq 'request') {
	    ### if job has status 'request' but is really running cancel the qsub job
	    ### if it is still running, cancel it 
	  }
	  if ($listener_type eq 'running') {
	    ### if job has status 'running' but is not really running rename its semaphore to 'error', 'complete' or 'request'?
	    ### if it is still running, cancel it 
	  }
	  if ($listener_type eq 'complete') {
	    ### if job has status 'running' but is not really running rename its semaphore to 'error', 'complete' or 'request'?
	    ### if it is still running, cancel it
                  
	   
	  }
	}
	### if listener_type eq 'complete' check whether results cached and if so move it on to trash
      }
      my $timestamp = $GU->timestamp;
      #$GU->user_info(3,$GU->timestamp." : $listener_types listener is listening to queue $queue.\n"); ### does this need to be print?
      $GU->user_info(3,"$timestamp : $listener_types listener is listening to queue $queue. Started by $username.\n");
      # ADD THIS INFO TO job_errors.log AS WELL:
      my $APPLES_management_home = ${ $self->settings }{APPLES_management_home};
      open (REPORTERRORS, ">>$APPLES_management_home/job_errors.log") or die 'Cannot open job_errors.log\n';
      print REPORTERRORS "$timestamp : $listener_types listener is listening to queue $queue. Started by $username.\n";
      close (REPORTERRORS);
      # keep listening until told to stop
      while (!(-e "$queue/$listener_types.stop")) {
		
	my $jobs_found = 0;
                
	# go through each listener type and process its queue    
	foreach my $listener_type (@listener_types) {

	  # get list of jobs of my type in my queue, and also other jobs already running on cluster if applicable
	  $self->private_get_cluster_job_status(\%jobs, \%other_cluster_jobs, $listener_type, TRUE);
                    
	  # loop through my queued jobs
	  foreach my $job (sort keys %jobs) {
	    $jobs_found++;

	    # try to get a lock on the job if it is still waiting for my $listener_type
	    if ((-e "$queue/$job.$listener_type") && (open(JOBLOCK, ">$queue/$job.$listener_type"))) {
	      #$GU->user_info(3,"$listener_type $job\n");
	      # set up variable to move job on to next listener once this listener done
	      my $next_step;

	      # Do the specific action for the relevant $listener_type
	      if ($listener_type eq 'request') {
		my $started = 0;
		while (!$started) { # keep going round loop until a slot has been found and job submitted to grid, or if not grid, job has been run
		  if ($grid) {
				
		    # count $my_running_jobs owned by $scheduler_owner, and derive $no_free_nodes from $no_available_nodes:

		    # retrieve current qstat details
		    my $my_running_jobs = $self->private_get_cluster_job_status(\%jobs, \%other_cluster_jobs, $listener_type, FALSE);# to restore original behaviour set to TRUE
		    my $no_free_nodes = $no_available_nodes;
			   
		    foreach my $my_job (keys %jobs) {
		      if ($jobs{$my_job}{active}) {
			$no_free_nodes--;
		      }	
		    }
		    foreach my $cluster_job (keys %other_cluster_jobs) {
		      $no_free_nodes--;
		    }
		    if ($no_free_nodes < 0) {
		      $no_free_nodes = 0;
		    }
		    # test whether there are free slots on the cluster, and we haven't exceeded our quota
		    if (($max_allowed_jobs > $my_running_jobs) && ($no_free_nodes >= $min_free_nodes)) {

		      # create a temporary file to capture cluster job no
		      my ($temp_filehandle, $temp_filename);
		      ($temp_filehandle, $temp_filename) = tempfile();
		      close $temp_filehandle;

		      # initiate cluster job
		      my $this_dir = getcwd();
		      chdir "$queue/$job" or die "Cannot change directory to $queue/$job.\n";
		      chmod (0755, $executable);
		      if ($grid_type eq 'PBS') {
			chmod ( 0777, $executable );
		      }
		      if ($grid_type eq 'SGE') {
			system ("qsub -cwd $executable > $temp_filename"); 
		      }
		      elsif ($grid_type eq 'PBS') {
			system ("qsub $executable > $temp_filename"); 
		      }
		      chdir $this_dir or die "Cannot change directory back to $this_dir.\n";

		      # capture cluster job number from unique filename
		      open (JOBNOFILE, $temp_filename);
		      my @jobno = <JOBNOFILE>;
		      close (JOBNOFILE); 
		      unlink ($temp_filename);
 
		      my @line = split / /, $jobno[0];
		      my $cluster_job = $line[2];
		      if ($grid_type eq 'PBS') {
			$cluster_job = $line[0];
			$cluster_job =~ s/\..*//; # remove .* tag (e.g. .fe1), as we're expecting numerical job numbers
		      }
		      # write cluster job number and timestamp to file in job directory
		      open (JOBFILE, ">>$queue/$job/$job_parameters_file");
		      print JOBFILE "\n$cluster_job\n".$GU->timestamp; # for PBS $cluster_job is currently undef - needs fixing
		      close (JOBFILE);

		      # save cluster job number and timestamp to %jobs
		      $jobs{$job}{cluster_job} = $cluster_job;
		      $jobs{$job}{started_timestamp} = $GU->timestamp;                    
		      # get ready to move on, and pass job to next listener in chain
		      $next_step = 'running';
		      $started = 1; # job has been started so we can move on
		    }		# end of there is a free slot on grid
		  }		#Êend of grid
		  else { # queue mode, not on grid
		    my $cwd = getcwd();
		    chdir "$queue/$job" or die "Cannot change directory to $queue/$job.\n";
		    chmod (0755, $executable);
		    if ($grid_type eq 'PBS') {
		      chmod ( 0777, $executable );
		    }
		    system ("./$executable > $executable.log 2> $executable.err") == 0 or die "Cannot execute $queue/$job/$executable.\n";
		    chdir $cwd or die "Cannot change directory back to $cwd.\n";
		    # get ready to move on, and pass job to next listener in chain
		    $next_step = 'running';
		    $started = 1; # job has been started (finished actually) so we can move on
		  }		  # end of not grid
		  sleep $listen_loop_delay; # sleep 
		}	     # end of loop to keep trying to start job
	      }		     # end of listener_type request
	      elsif ($listener_type eq 'running') {
		my $job_status = 'running';
		if ($grid) {
		  # retrieve current qstat details
		  $self->private_get_cluster_job_status(\%jobs, \%other_cluster_jobs, $listener_type, FALSE);# to restore original behaviour set to TRUE
		  if ($jobs{$job}{active} == FALSE) { # job has finished                           
		    # decide whether job succeeded or failed
		    $GU->user_info(3,"looks like job has finished running. Examining $queue / $job / $executable\n");
		    my $perlscript = 'false';
		    # for perl scripts, success=presence of resultfile and objectfile (regardless of "errors" (warnings) in the *.e* file)
		    # 1 - determine if job is a perl script (e.g. check if $queue/$job/$executable file contains 'perlscript.pl'?):
		    
		    open(EXECUTABLEFILE, "$queue/$job/$executable") or die "Cannot open $queue/$job/$executable\n";
		    if (grep{/perlscript/} <EXECUTABLEFILE>){
		      $perlscript = 'true';
		      $GU->user_info(3,"job was a perl script\n");
		    }
		    close EXECUTABLEFILE;
		    # 2 - check for presence of file called resultfile in $queue/$job directory
		    if ($perlscript eq 'true') {
		      if (-e "$queue/$job/resultfile") {
		        $GU->user_info(3,"$queue / $job / resultfile exists!\n");
			$job_status = 'complete';
		      }
		      else {
			$GU->user_info(1,"resultfile $queue / $job / resultfile does not exist!\n");
			$job_status = 'failed';
		      }
		    }
		    
		    else { # else (not a perlscript), check for an empty qsub error log
		      $job_status = 'failed';
		      ### begin temporary seaweed-fix:
		      # Seaweed binary is currently outputting harmless progress statements to the *.e* file. Until this is fixed, we cache the result anyway, and rely on the profile length sanity check in the job handler being sufficient to deal with erroneous results (LB 30/03/2010)
		      my $seaweed = 'false';
		      #determine if job was a seaweed job:
		      open(EXECUTABLEFILE, "$queue/$job/$executable") or die "Cannot open $queue/$job/$executable\n";
		      if (grep{/Seaweeds/} <EXECUTABLEFILE>){
			$seaweed = 'true';
			$GU->user_info(3,"job was a seaweed alignment\n");
		      }
		      close EXECUTABLEFILE;
		      # see if result.txt exists, and assume computation was ok if so
		      if ($seaweed eq 'true') {
			if (-e "$queue/$job/result.txt") {
			  $GU->user_info(3,"$queue / $job / result.txt file exists!\n");
			  $job_status = 'complete';
			}
			else {
			  $GU->user_info(1,"resultfile $queue / $job / result.txt does not exist!\n");
			  $job_status = 'failed';
			}
		      } # end if seaweed
		      ### end of temporary seaweed-fix
		      else {
			# for all other non-perl script, non-seaweed binary results, check for an empty qsub error log
			if (-e "$queue/$job/$executable.e".$jobs{$job}{cluster_job}) { # if the result file exists, conclude that the job succeeded
			  # open, read into array, count rows, if zero row, no errors, all is well. else 'failed';
			    my $job_error_file = "$queue/$job/$executable.e".$jobs{$job}{cluster_job};
			  open ( ELOG, $job_error_file ) or die "Cannot open error ".$job_error_file."\n";
			  
			  my @error_lines=<ELOG>;
			  close (ELOG);
			  if (scalar(@error_lines == 0)) {
			    $job_status = 'complete'; # error file is empty: status->complete
			  }
			  else { # error file non-empty: assume job failed
			    $job_status = 'failed';
			  }
			}
			else { # file doesn't exist
			  $job_status = 'failed';
			}
		      } # end else (non-perl non-seaweed)
		    #$GU->user_info(3,"job status = $job_status\n");
		    } # end else not a perlscript
		  }
		  else {
		    $next_step = 'unchanged';
		  }
		} # end of grid
		else { # if status eq 'running' and we are not on a grid then job has already finished.  If job were really running on non-grid, status would still be 'request' and semaphore would be locked by listener			   
		  # decide whether job succeeded or failed 
		  # Check content of $queue/$job/$executable.err - if not empty $job_status = 'failed' else $job_status = 'complete';
		  open (ERRORLOG, "$queue/$job/$executable.err") or die "Cannot open error $queue/$job/$executable.err.\n";
		  my @errors = <ERRORLOG>;
		  if ($#errors > 0) { # there is an error
		    $job_status = 'failed';
		  } 
		  # if it's a perl script that's run, also need to check apples/input/tmp/errordumperlog? or store results with an error flag?
		  
		  else {
		    $job_status = 'complete';
		  }
		  open (JOBFILE, ">>$queue/$job/$job_parameters_file") or die "Cannot open $queue/$job/$job_parameters_file.\n";
		  print JOBFILE "\n\n".$GU->timestamp;
		  close (JOBFILE);
		} # end of not grid
		if ($job_status eq 'failed') {		  
		    $next_step = 'error';
		} elsif ($job_status eq 'complete') {
		  # write timestamp to file in job directory
		  open (JOBFILE, ">>$queue/$job/$job_parameters_file") or die "Cannot open $queue/$job/$job_parameters_file.\n";
		  print JOBFILE "\n".$GU->timestamp;
		  close (JOBFILE);
		  $next_step = 'complete';
		}
	      }			# end of listener_type running
	      elsif ($listener_type eq 'complete') {
		if ($jobs{$job}{cache_results}) {
		  my $cache = $self->cache;
		  my %all_cache_settings = %{$cache->cache_settings};
		  my $cache_name = ${ $self->settings }{scheduler_cache_name};

		  my %cache_settings = %{$all_cache_settings{$cache_name}};
	
		  my $tablename = $cache_settings{db_cache_table};
		  my $cwd = getcwd();
		  chdir "$queue/$job" or die "Cannot chdir into $queue/$job";
		  $GU->user_info(3,"now we're in $queue/$job\n");
		  $GU->user_info(3,"our job: $job\n");

		  # remove stdout and stderr output before caching
		  my $stdout_filename = "$executable.o".$jobs{$job}{cluster_job};
		  my $stderr_filename = "$executable.e".$jobs{$job}{cluster_job};
		  unlink($stdout_filename);
		  unlink($stderr_filename);

		  if ($^O eq 'MSWin32'){
		    system ("tar -cf ../$job.results *");
		    system ("gzip ../$job.results");
		  }
		  else {
		    system ("tar -zcf ../$job.results.gz *");
		  }
		  chdir $cwd or die "could not change directory to $cwd";
		  my $result = "$queue/$job.results.gz";
		  $cache->store_result($db, $cache_name, $jobs{$job}{keyfile}, $jobs{$job}{calling_method}, $jobs{$job}{user_id}, $jobs{$job}{queued_timestamp}, $jobs{$job}{started_timestamp}, $jobs{$job}{ended_timestamp}, $jobs{$job}{cache_duration}, $result, $tablename);
		  if (-e "$queue/$job.session") {
		    my $destination_file = $jobs{$job}{result_file};
		    move ($result, $destination_file) or die "Copy failed for $queue $result to $destination_file: $!"; # copy the result to $result_file location - this is provided by user.;
		    $GU->user_info(3,"moved $queue $result to $destination_file\n");
		  }
		  else {
		    unlink ($result); # deletes zip file once cached #
		  }
		  # check result is now in cache
		  # call $cache->check_cache();  # returns a string/int
		  if ($cache->check_cache($db, $cache_name, $jobs{$job}{keyfile}, $jobs{$job}{calling_method}, $jobs{$job}{user_id}, $tablename) > 0){
		    $next_step = 'trash';
		  }
		  else {
		   $next_step = 'error';
		  }
		 
		} else { # no need to cache, but still need to pass results back to where user asked if in session mode
		  if (-e "$queue/$job.session") {
		    my $cwd = getcwd();
		    chdir "$queue/$job" or die "Cannot chdir into $queue/$job";
		    $GU->user_info(3,"now we're in $queue/$job\n");
		    $GU->user_info(3,"our job: $job\n");
		    if ($^O eq 'MSWin32') {
		      system ("tar -cf ../$job.results *");
		      system ("gzip ../$job.results");
		    }
		    else {
		      system ("tar -zcf ../$job.results.gz *");
		    }
		    chdir $cwd or die "change directory to $cwd";
		    my $result = "$queue/$job.results.gz";
		    my $destination_file = $jobs{$job}{result_file};
		    move ($result, $destination_file) or die "Copy failed for $queue $result to $destination_file: $!"; # copy the result to $result_file location - this is provided by user.;
		  }
		}
		if (open (SESSIONLOCK, ">$queue/$job.session")) { # job is not locked by a session, so pass it on to next listener
		  close (SESSIONLOCK);
		  unlink "$queue/$job.session";
		  $next_step = 'trash';
		}
	      }			#Êend of listener_type complete
	      elsif ($listener_type eq 'error') {
		### log error: "Job $job in queue $queue had error $error.\n";
		if (open (SESSIONLOCK, ">$queue/$job.session")) { # job is not locked by a session, so pass it on to next listener
		  close (SESSIONLOCK);
		  unlink "$queue/$job.session";
		  #$next_step = 'trash'; # for now, leave undefined, so we can diagnose errors
		  # write $job to an error log here
		  my $APPLES_management_home = ${ $self->settings }{APPLES_management_home};
		  open (REPORTERRORS, ">>$APPLES_management_home/job_errors.log") or die 'Cannot open job_errors.log\n'; 
		  print REPORTERRORS "$job";
		  print REPORTERRORS "\n";
		  close (REPORTERRORS);
		}
	      } # end of listener_type error
	      elsif ($listener_type eq 'trash') {
		# remove the job
		rmtree ("$queue/$job") or die "Cannot remove $queue/$job.\n";
		### also remove any semaphores associated with job, except the current one ($job.trash) which we have lock on - this will be dealt with later anyway by the 'no next_step' code
	      } # end of listener_type trash

	      # end of actions specific to each listener_type, go on to do stuff to the job which is common to all listener_types
		       
	      # if there is a next step, try to pass the job to the next listener
	      if (defined $next_step) {
		if ($next_step ne 'unchanged') {
		  if (open (JOBNEXT, ">$queue/$job.$next_step")) { # try to get a lock on the job's semaphore for the next listener
		    close (JOBLOCK); # release the current listener's lock on the job
		    unlink "$queue/$job.$listener_type"; # remove the current listener's semaphore
		    close (JOBNEXT); # release the job to the next listener
                    chmod (0775, "$queue/$job.$next_step");
		    if ($grid_type eq 'PBS') {
		      chmod ( 0777, "$queue/$job.$next_step" );
		    }
		  } else {	# Can't pass job on.  Abort mission.
		    close (JOBLOCK);
		    unlink "$queue/$job.$listener_type";
		    die "Cannot move job $job in queue $queue from $listener_type to $next_step.\n";
		  }
		}
	      }			# end of if there is a next step
	      else {
		# No next step, so release the job and remove the current semaphore
		close JOBLOCK;
		unlink "$queue/$job.$listener_type";
	      }
	    }		  # end of listener managed to lock the job
	    sleep $listen_loop_delay;
	  }		  # end of for each job of this $listener_type
	}		  # end of for each $listener_type 

	#if ($jobs_found == 0) {	# if no jobs were found of any relevant $listener_type, sleep to save CPUs before we try again
	  sleep $listen_loop_delay; # varies by listener type
	#}
      } # end of listener loop - if we pass here it means a stopfile was found
      foreach my $listener_type (@listener_types) {
	if ($listener_type eq 'complete') {
	  # close the connection to the cache
	  my $cache_name = ${ $self->settings }{scheduler_cache_name};
	  $cache->disconnect($db, $cache_name);
	}
      }
      $timestamp = $GU->timestamp;

      $GU->user_info(2,"$timestamp : $listener_types listener told to stop listening to queue $queue by $username.\n");
      # ADD THIS INFO TO job_errors.log AS WELL:
      #my $APPLES_management_home = ${ $self->settings }{APPLES_management_home};
      open (REPORTERRORS, ">>$APPLES_management_home/job_errors.log") or die 'Cannot open job_errors.log\n'; #
      print REPORTERRORS "$timestamp : $listener_types listener told to stop listening to queue $queue by $username.\n";
      close (REPORTERRORS);
    } # listener managed to get lock on $queue for the relevant listener_types
    else {
      die "Cannot start $listener_types listener for queue $queue.\n";
    }
  } # listen #    
       
  method stop_listening(ListenerType $listener_type, Str $username) {
    # create file to signal listener to stop
    my $queue = ${$self->settings}{queue};
    open (STOPFILE, ">$queue/$listener_type.stop") or die "Cannot stop $listener_type listener in queue $queue.\n";
    close (STOPFILE);
    chmod (0664, "$queue/$listener_type.stop");
    my $grid_type = ${$self->settings}{grid_type};
    if ($grid_type eq 'PBS') {
      chmod ( 0666, "$queue/$listener_type.stop" );
    }
    # get username
    $GU->user_info(3,"$username stopped listener for $listener_type\n");
  } # stop_listening #
	
  method get_listener_status() {
    my $queue = ${$self->settings}{queue};

    # find lockfiles for any listeners
    my @allfiles;
    if (-d "$queue") {
      opendir(QUEUE, "$queue");
      @allfiles = readdir QUEUE;
      closedir(QUEUE);
    } else {
      die "Cannot find $queue.\n";
    }
    @allfiles = sort @allfiles;

    my %listeners;
    foreach my $lockfile (@allfiles) {
      if (substr($lockfile,-5) eq '.lock') { # file is a listener semaphore
	my $listener_type = substr($lockfile, (length $lockfile)-5); # extract the listener_type from the lockfile name
	if (open (QUEUELOCK, ">$queue/$lockfile")) { # if we are able to open lockfile for write access then the listener is not active
	  close (QUEUELOCK);
	  $listeners{$listener_type} = 'inactive';
	} else {
	  $listeners{$listener_type} = 'active';
	}
      }
    }
    return %listeners; # returns hash of listeners where key is listener_type and value is listener status
  } # get_listener_status #
	
  method is_job_in_queue(Str $method, Str $keyfile) {
    # my %jobs; my %other_cluster_jobs;
    # loop over listener types ->private_get_cluster_job_status():populates hashes
    # loop over jobs within %jobs and test for $method and $keyfile
    # return T/F
  } # is_job_in_queue #

  method private_get_cluster_job_status(HashRef $jobs_ref, HashRef $other_cluster_jobs_ref, ListenerType $listener_type, Boolean $refresh_hashes) {
      
      # Sascha: added parameter $refresh_hashes in an attempt to increase listeners' performance under high job load
      #         set to TRUE to restore orginal behaviour of this function
      
      my $currently_running_application_queue_jobs = 0;
      my $queue = ${$self->settings}{queue};
      my $job_parameters_file = ${$self->settings}{job_parameters_file};
      my $grid = ${$self->settings}{grid};
      my $cluster_dir = ${$self->settings}{cluster_dir};
      my $executable = ${$self->settings}{application}.'.queue.'.${$self->settings}{queue_id}; # executable gets assigned this name so it can be identified in qstat
      my $grid_type = ${$self->settings}{grid_type};

      %$other_cluster_jobs_ref = (); # moved out of if () clause below

      
      if ($refresh_hashes) {
	  # clear out %jobs and %other_cluster_jobs
	  %$jobs_ref = ();
	  
	  # look in $queue, find subdirectories (i.e. jobs) associated with $listener_type, and populate %jobs
	  my @allfiles = ();
	  
	  if (-d "$queue") {
	      opendir(JOBDIR, "$queue");
	      @allfiles = readdir JOBDIR;
	      closedir(JOBDIR);
	  } else {
	      die "Cannot find $queue.\n";
	  }
	  @allfiles = sort @allfiles;
	  
	  foreach my $subdir (@allfiles) { 
	      if ((-d "$queue/$subdir") && ($subdir ne '.') && ($subdir ne '..')) # it is a directory
	      {
		  # test whether this job belongs to $listener_type
		  if (-e "$queue/$subdir.$listener_type") {
		      # add the job to %jobs
		      $jobs_ref->{$subdir}{status} = $listener_type;
		      
		      # Capture the job status and add it to %jobs
		      my $job_parameters_filename = "$queue/$subdir/$job_parameters_file";		      
		      open(JOBPARAMS, $job_parameters_filename) or die "Cannot open ".$job_parameters_filename.".\n";
		      my @read_params = <JOBPARAMS>;
		      my @params;
		      my $i=0;
		      foreach my $param(@read_params) {
			  chomp ($read_params[$i]); # remove newline characters
			  if ($read_params[$i] ne '') {
			      push(@params,$read_params[$i]);
			  }
			  $i++;
		      }
		      close (JOBPARAMS);
		      $jobs_ref->{$subdir}{active} = FALSE;
		      $jobs_ref->{$subdir}{keyfile} = $params[0];
		      $jobs_ref->{$subdir}{calling_method} = $params[1]; 
		      $jobs_ref->{$subdir}{result_file} = $params[2];
		      $jobs_ref->{$subdir}{user_id} = $params[3];	
		      $jobs_ref->{$subdir}{queued_timestamp} = $params[4];
		      $jobs_ref->{$subdir}{cache_results} = $params[5];
		      $jobs_ref->{$subdir}{cache_duration} = $params[6];
		      $jobs_ref->{$subdir}{cluster_job} = '';	
		      $jobs_ref->{$subdir}{started_timestamp} = '';
		      $jobs_ref->{$subdir}{ended_timestamp} = '';

		      if (defined $params[7]) { # job has been started
			  my $offset; # found that occasionally there is an extra newline in the job parameters file
			              # detect it by checking where the time stamp is
			  if ($params[8] =~ m/_/) {
			      $offset = 0;
			  } else {
			      $offset = 1;
			  }
			  $jobs_ref->{$subdir}{cluster_job} = $params[7+$offset];	
			  $jobs_ref->{$subdir}{started_timestamp} = $params[8+$offset];
			  if (defined $params[9+$offset]) {
			      $jobs_ref->{$subdir}{ended_timestamp} = $params[9+$offset];
			  }
		      }
		  }
	      }
	  }
	  # %jobs now contains the details of all the jobs belonging to $listener_type

      } # if ($refresh_hashes)
      else {
	  foreach my $job (keys %{ $jobs_ref }) {
	      $jobs_ref->{$job}{qstat} = undef;
	      $jobs_ref->{$job}{active} = FALSE;
	  }
      }

      if ($grid) {
	  # make a temporary file to contain the qstat output
	  my ($temp_filehandle, $temp_filename);
	  ($temp_filehandle, $temp_filename) = tempfile();
	  close $temp_filehandle;
	  # do qstat - and make sure this command is completed
	  my $qstat_success = FALSE;
	  my $loop_counter = 0;
	  while (!$qstat_success) {
	      if ($loop_counter > 1000) {
		  die '1000 consecutive qstat attempts failed\n';
	      }
	      $loop_counter++;
	      my $exit_status = system ("qstat > $temp_filename"); # removed $cluster_dir/ from before qstat - or is it really needed?
	      if ($exit_status == 0) {
		  $qstat_success = TRUE;
	      }
	  }
	  # capture qstat output	
	  open(QSTATFILE, $temp_filename) or die "Cannot open $temp_filename.\n";
	  my @qstat = <QSTATFILE>;
	  close QSTATFILE;
	  unlink ($temp_filename);
	  # go through qstat output updating the %jobs and %other_cluster_jobs objects, and calculating $currently_running_application_queue_jobs
	  my $line_no = 0;
	  foreach my $line (@qstat) {
	      $line_no++;
	      if ($line_no > 2) {
		  chomp $line;
		  my @line = split / +/, $line;
		  
		  my $cluster_job_no = $line[0];
		  if ($grid_type eq 'PBS') {
		      $cluster_job_no = $line[0]; # value held in a different line in PBS qstat output
		      $cluster_job_no =~ s/\..*//; # remove .* tag (e.g. .fe1), as we're expecting numerical job numbers
		  }
		  # add active status to %jobs because job is still active
		  my $this_cluster_job_done = 0;
		  foreach my $job (keys %{ $jobs_ref }) {
		      if (($jobs_ref->{$job}{cluster_job} eq '')||
			  ($jobs_ref->{$job}{cluster_job} == 0)) {
			  if ($listener_type eq 'running') {
			      $jobs_ref->{$job}{active} = TRUE; # not yet managed to determine the queue number of this job, so assume it is running
			  }
		      } else {
			  if ($jobs_ref->{$job}{cluster_job} == $cluster_job_no) { # should "==" be "eq"?
			      $this_cluster_job_done = 1;
			      $jobs_ref->{$job}{qstat} = $line;
			      $jobs_ref->{$job}{active} = TRUE;
			  }
		      }
		  }
		  
		  # if the job isn't in our list of %jobs, add it to %other_cluster_jobs instead
		  if ($this_cluster_job_done == 0) {
		      $other_cluster_jobs_ref->{$cluster_job_no} = $line;
		  }

		  my $compare_to = $line[3];
		  if ($grid_type eq 'PBS') {
		      $compare_to = $line[1];
		  }
		  my $prefix = substr("$executable",0,length($compare_to));
		  if ($prefix eq $compare_to) { # job is owned by this application and this queue
		      $currently_running_application_queue_jobs++;
		  }	    
	      }			# end of line in qstat
	  }				# end of for each line in qstat 
      }				# end of if ($grid)
      return $currently_running_application_queue_jobs; # Int
  } # private_get_cluster_job_status #
  
  method private_add_job_information_line_to_shell_script(Job_Parameters $job_parameters) {
    # method to add a line to the shell script, e.g.
    #PBS -l select=1:ncpus=1:mem=1900mb,walltime=08:00:00
    #low=500mb, high=1900

    my $mem;
    if ($job_parameters->memory_requirement_high) {
      $mem = 1900;
    }
    else {
      $mem = 500; 
    }

    my $walltime = "48:00:00";
    if ($job_parameters->wall_time_estimate < 3600) {
      $walltime = "08:00:00"; # < 1hour estimate, -> 8 hour job walltime
    }
    if ($job_parameters->wall_time_estimate < 60) {
      $walltime = "00:01:00"; # < 60s estimate, -> 1 hour job walltime
    }
    my $line = "#PBS -l select=1:ncpus=1:mem=".$mem."mb,walltime=".$walltime."\n";
    
    return $line;
  } # private_add_job_information_line_to_shell_script #

  

} # Scheduler #
