### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Cache Class ###

# Purpose:  Implements interface with database of previously-calculated results (the cache) 
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

# each job in the cache has:
# $calling_object - MD5 hash of file containing serialised calling object
# $method - name of calling object method
# $keyfile - MD5 hash of file containing parameters and optionally the calling object
# $user_id - user submitting the job
# $queued_timestamp - timestamp the job was queued
# $started_timestamp - timestamp the job was started
# $ended_timestamp - timestamp the job finished
# $cache_duration - time to store the result (days?)
# $results - large text object containing serialised results

use MooseX::Declare;

class Cache {
  use APPLES_Datatypes qw (Status);
  use General_Utilities;
  use Running_Perseverance;
  use DBI;
  use Data::Dumper;
  use File::Temp qw (tempdir);
  use File::Path qw (rmtree);
  use Cwd;
  use Digest::MD5 qw(md5_hex);
  use constant {FALSE => 0,
		TRUE  => 1};

  my $GU = General_Utilities->new;
  my $md5_generator = Digest::MD5->new();

  #my $db_engine; ### Name of database engine read from config (mysql)
  #my $db_name; ### Name of database engine read from config (apples_cache_dev)
  #my $db_server; ### Name of server hosting database read from config (www.wsbc.warwick.ac.uk)
  #my $db_username; ### Login name for database read from config (ensembl)
  #my $db_password; ### Password for database read from config (ens3mbl)
  #my $db_cache_table; ### Database table containing cache items read from config (plants)
  
  has 'cache_settings' => (is => 'rw', isa => 'HashRef');#, required => 1);# obtained from Job_Handler, usage: ${$self->cache_settings}{}
  
  ### permissions:  read/write, ownership (user, group)
  ### MD5 here
  
  method connect (Str $cache_name) {

    # Connect to a database and return $db (db handle)
    my %all_cache_settings = %{$self->cache_settings};
    my %cache_settings = %{$all_cache_settings{$cache_name}};
   
    my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";
    my $db;
    while (!$main::global_perseverance->stop_trying) {
	eval {
	    $db = DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } );
  	      # DBI source specification, Username and password, Additional options
	};
	
	
	# This should resolve bug with linne 139, by reconnecting when connection timesout.
	$db->{mysql_auto_reconnect} = 1;
	
	my $error = $@;
	$main::global_perseverance->decide_on_rerun('CacheDB', TRUE, $error);
    };
    $main::global_perseverance->stop_trying(FALSE);
    return $db;
  } # connect #

  method connect_all_caches (ArrayRef[Str] $caches_ref) {
    my %db_handles;
    foreach my $cache_name (@{ $caches_ref }) {
      $db_handles{$cache_name} = $self->connect($cache_name);
    }
    return %db_handles; #return hash of cache db_handles
  } # connect_all_caches #
  
  
  method disconnect (Any $db, Str $cache_name) { # stricter datatype? Object?
    # Disconnecting from a database
    # reconnect if we have lost the connection
    my %all_cache_settings = %{$self->cache_settings};
    my %cache_settings = %{$all_cache_settings{$cache_name}};
   
    my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";

    #$db ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options

    $db->disconnect or die "Could not disconnect cache database.\n";
  } # disconnect #
    
  method store_result(Any $db_handle, Str $cache_name, Str $keyfile, Str $method, Str $user_id, Str $queued_timestamp, Str $started_timestamp, Str $ended_timestamp, Int $cache_duration, Str $result, Str $tablename) { 
    ### Str $result is the gzip filepath 
    $GU->user_info(3,"ATTEMPTING TO STORE RESULT $result\n");
    my $success = 0;
    open(RESULTFILE, $result) or die "Cannot open '$result': $!";
    binmode(RESULTFILE);
   
    my $resultstring = do { local $/; <RESULTFILE> };
    close(RESULTFILE);
    $resultstring =~ s/\\/\\\\/ogis; #escapes backslashes
    $resultstring =~ s/\'/\\\'/ogis; #escapes single quotes
    $resultstring =~ s/\"/\\\"/ogis; #escapes double quotes
    $resultstring =~ s/\0/\\\0/ogis; #escapes null character
    ### check permissions?	  

    ### NEED TO FIDDLE ABOUT WITH THIS BIT TO GET IT TO WORK!! -V:2
    open(FILE, $result) or die "Can't open '$result': $!";
    binmode(FILE);
    my $new_resultMD5 = $md5_generator->addfile(*FILE)->hexdigest or die "Cannot make md5 of $result.\n"; # obtain MD5 of $result
    close(FILE);
    $GU->user_info(3,"new result MD5: ".$new_resultMD5."\n");
    # get results from the cache or return null/error
    my $cached_result;
    my $no_rows = 0;
    # reconnect if we have lost the connection
    my %all_cache_settings = %{$self->cache_settings};

    my %cache_settings = %{$all_cache_settings{$cache_name}};
   
    my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";
    
    # THE OPTIONAL RECONNECT BELOW DOES NOT PREVENT THE FOLLOWING CRASH IN THE COMPLETE LISTENER:
    # DBD::mysql::st execute failed: MySQL server has gone away at /cluster/laurabaxter/APPLES_WORKING_COPY/Main//Cache.pm line 126
    #$db_handle ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options
	$db_handle = $self->connect($cache_name);
    for my $row ($db_handle->query("SELECT results FROM $tablename WHERE keyfile='$keyfile' AND method='$method';")->hashes) {
      $cached_result = $row->{results};  ### what about timestamps and all the rest?
      $no_rows++;      
    }
    my $cached_resultMD5;
    if ($no_rows == 1) {
      $GU->user_info(3,"$no_rows row found in cache\n");  
      # $cached_result was found so check md5 of $cached_result against md5 of $result
      # get MD5 of cached result
      my $tempdir = tempdir();
      my $tempfile = "tempfile.gz"; # file to write zip to
      open (ZIPFILE, ">$tempdir/$tempfile") or die "Cannot create $tempdir/$tempfile.\n";
      $GU->user_info(3,"zipped file temporarily saved into $tempdir/$tempfile\n");
      binmode (ZIPFILE);
      print ZIPFILE $cached_result;
      close (ZIPFILE);
      open(FILE, "$tempdir/$tempfile") or die "Can't open '$tempdir/$tempfile': $!";
      binmode(FILE);
      
      $cached_resultMD5 = $md5_generator->addfile(*FILE)->hexdigest or die "Cannot make md5 of $tempdir/$cached_result.\n"; # obtain MD5 of cached_result
      close(FILE);
      $GU->user_info(3,"cached result MD5: ". $cached_resultMD5."\n");
      #my $cwd = getcwd();
      #chdir $tempdir;
      #chdir $cwd;
      
      
      if ($new_resultMD5 eq $cached_resultMD5) {
	# cached result is fine (md5 matches) so leave it and move on
	$success = 1;
	$GU->user_info(3,"MD5 of both results is the same (keep original result in cache)\n");
	unlink "$tempdir/$tempfile";
	rmtree("$tempdir");
      }
      else {
	$GU->user_info(1,"Result was previously cached, but new result does not match cached one!\n");
      }	  
    }
    else {
      # result has not been previously cached so cache it
      ### include last time retrieved, total number of times retrieved, size of cached result file??
      my $values = "'$method', '$keyfile', '$user_id', '$queued_timestamp', '$started_timestamp', '$ended_timestamp', '$cache_duration'"; 
      my %all_cache_settings = %{$self->cache_settings};

      my %cache_settings = %{$all_cache_settings{$cache_name}};

      my $tabletoinsert = $cache_settings{db_cache_table};
      $GU->user_info(3,"values: $values\n");
      $GU->user_info(3,"$tabletoinsert\n");
      # reconnect if we have lost the connection
   
      my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";

      #$db_handle ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options
	
		$db_handle = $self->connect($cache_name);
		
      $db_handle->query("INSERT INTO $tabletoinsert (method, keyfile, user_id, queued_timestamp, started_timestamp, ended_timestamp, cache_duration, results) VALUES ($values, '$resultstring');") or die $db_handle->error;
      $success = 1;
    }	  
    return $success;
  } # store_result #
  
  method retrieve_result(HashRef $db_handles, Str $keyfile, Str $method, Str $user_id) {
	  $GU->user_info(3, "checking cache for result (".$method.")\n");

	  ### got to here!! do we need to 'connect_all_caches'?

	  #die "method not implemented yet!";
	  # $db is a hash ref of db handles
	  ### check permissions?	  
	  ### check each cache (in the HashRef of db that gets passed in) for the result? (loop here or outside?)
          ### $keyfile is an MD5 sum

	  ### include last time retrieved, total number of times retrieved, size of cached result file??
	  # need to dereference the $db HashRef and check each for the result, or only pass in one at a time
          # get results from the cache or return null('NULL')/error 
	  my $this_cache;
	  my %all_cache_settings = %{$self->cache_settings};

	  foreach my $cache (keys %{ $db_handles }) {
	    my %cache_settings = %{$all_cache_settings{$cache}};
	    if ($self->check_cache($db_handles->{$cache}, $cache, $keyfile, $method, $user_id, $cache_settings{db_cache_table}) > 0) {
	      $this_cache = $cache;
	    }
	  }
	  if (defined $this_cache) {
	    
	    my %cache_settings = %{$all_cache_settings{$this_cache}};
	    my $result;
	    my $no_rows = 0;
		  #$GU->user_info(3,Dumper($keyfile));
		  #$GU->user_info(3,($method));

	    # reconnect if we have lost the connection
	    
	    #$db_handles->{$this_cache} ||= DBIx::Simple->connect("DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}", $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to database.\n";  # DBI source specification, Username and password, Additional options

	    for my $row ($db_handles->{$this_cache}->query("SELECT results FROM $cache_settings{db_cache_table} WHERE keyfile='$keyfile' AND method='$method';")->hashes ) { #or die $db_handles->{$this_cache}->error
	      $result = $row->{results};  ### what about timestamps and all the rest?
	      $no_rows++;      
	    } 
	    if ($no_rows == 1) {
	    
	      my $tempdir = tempdir();
	      #$GU->user_info(3,$tempdir."\n");
	      my $tempfile = "tempfile.gz"; # to write zip to
	      open (ZIPFILE, ">$tempdir/$tempfile") or die "Cannot create $tempdir/$tempfile.\n";
	      binmode (ZIPFILE);
	      print ZIPFILE $result;
	      close (ZIPFILE);
	      my $cwd = getcwd();
	      chdir $tempdir;
	      if ($^O eq 'MSWin32'){
	      	 my $intermediateFile = substr($tempfile,0,-3); 
		     system ("gunzip $tempfile") == 0 || die "System error!";
		     system ("tar -xf $intermediateFile") == 0 || die "System error!";
	      }else {	
		  	$GU->user_info(3,"trying to extract: $tempdir/$tempfile\n");
		  	my $return_code = system ("tar -zxf $tempfile 2> noisy.txt"); # avoid redundant warnings
		  	if ($return_code != 0 ) {
		      	if ($return_code != 512 ) {  # occurs frequently, but all cases checked so far were fine
			  		$GU->user_info(3,"return code: ".$return_code);
			  		die 'System error!';
		      	}
		  	}
	      }
	      chdir $cwd;

	      unlink "$tempdir/$tempfile";
	      return $tempdir;
	    }
	    else {
	      die "$no_rows were found when 1 was expected.\n"; 
	    }
	  }
	  else {
	    return 'NULL';
	  }
	} # retrieve_result #

  method delete_results(HashRef $db_handles, Str $keyfile, Str $method, Str $user_id) {
	  $GU->user_info(3, "deleteing result in cache\n");
	  my $this_cache;
	  my %all_cache_settings = %{$self->cache_settings};

	  foreach my $cache (keys %{ $db_handles }) {
	      my %cache_settings = %{$all_cache_settings{$cache}};
	      $self->delete_result($db_handles->{$cache}, $keyfile, $method, $user_id, $cache_settings{db_cache_table});
	  }
  } # delete_results #

  method delete_result(Any $db, Str $keyfile, Str $method, Str $user_id, Str $tablename) {
          my $success = 0;

	  ### check permissions?

      ### take md5 of $calling_object_file and $parameters_file
	  #my $calling_object = md5($calling_object_file);
	  #$keyfile = md5($keyfile); # MD5 already passed through!
	  # $user_id is not used - doesn't need to be passed through
	  # delete results from the cache
	  $db->query("DELETE FROM $tablename WHERE keyfile='$keyfile' AND method='$method';") or die $db->error;
	  $success = 1;
	  return $success;
	} # delete_result #

  method check_all_caches () {
    die 'method not implemented yet';
  } # check_all_caches #

        method check_cache(Any $db_handle, Str $cache_name, Str $keyfile, Str $method, Str $user_id, Str $tablename) {
          my $success = 'null';
	  #$GU->user_info(3,Dumper ($db_handle));
	  ### check permissions?
          ### should we do the MD5 here?
	  # $db is a db handle
          ### take md5 of $calling_object_file and $parameters_file
	  #my $calling_object = md5($calling_object_file);
	  #$keyfile = md5($keyfile); #  MD5 already passed through!
	  
	  # query cache for existence of matching result ### NOT WORKING?!
          my $result;
	  my $no_rows = 0;

	  # reconnect if we have lost the connection
	  my %all_cache_settings = %{$self->cache_settings};
	  my %cache_settings = %{$all_cache_settings{$cache_name}};
   
	  my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";

	  #$db_handle ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options
	  $db_handle = $self->connect($cache_name);
          for my $row ($db_handle->query("SELECT COUNT(*) no_cached_results FROM $tablename WHERE method='$method' AND keyfile='$keyfile';")->hashes ) {#or die $db_handle->error
            $result = $row->{no_cached_results}; # result is one row containing the number of matching results
            $no_rows++;      
          }
	  if ($no_rows == 1) {
	    return $result;
	  }
	  else {
	    die "$no_rows were found when 1 was expected.\n";
	  }
	} # check_cache #

        method tidy_cache_age (Str $db, Str $oldest_to_keep, Str $user_id) {
          die 'method not implemented yet';
	  my $success = 0;
	  ### check permissions?
          ### delete results from the cache where ended_timestamp < $oldest_to_keep  
	  return $success;
	} # tidy_cache_age #

        method tidy_cache_unused(Str $db, Str $oldest_to_keep, Str $user_id) {
          die 'method not implemented yet';
	  my $success = 0;
	  ### check permissions?
          ### delete results from the cache where last_retrieved > $oldest_to-keep
	  return $success;
	} # tidy_cache_unused #

        method tidy_cache_method(Str $db, Str $method, Str $user_id) {
          die 'method not implemented yet';
          my $success = 0;
	  ### check permissions?
          ### delete results from the cache where method = $method
	  return $success;
	} # tidy_cache_method # 

        method tidy_cache_expired(Str $db, Str $user_id) {
	  die 'method not implemented yet';
          my $success = 0;
	  ### check permissions?
          ### delete results from the cache where insertion timestamp < $GU->timestamp - cache_duration	  
	  return $success;
	} # tidy_cache_expired #

        method reset_cache(Any $db, Str $user_id, Str $cache_name) {
	  ### THIS FUNCTION WILL DELETE ENTIRE CONTENTS OF CACHE DATABASE: USE WITH EXTREME CAUTION! ###
          my $success = 0;
	  my %all_cache_settings = %{$self->cache_settings};
	  my %cache_settings = %{$all_cache_settings{$cache_name}};
	  ### check permissions VERY CAREFULLY?
	  if ($user_id ne 'laurabaxter') {
	    die 'you do not have sufficient privileges to reset the cache!';
	  }
	  else {
	    $GU->user_info(3,"resetting cache\n");
	  }
	  # NB need to be connected to $db first
          # Destroy the cache table and recreate it empty 

	  # reconnect if we have lost the connection
   
	  my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";

	  #$db ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options
	$db = $self->connect($cache_name);
	
	  $db->query("DROP TABLE IF EXISTS $cache_settings{db_cache_table}") or die $db->error;
	  $db->query("CREATE TABLE $cache_settings{db_cache_table} (keyfile VARCHAR(255), method VARCHAR(255), user_id VARCHAR(255), queued_timestamp VARCHAR(255), started_timestamp VARCHAR(255), ended_timestamp VARCHAR(255), cache_duration INT, results LONGBLOB)") or die $db->error;
	  $db->query("ALTER TABLE $cache_settings{db_cache_table} ADD PRIMARY KEY (keyfile, method)") or die $db->error;
          $success = 1;
	  return $success;
	} # reset_cache #

  method build_cache (Any $db, Str $user_id, Str $cache_name) {
    die 'method not tested yet';
    my $success = 0;
    my %all_cache_settings = %{$self->cache_settings};
    my %cache_settings = %{$all_cache_settings{$cache_name}};
    ### check permissions VERY CAREFULLY?
    if ($user_id ne 'laurabaxter') { # simple check for now, edit as necessary
      die 'you do not have sufficient privileges to build a cache!';
    }
    else {
      $GU->user_info(3,"building cache\n");
    }
    # NB need to be connected to $db first
    # Build empty cache 
    # reconnect if we have lost the connection
   
    my $dbi_string = "DBI:$cache_settings{db_engine}:database=$cache_settings{db_name};host=$cache_settings{db_server}";

    #$db ||= DBIx::Simple->connect($dbi_string, $cache_settings{db_username}, $cache_settings{db_password}, { RaiseError => 1 } ) or die "Could not connect to $dbi_string.\n";  # DBI source specification, Username and password, Additional options
	$db = $self->connect($cache_name);

    $db->query("CREATE TABLE $cache_settings{db_cache_table} (keyfile VARCHAR(255), method VARCHAR(255), user_id VARCHAR(255), queued_timestamp VARCHAR(255), started_timestamp VARCHAR(255), ended_timestamp VARCHAR(255), cache_duration INT, results LONGBLOB)") or die $db->error;
    $db->query("ALTER TABLE $cache_settings{db_cache_table} ADD PRIMARY KEY (keyfile, method)") or die $db->error;
    $success = 1;
    return $success;
  } # build_cache #

} # Cache # 
