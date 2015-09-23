### (c) copyright University of Warwick 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Running_Perseverance_Maker Class ###
# has methods that are default constructors for Running_Perseverance objects

use MooseX::Declare;

class Running_Perseverance_Maker {
  use Running_Perseverance;
  use constant {FALSE => 0,
                TRUE  => 1};

  method make_running_perseverance_try_hard() {
    # fill out hashrefs, with highly persistent settings
    my %max_attempts = ('Ensembl' => 10,
			'BiFa' => 10,
			'TempDirs' => 10,
			'CacheDB' => 10,
			'Unknown' => 10
		       );
    my %waiting_time = ('Ensembl' => 10,
			'BiFa' => 5,
			'TempDirs' => 10,
			'CacheDB' => 10,
			'Unknown' => 5
		       );
    my %infinite_tries = ('Ensembl' => FALSE,
			  'BiFa' => FALSE,
			  'TempDirs' => FALSE,
			  'CacheDB' => FALSE,			  
			  'Unknown' => FALSE
			 );
    my %try_counter = ('Ensembl' => 1,
		       'BiFa' => 1,
		       'TempDirs' => 1,
		       'CacheDB' => 1,
		       'Unknown' => 1
		      );
    
    my $running_perseverance = Running_Perseverance->new(
							 max_attempts => \%max_attempts,
							 waiting_time => \%waiting_time,
							 infinite_tries => \%infinite_tries,
							 try_counter => \%try_counter
							);
    return $running_perseverance;
  } # make_running_perseverance_try_hard #

  method make_running_perseverance_lightweight() {
  # fill out hashrefs, with minimum effort settings

    my %max_attempts = ('Ensembl' => 3,
			'BiFa' => 3,
			'TempDirs' => 3,
			'CacheDB' => 3,
			'Unknown' => 3
		       );
    my %waiting_time = ('Ensembl' => 1,
			'BiFa' => 1,
			'TempDirs' => 5,
			'CacheDB' => 5,
			'Unknown' => 1,			
		       );
    my %infinite_tries = ('Ensembl' => FALSE,
			  'BiFa' => FALSE,
			  'TempDirs' => FALSE,
			  'CacheDB' => FALSE,
			  'Unknown' => FALSE			 
			 );
    my %try_counter = ('Ensembl' => 1,
		       'BiFa' => 1,
		       'TempDirs' => 1,
		       'CacheDB' => 1,
		       'Unknown' => 1		       
		      );
    
    my $running_perseverance = Running_Perseverance->new(
							 max_attempts => \%max_attempts,
							 waiting_time => \%waiting_time,
							 infinite_tries => \%infinite_tries,
							 try_counter => \%try_counter
							);
    return $running_perseverance;
  } # make_running_perseverance_lightweight #

 method make_running_perseverance_try_infinitely_hard() {
    # fill out hashrefs, with highly persistent settings
    my %max_attempts = ('Ensembl' => 10,
			'BiFa' => 10,
			'TempDirs' => 10,
			'CacheDB' => 10,
			'Unknown' => 10
		       );
    my %waiting_time = ('Ensembl' => 5,
			'BiFa' => 5,
			'TempDirs' => 5,
			'CacheDB' => 5,
			'Unknown' => 5
		       );
    my %infinite_tries = ('Ensembl' => TRUE,
			  'BiFa' => TRUE,
			  'TempDirs' => TRUE,
			  'CacheDB' => TRUE,
			  'Unknown' => TRUE
			 );
    my %try_counter = ('Ensembl' => 1,
		       'BiFa' => 1,
		       'TempDirs' => 1,
		       'CacheDB' => 1,
		       'Unknown' => 1
		      );
    
    my $running_perseverance = Running_Perseverance->new(
							 max_attempts => \%max_attempts,
							 waiting_time => \%waiting_time,
							 infinite_tries => \%infinite_tries,
							 try_counter => \%try_counter
							);
    return $running_perseverance;
  } # make_running_perseverance_try_infinitely_hard #

   method make_running_perseverance_disabled() {
    my %max_attempts = ('Ensembl' => 0,
			'BiFa' => 0,
			'TempDirs' => 0,
			'CacheDB' => 0,
			'Unknown' => 0
		       );
    my %waiting_time = ('Ensembl' => 1,
			'BiFa' => 1,
			'TempDirs' => 1,
			'CacheDB' => 1,
			'Unknown' => 1,			
		       );
    my %infinite_tries = ('Ensembl' => FALSE,
			  'BiFa' => FALSE,
			  'TempDirs' => FALSE,
			  'CacheDB' => FALSE,
			  'Unknown' => FALSE			 
			 );
    my %try_counter = ('Ensembl' => 1,
		       'BiFa' => 1,
		       'TempDirs' => 1,
		       'CacheDB' => 1,
		       'Unknown' => 1		       
		      );
    
    my $running_perseverance = Running_Perseverance->new(
							 max_attempts => \%max_attempts,
							 waiting_time => \%waiting_time,
							 infinite_tries => \%infinite_tries,
							 try_counter => \%try_counter
							);
    return $running_perseverance;
  } # make_running_perseverance_disabled #

} # Running_Perseverance_Maker #
