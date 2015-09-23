### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Running_Perseverance Class ###
# class to handle failures e.g. errors caused by remote Ensembl DB connection (or lack thereof)
use MooseX::Declare;

class Running_Perseverance {
  use General_Utilities;
  use APPLES_Datatypes qw (Boolean PerseveranceCategory);
  use constant {FALSE => 0,
		TRUE  => 1};
  use Data::Dumper;

  use Moose;
  
  has 'max_attempts' => (is => 'rw', isa => 'HashRef', required => 1);
  has 'waiting_time' => (is => 'rw', isa => 'HashRef', required => 1); # seconds between tries
  has 'infinite_tries' => (is => 'rw', isa => 'HashRef', required => 1);
  has 'try_counter' => (is => 'rw', isa => 'HashRef', required => 1);
  has 'stop_trying' => (is => 'rw', isa => Boolean, required => 1, default => FALSE);
  
  my $GU = General_Utilities->new();
 
  method decide_on_rerun(PerseveranceCategory $perseverance_category, Boolean $die_on_job_information_exception, Any $error) { # this is $@
# case 1: no error -> stop_trying to TRUE, reset counter for this category
# case 2: an error, but it is Job_Information -> same behaviour as case 1, unless $die_on_job_information == TRUE
# case 3: below ...
    if ($error eq '') {
      $self->stop_trying (TRUE);
      ${$self->try_counter}{$perseverance_category} = 1;
      return;
    }
    elsif (!UNIVERSAL::can($error, 'isa')) {
      $GU->user_info(1, " *** global perseverance in use\n *** category: ". $perseverance_category. "\n");
      $GU->user_info(1, " *** try counter: ".${$self->try_counter}{$perseverance_category}."\n");
      ${$self->try_counter}{$perseverance_category} = ${$self->try_counter}{$perseverance_category} + 1;
      if (${$self->try_counter}{$perseverance_category} >= ${$self->max_attempts}{$perseverance_category}) {
	unless (${$self->infinite_tries}{$perseverance_category}) { 
	  $GU->user_info(1, "maximum attempts reached\n");
	  $self->stop_trying (TRUE);
	  die $error;
	}
      }
      else {
	sleep (${$self->waiting_time}{$perseverance_category});
	return;
      }

    }
    else {
      if (($error->isa('Job_Information_Exception'))||
	  ($error->isa('Aggregate_Exception'))){
	if ($die_on_job_information_exception) {
	  die ($error); # we have been asked to die to throw the Job_Information_Exception error
	}
	else {
	  ${$self->try_counter}{$perseverance_category} = 1; #  reset try counter
	  $self->stop_trying (TRUE);
	  return;
	}
      }
      else { # error is an object but not of Job_Information_Exception type 
	$GU->user_info(1, "global perseverance in use: ". $perseverance_category. " try counter: ".${$self->try_counter}{$perseverance_category}."\n");
	${$self->try_counter}{$perseverance_category} = ${$self->try_counter}{$perseverance_category} + 1;
	$GU->user_info(1, $perseverance_category. " try counter: ".${$self->try_counter}{$perseverance_category}."\n");
	if (${$self->try_counter}{$perseverance_category} >= ${$self->max_attempts}{$perseverance_category}) {
	  unless (${$self->infinite_tries}{$perseverance_category}) { 
	    $GU->user_info(1, "maximum attempts reached\n");
	    $self->stop_trying (TRUE);
	    die $error;
	  }
	}
	else {
	  sleep (${$self->waiting_time}{$perseverance_category});
	  return;
	}
      }
    }
  } # decide_on_rerun #

} # Running_Perseverance #
