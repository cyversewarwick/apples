#!/usr/bin/perl

use MooseX::Declare;

=head1 Parameterized Computation Interface

This module specifies the interface for any computation for 
which we can cache the result based on an assignment of 
parameters.


=cut

class Runtime::Computation extends Datatypes::Datatype {
	use Runtime;
	use Carp;
	use POSIX qw(strftime);

	use Runtime::Cache;

## This is set to the date/time the computation was started
	has 'ns__started' => (
						   is      => 'rw',
						   isa     => 'Undef|Str',
						   default => undef,
	);

## This is set to the date/time the computation was started
	has 'ns__runningtime' => (
							   is      => 'rw',
							   isa     => 'Undef|Num',
							   default => 0,
	);

=head2 Reset the computation to runnable state

Needs to be called once a computation has been run, but failed
and wants to be run again.

=cut

	method reset () {
		$self->ns__runningtime(undef);
		$self->ns__started(undef);
	} 

=head2 This method will be run by the cache after the result has 
been retrieved

 Parameters:
 $result : the computation result (cached or just computed)

 Returns:
 a postprocessed version of the result

=cut

#	method postprocess ( Serialization::Serializable $result ) {
#		## by default, no action is necessary
#		return $result;
#	}

=head2  Run the computation

This method will run the _run method of the subclass, and 
validate the output. Furthermore, it will check the cache.

 Returns:
 
 Subclasses of this class need to return a Serializable object
 in _run, which is verified and returned here.

=cut

	method run () {
		my $result = undef;
		## this means our computation has been scheduled to run.
		## here, we see
		if ( !defined( $self->ns__started ) ) {
			$self->ns__started(
							"C " . strftime( '%d-%b-%Y %H:%M:%S', localtime ) );
			$result = cache->check_cache( $self, 1 );
		} elsif ( $self->ns__started =~ m/^C/ ) {
			## this means our result wasn't in the cache.
			## here, we see
			$self->ns__started(
							"R " . strftime( '%d-%b-%Y %H:%M:%S', localtime ) );
			$result = jobservice->run_computation($self);
		} else {
			## check if we have to load data
			if ( UNIVERSAL::can( $self, 'load_data' ) ) {
				$self->load_data();
			}
			if ( UNIVERSAL::can( $self, 'validate' ) ) {
				$self->validate();
			}

			use Time::HiRes qw(time);
			my $t0 = time;
			$result = $self->_run;
			$self->ns__runningtime( time() - $t0 );
		}
		if (    !defined($result)
			 || !UNIVERSAL::isa( $result, "Serialization::Serializable" ) )
		{
			confess(
"Computation result was invalid, must be defined and serializable." );
		}

		return $result;
	}

=head2 (virtual) Run the computation

 Parameters:
  
 
 Returns:
 
 Subclasses of this class need to return a Serializable object
 here.

=cut

	method _run() {
		confess("Virtual function called");
	}

=head2 Return the validity time for caching

By default, computational results never expire. 

You can overload this to change the default behaviour.
	
=cut

	method expiry_time() {
		return 'never';
	}

}
