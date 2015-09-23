#!/usr/bin/perl

=head1 Computation result cache

=cut

package Runtime::Cache;
use strict;
no warnings;

use CHI;
use DBI;
use Scalar::Util qw(blessed);
use JSON;
use Digest::MD5;
use Digest::SHA qw(sha512_hex);
use Data::Dumper;

use Runtime;
use Carp;

=head2 Constructor

 Parameters:
 $class : Runtime::Cache

 Returns
 a new Runtime::Cache

=cut

sub new {
	my $class = shift;

	my $conf = config;

	my $cache_driver = $conf->{'cache_driver'};

	my $cache_config_file = $conf->{'cache_config'};

	my $cache_conf = new Config::General($cache_config_file);
	my %all_caches = $cache_conf->getall();

	my $chi = undef;
	my $db  = undef;

	if ( $all_caches{disable_caching} ) {
		$chi = CHI->new( driver => 'RawMemory', );
	} else {
		my %cache_settings = %all_caches;

		foreach my $cc ( keys %all_caches ) {
			unless ( $cc =~ m/^db_/
				|| $cc eq 'disable_caching'
				|| $cc eq 'cache_writeonly'
				|| $cc eq 'cache_readonly' )
			{
				%cache_settings = %{ $all_caches{$cc} };
			}
		}
		
		if (lc ($cache_settings{driver}) eq "file") {
			## we need to change umask if we want the cache to be accessible 
			## to both www and other users
			umask 0;
			$chi = CHI->new(
			        driver         => 'File',
			        root_dir       => $cache_settings{directory},
			        depth          => 3,
			        max_key_length => 64,
					dir_create_mode => 0777,
					file_create_mode => 0666
			    );
		} elsif ( lc ( $cache_settings{driver} ) eq "dbi" ){
			my $dbi_string = "DBI:"
			  . "$cache_settings{db_engine}:"
			  . "database=$cache_settings{db_name};"
			  . "host=$cache_settings{db_server}";

			eval {

				$db = DBI->connect(
					$dbi_string,
					$cache_settings{db_username},
					$cache_settings{db_password},
					{ RaiseError => 1 }
				);

				$db->{mysql_auto_reconnect} = 1;

				$chi = CHI->new(
					driver   => 'DBI',
					dbh      => $db,
					l1_cache => {
						driver => 'RawMemory',
						global => 1,
					},
					compress_threshold => 2 * 1024*1024,
					create_table => 1,
					## FIXME Important: after creating the table, change
					## the field type to LONGTEXT, otherwise, the cache will
					## silently fail for everything that's bigger.
					##
					## This is a bug in CHI::Driver::DBI when using MySQL:
					## they use a TEXT field which is limited to 64kb
					##
					table_prefix => $cache_settings{'db_cache_table'},
				);

			};

			if ($@) {
				my $err = $@;
				$chi = CHI->new(
					driver => 'RawMemory',
					global => 1,
				);
				warn("Could not access MySQL ($err). Caching to memory.");
			}
			
		} else {
			die "Unknown cache driver : $cache_settings{driver}";
		}
	}

	return bless {
		'dbh' => $db,
		'chi' => $chi,
		'ro'  => $all_caches{cache_readonly} || 0,
		'wo'  => $all_caches{cache_writeonly} || 0,
	}, $class;
}

=head2 Make a cache id (private)

Parameter values can be marked to be not for caching
(i.e. serialized, but not included in the cache id)
using the prefix 'nc__'). Furthermore, if the parameter
class provides a method load_data ( $p ) [p is the
name of the key of the parameter], it can dynamically
load data just before we check the cache.

This provides a mechanism to change parameter
values just before checking whether a computation
result is in the cache. Use cases for this are:
 1) dynamic loading of input data from ResultDB
 2) random variables as input
 3) advanced parameters which we don't want the
    user to see, but which can influence computation
    outcomes

 Parameters:
 $self : a Runtime::Cache
 $type : a type
 $parameters : a Serialization::Serializable object
 $version : a version number (defaults to 1)

 Returns:
 a HASHREF like this : {
 	json => <JSON encoded id>,
 	hashstring => a hex-encoded hash value,
 }

=cut

sub _make_id {
	my $self       = shift;
	my $type       = shift;
	my $parameters = shift;
	my $version    = shift || 1;

	if ( UNIVERSAL::can( $parameters, 'load_data' ) ) {
		$parameters->load_data();
	}

	my $para = Serialization::Serializable::to_hash($parameters);

	## filter nc__ keys
	foreach my $k ( keys(%$para) ) {
		if ( $k =~ m/^nc__.*/
		|| ($k eq 'output_dataset' && UNIVERSAL::isa($parameters, 'Jobs::Job'))
		) {
			delete $para->{$k};
		}
	}

	#use Data::Dumper;
	#debug (Dumper ( $para ) );

	my $computation_id = {
		type       => $type,
		version    => $version,
		parameters => $para,
	};

	my $json_string =
	  to_json( $computation_id, { allow_blessed => 1, canonical => 1 } );

	my $digest = new Digest::MD5;
	$digest->add($json_string);
	my $hex_string = $digest->hexdigest();

	return {
		json       => $json_string,
		hashstring => $hex_string,
	};
}

=head2 Check if a cache id exists, and return its raw (serialized) record (private)

 Parameters:
 $self : a Runtime::Cache
 $cid : a cache id as returned by L<make_id>

 Returns:
 the record if the item exists, 0 otherwise

=cut

sub _raw_get_id {
	my $self = shift;
	my $cid  = shift;
	my $rv   = 0;

	return 0 if !defined($cid);

	eval {
		my $current_val = $self->{chi}->get( $cid->{hashstring} );

		if ( defined($current_val) ) {
			if ( $cid->{json} eq "*" ) {
				$rv = $current_val;
			}
			foreach my $result (@$current_val) {
				if ( $result->{input} eq $cid->{json} ) {
## This code can be used for debugging purposes.
## it transforms the output to JSON before it
## stores it.
## very slow.
					#					my $out_json = $result->{output};
					#					my $output   = from_json($out_json);
					$rv = $result;
					return;
				}
			}
		}
	};
	if ($@) {
		our $AR;
		warn("Cache test operation failed : $@");
		## Remove this hash key to prevent causing damage again
		$self->{chi}->remove( $cid->{hashstring} );
	}

	return $rv;
}

=head2 Get data from a cache id (private)

 Parameters:
 $self : a Runtime::Cache
 $cid : a cache id as returned by L<make_id>

 Returns:
 the cached data item for $cid, or undef

=cut

sub _get_id {
	my $self = shift;
	my $cid  = shift;
	my $rv   = undef;

	## if the cache is in write only mode, we don't return stuff here.
	if ( $self->{wo} != 0 ) {
		return undef;
	}

	eval {

		my $current_val = $self->{chi}->get( $cid->{hashstring} );
		
		# if (defined ($current_val)) {
		# 	debug ("Cache get FAILED for $cid->{hashstring} / $cid->{json}");			
		# } else {
		# 	debug ("Cache get SUCCESSFUL for $cid->{hashstring} / $cid->{json}");			
		# }

		if ( defined($current_val) ) {
			if ( $cid->{json} eq "*" ) {
				$rv = $current_val;
			}
			foreach my $result (@$current_val) {
				if ( $result->{input} eq $cid->{json} ) {
## This code can be used for debugging purposes.
## it transforms the output to JSON before it
## stores it.
## very slow.
					#					my $out_json = $result->{output};
					#					my $output   = from_json($out_json);
					my $output = $result->{output};
					$rv = Serialization::Serializable::from_hash($output);
					return;
				}
			}

			warn(   " HASH key matches, but JSON key doesn't: "
				  . Dumper($cid)
				  . Dumper($current_val) );
		}
	};
	if ($@) {
		our $AR;
		warn("Cache get operation failed : $@");
		## Remove this hash key to prevent causing damage again
		$self->{chi}->remove( $cid->{hashstring} );
	}

	return $rv;
}

=head2 Put data to a cache id (private)

 Parameters:
 $self : a Runtime::Cache
 $cid : a cache id as returned by L<make_id>
 $value : a Serialization::Serializable
 $expiry_time : a CHI expiry time (defaults to 'never')

 Returns:
 Nothing

=cut

sub _put_id {
	my $self        = shift;
	my $cid         = shift;
	my $value       = shift;
	my $expiry_time = shift || 'never';

	return if $expiry_time eq "0";

	## if the cache is in write only mode, we don't return stuff here.
	if ( $self->{ro} != 0 ) {
		return;
	}

	eval {
		# debug ("Cache put for $cid->{hashstring} / $cid->{json}");
		my $current_val = $self->{chi}->get( $cid->{hashstring} );

		if ( defined($current_val) ) {
			my $found = 0;

			## undef: remove the item
			if ( !defined($value) ) {
				@$current_val =
				  grep { $_->{input} ne $cid->{json} }
				  	@$current_val;
			} else {
				foreach my $result (@$current_val) {
					if ( $result->{input} eq $cid->{json} ) {
						$result->{output} =
						  Serialization::Serializable::to_hash($value);
## This code can be used for debugging purposes.
## it transforms the output to JSON before it
## stores it.
## very slow.
						#					$result->{output} = to_json(
						#						Serialization::Serializable::to_hash($value),
						#						{ allow_blessed => 1, canonical => 1 }
						#					);
						$found = 1;
						last;
					}
				}
				unless ($found) {
					push @$current_val, {
						input  => $cid->{json},
						output => Serialization::Serializable::to_hash($value),
## This code can be used for debugging purposes.
## it transforms the output to JSON before it
## stores it.
## very slow.
						#					output => to_json(
						#						Serialization::Serializable::to_hash($value),
						#						{ allow_blessed => 1, canonical => 1 }
						#					),
					};
				}

			}

		} elsif ( defined ($value) ) {
			$current_val = [
				{
					input  => $cid->{json},
					output => Serialization::Serializable::to_hash($value),
## This code can be used for debugging purposes.
## it transforms the output to JSON before it
## stores it.
## very slow.
					#					output => to_json(
					#						Serialization::Serializable::to_hash($value),
					#						{ allow_blessed => 1, canonical => 1 }
					#					),
				},
			];
		}
		if ( scalar @$current_val < 1 && defined ($value) ) {
			our $AR;
			confess( 'Could not put item to cache: ' . Dumper($cid) );
		}
		if ( scalar @$current_val > 1 ) {
			our $AR;
			warn( 'Cache hash collision for item ' . Dumper($cid) );
		}
		if (defined ($current_val)) {
			$self->{chi}->set( $cid->{hashstring}, $current_val, $expiry_time );
		}
	};
	if ($@) {
		our $AR;
		warn("Cache put operation failed : $@");
	}
}

=head2 Check if a cached result is available

Check if a computation result is available in the cache. If so,
return the result. If not, return undef if $run_if_na == 0, or
run the computation and cache+return it's result.

 Parameters:
 $self : a Runtime::Cache
 $comp : a computation

 Returns:
 1 if the result of the computation is in the cache,
 or 0 if not.

=cut

sub test_check_cache {
	my $self = shift;
	my $comp = shift;

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $comp, 'Runtime::Computation' ) );

	my $cid = $self->_make_id( blessed($comp), $comp, $comp->version() );
	return $self->_raw_get_id($cid);
}

=head2 Check if a cached result is available

Check if a computation result is available in the cache. If so,
return the result. If not, return undef if $run_if_na == 0, or
run the computation and cache+return it's result.

 Parameters:
 $self : a Runtime::Cache
 $comp : a computation
 $run_if_na : run the computation if the result is not available

 Returns:
 The result of the computation from the cache (if it's there),
 or undef if not.

=cut

sub check_cache {
	my $self      = shift;
	my $comp      = shift;
	my $run_if_na = shift || 1;

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $comp, 'Runtime::Computation' ) );

	my $cid = $self->_make_id( blessed($comp), $comp, $comp->version() );

	my $data = $self->_get_id($cid);

	if ( !defined($data)
		&& $run_if_na )
	{

		# run and store
		$data = $comp->run;

		$self->_put_id( $cid, $data, $comp->expiry_time() );
		$self->{chi}->set(
			"PROFILING_" . blessed($comp) . "_" . $comp->ns__started,
			{
				type    => blessed($comp),
				runtime => $comp->ns__runningtime,
				started => $comp->ns__started,
				cid     => $cid,
			},
			'3months'
		);
	}

	## cached or not, the computation object gets a chance here to
	## do some postprocessing
	if ( UNIVERSAL::can( $comp, 'postprocess' ) ) {
		$data = $comp->postprocess($data);
	}

	return $data;
}

=head2 Remove a computation from the cache

 Parameters:
 $self : a Runtime::Cache
 $comp : a computation

=cut

sub uncache {
	my $self = shift;
	my $comp = shift;

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $comp, 'Runtime::Computation' ) );

	my $cid = $self->_make_id( blessed($comp), $comp, $comp->version() );

	$self->_put_id( $cid, undef, $comp->expiry_time() );
}

=head2 Get a Serializable item from cache.

 Parameters:
 $self : a Runtime::Cache
 $parameters : the parameters (isa Serialization::Serializable)

 Returns:
 a Serialization::Serializable stored for $parameters, or undef

=cut

sub cache_get {
	my $self       = shift;
	my $parameters = shift;

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $parameters, 'Serialization::Serializable' ) );

	my $cid =
	  $self->_make_id( blessed($parameters), $parameters,
		$parameters->version );
	return $self->_get_id($cid);
}

=head2 Check if a Serializable item is in the cache.

 Parameters:
 $self : a Runtime::Cache
 $parameters : the parameters (isa Serialization::Serializable)

 Returns:
 1 if something is stored for $parameters, or 0 otherwise

=cut

sub cache_test {
	my $self       = shift;
	my $parameters = shift;

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $parameters, 'Serialization::Serializable' ) );

	my $cid =
	  $self->_make_id( blessed($parameters), $parameters,
		$parameters->version );
	return $self->_raw_get_id($cid);
}

=head2 Put a Serializable item to cache.

 Parameters:
 $self : a Runtime::Cache
 $parameters : the parameters (isa Serialization::Serializable)
 $value : the value to put (isa Serialization::Serializable)
 $expiry_time : a CHI expiry time (defaults to 'never')

 Returns:
 nothing

=cut

sub cache_put {
	my $self        = shift;
	my $parameters  = shift;
	my $value       = shift;
	my $expiry_time = shift || 'never';

	confess "Parameter check failed."
	  unless ( UNIVERSAL::isa( $self, 'Runtime::Cache' )
		&& UNIVERSAL::isa( $value,      'Serialization::Serializable' )
		&& UNIVERSAL::isa( $parameters, 'Serialization::Serializable' ) );

	my $cid =
	  $self->_make_id( blessed($parameters), $parameters,
		$parameters->version );
	$self->_put_id( $cid, $value, $expiry_time );
}

=head2 Get any type of item from cache using SHA hashing

 This function uses Storable to serialize. We don't check for
 double hash entries, so beware!

 Parameters:
 $self : a Runtime::Cache
 $key  : any Perl object that can be dumped using Data::Dumper

 Returns:
 The cached value.

=cut

sub cache_get_hash {
	my $self = shift;
	my $key  = shift;

	$key = "SHA:"
	  . sha512_hex( to_json( $key, { allow_blessed => 1, canonical => 1 } ) );

	my $val = $self->{chi}->get($key);

	return $val;
}

=head2 Put any type of item to cache using SHA hashing

 This function uses Storable to serialize. We don't check for
 double hash entries, so you need to check if what you get is actually
 valid!

 Parameters:
 $self  : a Runtime::Cache
 $key   : any Perl object that can be dumped using Data::Dumper
 $value : any Perl object that can be stored using Storable

 Returns:
 Nothing

=cut

sub cache_put_hash {
	my $self  = shift;
	my $key   = shift;
	my $value = shift;

	$key = "SHA:"
	  . sha512_hex( to_json( $key, { allow_blessed => 1, canonical => 1 } ) );

	$self->{chi}->set( $key, $value );
}

=head2 Dump all the cache contents

If there is lots of stuff in the cache, this will generate a lot of
output.

 Returns:
 A hashref with everything that's in the cache

=cut

sub get_cache_dump {
	my $self = shift;
	my $hr   = $self->{chi}->dump_as_hash();

	return $hr;
}

=head2 Purge expired stuff

=cut

sub purge_cache {
	my $self = shift;

	$self->{chi}->purge();
}

1;
