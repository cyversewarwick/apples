#!/usr/bin/perl

=head1 Class Remote_Job_Service

Here we implement a job submission service that will connect to
ORANGES' jobservice.pl to submit jobs e.g. to a cluster.

=cut

package Runtime::Remote_Job_Service;

use strict;

our @ISA = qw(Runtime::Job_Service Utilities::JSON_Client);
require Runtime::Job_Service;
require Utilities::JSON_Client;

use Configuration::AppleSeeds;

use JSON;

use HTTP::Request;
use LWP::UserAgent;
use File::Temp qw(tempfile);
use Sys::Hostname;

require Runtime;
use Carp;

use Serialization::Serializable;
use Scalar::Util qw(blessed);

## for importing the job states
use Runtime::Web_Job_Service;

=head2 Constructor

 Parameters:
 $class : Runtime::Remote_Job_Service
 $url   : The url where jobservice.pl is running (we will append jobservice.pl?)
 $user  : The Authservice username
 $pass  : The Authservice password

 $huser : http auth user name (optional)
 $hpass : http auth password (optional)

 Returns:

 a new Remote_Job_Service

=cut

sub new {
	my ( $class, $url, $user, $pass, $huser, $hpass ) = @_;
	my $self = Utilities::JSON_Client->new( $url, $huser, $hpass, );

	$self->{user} = $user;
	$self->{pass} = $pass;

	bless $self, $class;
	return $self;
}

=head2 Returns true if we schedule computations
rather than running them straight away.

To enable this, include a line like this:

 our $use_scheduler = 1;

in your main perl script.
=cut

sub _use_scheduler {
	return defined($main::use_scheduler)
	  && $main::use_scheduler;
}

=head2 Run a computation

 Parameters:
 $self : a self object
 $comp : a Runtime::Computation

 Returns :
 This function will die
 	1) with a Runtime::Job_Queued_Exception
 	   if all went well
 	2) die using confess in any other case

=cut

sub run_computation {
	my $self = shift;
	my $comp = shift;

	if ( $self->_use_scheduler ) {
		my $json = to_json( Serialization::Serializable::to_hash($comp),
			{ allow_blessed => 1 } );

		my $result = $self->request(
			"jobservice.pl",
			[
				action => 'submit',
				job    => $json,
			]
		);

		die Runtime::Job_Queued_Exception->new( job => $comp );

	} else {
		$self->SUPER::run_computation($comp);
	}
}

=head2 Get the computations that have been queued.

 Parameters:
 $self : a self object

 Returns:
 A list of the jobs currently queued as follows :
 A list of the queued computations like this:
 [ {
 	id_comp => <a unique id>,
 	computation => the computation
 	state => the state the computation is in
 	         (see Runtime::Job_Service)
 },
 ...
 ]

=cut

sub get_queued_computations {
	my $self = shift;

	my $ret = $self->request("jobservice.pl?");

	if ( defined( $ret->{jobs} ) && ref( $ret->{jobs} ) eq "ARRAY" ) {
		my $result = [];

		foreach my $r ( @{ $ret->{jobs} } ) {
			eval {
				my $comp = {
					id_comp => int( $r->{id_comp} ),
					'state' => int( $r->{'state'} ),
					comp    => Serialization::Serializable::from_hash(
						from_json( $r->{json} )
					),
					host    => $r->{host},
					threads => int( $r->{threads} ),
					procs   => int( $r->{procs} ),
				};

				push @$result, $comp;
			};
			if ($@) {
				## if we can't deserialize it, we probably can't do anything
				## with it either. Others might though, so let's ignore it.
				Runtime::warn( "Ignoring invalid computation : " . $@ );
			}
		}
		return $result;
	} else {
		confess("Could not get queued computations from job service.");
	}
}

=head2 Get information about a queued/failed computation

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation

 Returns:
 Information about the queued computations like this:
 {
 	id_comp => <a unique id>,
 	computation => the computation,
 	output => the output
 	state => the state the computation is in
 }

=cut

sub get_computation_info {
	my $self    = shift;
	my $id_comp = shift;

	return $self->request(
		"jobservice.pl?",
		[
			action  => 'jobinfo',
			id_comp => $id_comp,
		]
	);
}

=head2 Remove a computation from the queue

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation

 Returns:
 Nothing if successful, dies otherwise

=cut

sub remove_computation {
	my $self    = shift;
	my $id_comp = shift;

	return $self->request(
		"jobservice.pl?",
		[
			action  => 'remove',
			id_comp => $id_comp,
		]
	);
}

=head2 Mark a computation as running

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation

 Returns:
 Nothing if successful, dies otherwise

=cut

sub mark_computation_running {
	my $self    = shift;
	my $id_comp = shift;

## jobs_hostid is necessary if we submit stuff via
## PBS, since the hostname might be different
## on different nodes.
## we use it to reset jobs that might have been
## aborted when a PBS job was killed.
	my $host = get_config_key('jobs_hostid') || hostname;

	return $self->request(
		"jobservice.pl?",
		[
			action  => 'markrunning',
			id_comp => $id_comp,
			host    => $host,
		]
	);
}

=head2 Mark a computation as failed

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation
 $message : a message

 Returns:
 Nothing if successful, dies otherwise

=cut

sub mark_computation_failed {
	my $self    = shift;
	my $id_comp = shift;
	my $message = shift;
	my $host    = get_config_key('jobs_hostid') || hostname;
	return $self->request(
		"jobservice.pl?",
		[
			action  => 'markfailed',
			id_comp => $id_comp,
			message => $message,
			host    => $host,
		]
	);
}

=head2 Mark a computation as complete

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation
 $message : a message

 Returns:
 Nothing if successful, dies otherwise

=cut

sub mark_computation_complete {
	my $self    = shift;
	my $id_comp = shift;
	my $message = shift;
	my $host    = get_config_key('jobs_hostid') || hostname;
	return $self->request(
		"jobservice.pl?",
		[
			action  => 'markcomplete',
			id_comp => $id_comp,
			message => $message,
			host    => $host,
		]
	);
}

=head2 Reset a computation

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation

 Returns:
 Nothing if successful, dies otherwise

=cut

sub reset_computation {
	my $self    = shift;
	my $id_comp = shift;
	my $host    = get_config_key('jobs_hostid') || hostname;
	return $self->request(
		"jobservice.pl?",
		[
			action  => 'reset',
			id_comp => $id_comp,
			host    => $host,
		]
	);
}

=head2 Get a data item from the result database

Jobs can retrieve data items from a result database. These data items
are associated with a dataset ID.

 Parameters:
 $self    : a self object
 $dataset : a dataset name (folder)
 $source  : the dataitem source name
 $version : (optional) the version of this object
            by default, we return the latest version
            version can be an identifier returned by
            get_dataitem_versions

 Returns:
 The latest version of that dataitem as {
		type => $type,
		data => $data,
 }

 If the result type is JSON, the data will be decoded.
=cut

sub get_dataitem {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;

	my $type;
	my $version = shift || "";

	my $didata = undef;
	if ( $version ne "" ) {
		$didata = $self->request(
			"dataservice.pl?",
			[
				action  => 'get_versioned_dataitem',
				dataset => $dataset,
				source  => $source,
				version => $version,
			]
		);
	} else {
		$didata = $self->request(
			"dataservice.pl?",
			[
				action  => 'get_versioned_dataitem',
				dataset => $dataset,
				source  => $source,
			]
		);
	}

	## check if we downloaded something not JSON
	if ( defined( $didata->{NONJSON} ) ) {
		$didata = $didata->{NONJSON};
	} else {
		## otherwise, JSON_Client will have decoded the JSON stuff and returned it as a HASHREF
		$didata = {
			type => "json",
			data => $didata,
		};
	}
	return $didata;
}

=head2 Put a data item to the result database

Jobs can retrieve data items from a result database. These data items
are associated with a dataset ID. Note that all conversions must be
carried out before calling this (like converting to JSON/etc). We will
expect a scalar item, which will be written as it is.

Dataset and source name will be transformed to make safe file names,
allowed characters are [A-Za-z0-9\._]

In this implementation, we will use a flat file system database for this.

 Parameters:
 $self    : a self object
 $dataset : a dataset name (folder)
 $source  : the dataitem source name
 $item    : the data item
 $type    : (optional) the type to use for storage (default : JSON)

 Returns:
 nothing, but might die over serialization

=cut

sub put_dataitem {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;
	my $item    = shift;
	my $type    = shift || 'JSON';

	## make hashes into JSON
	if ( $type =~ m/JSON/i ) {
		$type = "application/json";
		if ( Serialization::Serializable::is_hash($item) ) {
			$item = to_json( $item, { allow_blessed => 1 } );
		}
	}

	my ( $fh, $filename ) = tempfile();

	syswrite $fh, $item;
	close($fh);

	my $res = $self->request(
		"dataservice.pl?",
		[
			action  => "upload",
			item    => [$filename],
			dataset => $dataset,
			type    => $type,
			source  => $source,
		]
	);
}

=head2 Get all versions of a data item

 Parameters:
 $self    : a self object
 $dataset : a dataset name (folder)
 $source  : the dataitem source name

 Returns:
 an ARRAYREF with version identifiers

=cut

sub get_dataitem_versions {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;

	return $self->request(
		"dataservice.pl",
		[
			action  => 'listversions',
			dataset => $dataset,
			source  => $source,
		]
	);
}

=head2 Purge old versions of a data item

This function will remove versions of a data item.

When not specifying versions, use with caution:
If multiple jobs write to the same item, they might need
to be merged before calling this function, otherwise, data could be lost.
(not likely with the debug version of the job service, but possible
when jobs are run on a remote job service)

 Parameters:
 $self      : a self object
 $dataset   : a dataset name (folder)
 $source    : the dataitem source name
 $rversions : (optional) an ARRAYREF containing the versions to remove
              defaults to all but the most recent version

 Returns:
 nothing

=cut

sub purge_dataitem_versions {
	my $self      = shift;
	my $dataset   = shift;
	my $source    = shift;
	my $rversions = shift || [];

	return $self->request(
		"dataservice.pl",
		[
			action   => 'purgeversions',
			dataset  => $dataset,
			source   => $source,
			versions => join( ':', @$rversions ),
		]
	);
}

=head2 Perform JSON request, authenticate if necessary

 Parameters:
 $self : a self object
 $urn  : the service name to access
 $parameters : the POST parameters

 Returns:
 a translated HASHREF from the JSON call

=cut

sub request {
	my $self       = shift;
	my $urn        = shift;
	my $parameters = shift;

	if ( defined( $self->{session} ) ) {
		push @$parameters, ( CGISESSID => $self->{session} );
	}

	my $result = $self->SUPER::request( $urn, $parameters );

	## access denied? authenticate and try again
	if ( defined($result)
		&& Serialization::Serializable::is_hash($result) )
	{

		if ( defined( $result->{access} )
			&& $result->{access} eq 'denied' )
		{
			$self->_authenticate;
			$result = $self->SUPER::request( $urn, $parameters );
		}

		if ( defined( $result->{error} ) ) {
			confess("Could connect to job service: $result->{error}.");
		}
	}
	return $result;
}

=head2 Connect to authservice.pl and get a valid session

 Parameters:
 $self : a self object

 Returns:
 Nothing, but will die if access is denied

=cut

sub _authenticate {
	my $self = shift;

	my $result = $self->request(
		"authservice.pl",
		[
			action   => 'login',
			user     => $self->{user},
			password => $self->{pass},
		]
	);

	if ( $result->{access} ne "granted" ) {
		confess(
			"Could not log in to JSON authservice : " . $result->{access} );
	} else {
		$self->{session} = $result->{session};
	}
}

=head2 Get job type buckets

 Parameters:
 $self : a self object

 Returns:
 A bucket hash
=cut

sub _buckets {
	my $self = shift;

	return $self->request(
		"jobservice.pl",
		[ "action" => "parallelism_buckets" ]
	);
}

1;
