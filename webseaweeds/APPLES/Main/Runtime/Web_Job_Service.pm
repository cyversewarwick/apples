#!/usr/bin/perl

=head1 ORANGES Web job service

Everything that is needed in jobservice.pl to implement job handling.

=cut

package Runtime::Web_Job_Service;

use strict;
use feature ':5.10';

our @ISA = qw(Exporter);
require Exporter;
our @EXPORT = qw(STATE_QUEUED STATE_RUNNING STATE_FAILED STATE_COMPLETE);

use Runtime;
use Configuration::ReadConfig qw($basedir);
use Configuration::AppleSeeds;

use Sys::Hostname;
use JSON;

use Scheduler::PBS;

use Runtime::Web_Session_Helpers;

=head1 Job Codes

=cut

use constant {
	STATE_QUEUED   => 0,
	STATE_RUNNING  => 1,
	STATE_FAILED   => 100,
	STATE_COMPLETE => 200,
};

=head1 Interface compatible to Runtime::Remote_Job_Service

=cut

=head2 Constructor

Construct a web job service from the current session

 Parameters:
 $class : Runtime::Web_Job_Service
 $userid : the user id (defaults to the one in the current session)

 Returns:
 a new Runtime::Web_Job_Service

=cut

sub new {
	my $class  = shift;
	my $userid = shift;

	if ( !defined($userid) ) {
		my ( $cgi, $session ) = get_cgi_and_session;
		$userid = $session->param('userid');
	}

	my $self = {};
	$self->{userid} = $userid;

	bless $self, $class;
	return $self;
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

	my $use_scheduler = defined($main::use_scheduler)
	  && $main::use_scheduler;

	if ($use_scheduler) {
		my $json = "";
		if ( Serialization::Serializable::is_hash($comp) ) {
			$json = to_json( Serialization::Serializable::to_hash($comp),
				{ allow_blessed => 1 } );
			my $result = $self->_submit_job($json);
			if ($result->{submitted} == 1) {
				die Runtime::Job_Queued_Exception->new( job => $comp );
			} else {
				if ($result->{'state'} == STATE_QUEUED) {
					die Runtime::Job_Queued_Exception->new( job => $comp );
				} else {
					if ($result->{'state'} == STATE_COMPLETE) {
						my $result = cache->check_cache($comp, 0);
						if (!defined ($result)) {
							use Data::Dumper;
							die "Cached result for $result->{id_comp} not available; is the cache in write-only mode? " . Dumper ($comp);
						}
						return $result;
					} elsif ($result->{'state'} == STATE_FAILED) {
						use Data::Dumper;
						die "This job has failed. To re-run, you need to re-submit the job ($result->{id_comp}): " . Dumper ($comp);
					} else {
						$self->reset_computation($result->{id_comp});
						warn ("Computation $result->{id_comp} re-submitted but is in state $result->{state}. Computation was reset.");						
					}
				}
			}
		} else {
			$json = $comp;
			return $self->_submit_job($json);
		}
	} else {
		if ( !UNIVERSAL::isa( $comp, "Runtime::Computation" ) ) {
			$comp =
			  Serialization::Serializable::from_hash( from_json($comp) )
			  ;
		}
		if ( !UNIVERSAL::isa( $comp, "Runtime::Computation" ) ) {
			die "Invalid computation: $comp";
		}
		return $comp->run();
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

	my $ret = $self->_list_jobs();

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

=head2 Get job information

 Parameters:
 $self     : a self object
 $id       : the id
 $result   : (optional) hashref to append results to

 Returns:
 $result = {
 	killed => number of jobs killed
 }
=cut

sub get_computation_info {
	my $self   = shift;
	my $userid = $self->{userid};
	my $id     = shift;
	my $result = shift || {};

	my $dbh = get_dbh;

	my $sth = $dbh->prepare(
		"SELECT * FROM `computations` WHERE `id_user` = ? AND `id_comp` = ?")
	  or die "Cannot access the computation database.";
	$sth->execute( $userid, $id )
	  or die $sth->errstr;

	if ( my $ref = $sth->fetchrow_hashref() ) {
		$result->{jobinfo} = {
			id_comp   => $ref->{id_comp},
			json      => $ref->{json},
			'state'   => $ref->{'state'},
			'output'  => $ref->{'output'},
			'procs'   => $ref->{'procs'},
			'threads' => $ref->{'threads'},
			'output'  => $ref->{'output'},
		};
	} else {
		die "Computation not found.";
	}
	return $result;
}

=head2 Remove a computation from the table

 Parameters:
 $self     : a self object
 $id       : the id to remove
 $result   : (optional) hashref to append results to

 Returns:
 $result = {
 	killed => number of jobs killed
 }
=cut

sub remove_computation {
	my $self   = shift;
	my $userid = $self->{userid};
	my $id     = shift;
	my $result = shift || {};

	my $dbh = get_dbh;
	my $sth = $dbh->prepare(
		"DELETE FROM `computations` WHERE `id_user` = ? AND `id_comp` = ?")
	  or die "Cannot access the computation database.";
	$sth->execute( $userid, $id )
	  or die $sth->errstr;
	$result->{killed} = $sth->rows;
	return $result;
}

=head2 Mark a computation as running

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation
 $host    : (optional) the host it is running on

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
	my $host = shift || get_config_key('jobs_hostid') || hostname;

	return $self->_mark_job( $id_comp, STATE_RUNNING, $host );
}

=head2 Mark a computation as failed

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation
 $message : a message
 $host    : (optional) the host it is running on

 Returns:
 Nothing if successful, dies otherwise

=cut

sub mark_computation_failed {
	my $self    = shift;
	my $id_comp = shift;
	my $message = shift;
	my $host    = shift || get_config_key('jobs_hostid') || hostname;

	return $self->_mark_job( $id_comp, STATE_FAILED, $host, $message );
}

=head2 Mark a computation as complete

 Parameters:
 $self    : a self object
 $id_comp : the ID of the computation
 $message : a message
 $host    : (optional) the host it is running on

 Returns:
 Nothing if successful, dies otherwise

=cut

sub mark_computation_complete {
	my $self    = shift;
	my $id_comp = shift;
	my $message = shift;
	my $host    = shift || get_config_key('jobs_hostid') || hostname;

	return $self->_mark_job( $id_comp, STATE_COMPLETE, $host, $message );
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

	return $self->_mark_job( $id_comp, STATE_QUEUED, "", "" );
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

	my $rdb = get_resultdb( $self->{userid} );
	return $rdb->get_versioned_dataitem( $dataset, $source, $version );
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

	my $rdb = get_resultdb( $self->{userid} );

	## make hashes into JSON
	if ( $type =~ m/JSON/i ) {
		$type = "application/json";
		if ( Serialization::Serializable::is_hash($item) ) {
			$item = to_json( $item, { allow_blessed => 1 } );
		}
	}

	if ( $type eq 'application/json' ) {
		$type = 'JSON';
	}
	if ( $type eq 'text/xml' ) {
		$type = 'XML';
	}
	if ( !defined($source) ) {
		$source = 'unknown';
	}

	my $set   = $rdb->get_or_create_dataset($dataset);
	my $id_ds = $set->{'id_ds'};

	my $return = {};

	$return->{"$id_ds"} = $rdb->store_dataitem( $id_ds, $type, $source, $item );

	$return->{upload} = "ok";
	return $return;
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

	my $rdb = get_resultdb( $self->{userid} );

	return $rdb->list_versions( $dataset, $source );
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

	my $rdb = get_resultdb( $self->{userid} );

	if ( scalar @$rversions == 0 ) {
		$rversions = $rdb->list_versions( $dataset, $source );
		## get latest version
		my $item = $rdb->get_versioned_dataitem( $dataset, $source );
		if ( !defined( $item->{id_di} ) ) {
			die "Item not found";
		}
		my @versions = grep { int($_) != $item->{id_di} } @$rversions;
		$rversions = \@versions;
	}

	foreach (@$rversions) {
		$rdb->delete_dataitem($_);
	}

	return { purged => join( ':', @$rversions ), };
}

=head1 Database Operations

=cut

=head2 Get job bucket information

 Parameters:
 $self     : a self object
 $result       : (optional) hashref to append results to

 Returns:
 $result = {

 }
=cut

sub _buckets {
	my $self   = shift;
	my $userid = $self->{userid};

	my $result = shift || {};

	my $dbh = get_dbh;

	my $sth = $dbh->prepare("SELECT * FROM `computations` WHERE `id_user` = ?")
	  or die "Cannot access the computation database.";
	$sth->execute($userid)
	  or die $sth->errstr;

	my $buckethash = {};
	my $count      = 0;
	while ( my $ref = $sth->fetchrow_hashref() ) {
		my $overall = $ref->{procs} * $ref->{threads};

		if ( !defined( $buckethash->{"$ref->{procs}:$ref->{threads}"} ) ) {
			$buckethash->{"$ref->{procs}:$ref->{threads}"} = 1;
		} else {
			$buckethash->{"$ref->{procs}:$ref->{threads}"}++;
		}

		if ( !defined( $buckethash->{"$overall"} ) ) {
			$buckethash->{"$overall"} = 1;
		} else {
			$buckethash->{"$overall"}++;
		}
		++$count;
	}
	$buckethash->{count} = $count;
	$result = $buckethash;

}

=head2 Submit a job

 Parameters:
 $self     : a self object
 $jsontext : the JSON representation of the job
 $result   : (optional) hashref to append results to

 Returns:
 $result = {
 	id_comp => the id of the inserted computation
 	submitted => 1 if successful, 0 otherwise
 }

=cut

sub _submit_job {
	my $self     = shift;
	my $userid   = $self->{userid};
	my $jsontext = shift;
	my $result   = shift || {};

	my $jshash = from_json($jsontext);
	my $canon_json = to_json( $jshash, { allow_blessed => 1, canonical => 1 } );

	my $dbh = get_dbh;

## we don't validate here -- the server running the job service
## doesn't necessarily need to be able to validate all
## computations it stores, the job runner will be able to
## recognize the ones it can actually run
	#		my $comp = Serialization::Serializable::from_hash($jshash);
	#		$comp->validate;

	my $sth = $dbh->prepare(
"SELECT * FROM `computations` WHERE `json_hash` = UNHEX(MD5(?)) AND `id_user` = ?"
	) or die "Cannot access the computation database.";
	$sth->execute( $canon_json, $userid )
	  or die $sth->errstr;

	my $found = undef;
	while ( my $ref = $sth->fetchrow_hashref() ) {
		if ( $ref->{json} eq $canon_json ) {
			$found = $ref;
			last;
		}
	}

	if ( defined($found) ) {
		## reset the computation
		$result->{id_comp}   = $found->{id_comp};
		$result->{'state'}     = $found->{'state'};
		$result->{submitted} = 0;
	} else {
		$sth = $dbh->prepare(
"INSERT INTO `computations` (`json_hash`, `json`, `id_user`) VALUE (UNHEX(MD5(?)),?,?)"
		);
		$sth->execute( $canon_json, $canon_json, $userid )
		  or do {
			debug("Could not create computation : $canon_json");
			confess("Could not create computation.");
		  };

		my $id_comp =
		  $dbh->last_insert_id( undef, undef, qw(computations id_comp) );
		$result->{id_comp}   = $id_comp;
		$result->{submitted} = 1;
	}

	return $result;
}

=head2 Get job information

 Parameters:
 $self     : a self object
 $result       : (optional) hashref to append results to

 Returns:
 $result = [ info, ... ]
=cut

sub _list_jobs {
	my $self   = shift;
	my $userid = $self->{userid};
	my $result = shift || {};

	my $dbh = get_dbh;

	$result->{'jobs'} = [];

	## TODO perform maintenance : reset all computations to waiting which
	## have been running longer than max time

	my $sth = $dbh->prepare("SELECT * FROM `computations` WHERE `id_user` = ?")
	  or die "Cannot access the computation database.";
	$sth->execute($userid)
	  or die $sth->errstr;

	my $found = undef;
	while ( my $ref = $sth->fetchrow_hashref() ) {
		push @{ $result->{jobs} },
		  {
			'id_comp' => $ref->{'id_comp'},
			'json'    => $ref->{'json'},
			'state'   => $ref->{'state'},
			'host'    => $ref->{'host'},
			'procs'   => $ref->{'procs'},
			'threads' => $ref->{'threads'},
			'output'  => $ref->{'output'},
		  };
	}

	return $result;
}

=head2 Check if there are pending computations and submit a runner
 if necessary

 Parameters:
 $self     : a self object

=cut

sub _check_pending_computations {
	my $self   = shift;
	my $userid = $self->{userid};

	if (
		defined (get_config_key ('jobservice_config')->{no_runner}) &&
		get_config_key ('jobservice_config')->{no_runner} == 1
	)
	{
		## no runner -> no check
		return;
	}

	## run the checker
	my $runner_dir = File::Spec->catdir( $basedir, "Runner" );
	my $runner_file = File::Spec->catdir( $basedir, "Runner", "runner.pl" );

	if ( !-d $runner_dir || !-x $runner_file ) {
		die "Runner was not found in $runner_dir.";
	}

	my $jobs = $self->_list_jobs();

	eval {
		my $njobs = 0;
		foreach my $j ( @{ $jobs->{jobs} } ) {
			if ( $j->{'state'} == STATE_QUEUED ) {
				++$njobs;
			}
		}

		if ( $njobs > 0 ) {
			Scheduler::PBS::check_limits;

			eval {

				Scheduler::PBS::qsub_runner( 1, 1, "$basedir/Runner/runner.pl",
					"", { ORANGES_USERID => $userid, } );
			};
			if ($@) {
				warn $@;
			}
		}
	};
	if ($@) {
		unless ( $@ =~ m/exceeded for current user/ ) {
			die($@);
		}
	}
}

=head2 Mark job to have a specified state

 Parameters:
 $self         : a self object
 $id           : the id to remove
 $target_state : the target state
 $host         : the host
 $message      : (optional) a message
 $result       : (optional) hashref to append results to

 Returns:
 $result = {
 	mark_job => <state the job was set to>
 }
=cut

sub _mark_job {
	my $self         = shift;
	my $userid       = $self->{userid};
	my $id           = shift;
	my $target_state = shift;
	my $host         = shift;
	my $message      = shift || "";
	my $result       = shift || {};

	my $dbh = get_dbh;

	if ( !defined($host) ) {
		die "A host must be specified when marking a computation running.";
	}

	my $state_validation = "";
	given ($target_state) {
		when (STATE_RUNNING) {
			$state_validation = " AND `state` != " . int(STATE_RUNNING);
		}
		when (STATE_FAILED) {
			$state_validation = " AND `state` = " . int(STATE_RUNNING);
		}
		when (STATE_COMPLETE) {
			$state_validation = " AND `state` = " . int(STATE_RUNNING);
		}
	}

	my $sth = $dbh->prepare(
		"SELECT * FROM `computations` WHERE `id_user` = ? AND `id_comp` = ? ")
	  or die "Cannot access the computation database.";
	$sth->execute( $userid, $id )
	  or die $sth->errstr;

	if ( my $ref = $sth->fetchrow_hashref() ) {
		$sth = $dbh->prepare(
"UPDATE `computations` SET `state` = ?, `start` = NOW(), `host` = ?, `output` = ? WHERE `id_user` = ? AND `id_comp` = ? $state_validation"
		) or die "Cannot access the computation database.";
		$sth->execute( $target_state, $host, $message, $userid, $id )
		  or die $sth->errstr;
		if ( $sth->rows == 0 ) {
			die
"Computation was not in the correct state to go into state $target_state";
		}
	} else {
		die "Computation not found.";
	}
	$result->{mark_job} = $target_state;
	return $result;
}

=head2 Get a list of user ids for which there are waiting computations

 Parameters:
 $self         : a self object

 Returns:
 an ARRAYREF of userids

=cut

sub _get_users_with_pending_computations {
	my $self         = shift;

	my $dbh = get_dbh;

	my $sth = $dbh->prepare ("SELECT DISTINCT `id_user` FROM `computations` WHERE `state` = 0");

	debug ( "Rows : ".  $sth->execute ( ) );

	my $result = [];

	while (my $ref = $sth->fetchrow_hashref ()) {
		push @$result, $ref->{id_user};
	}

	return $result;
}

1;
