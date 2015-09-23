#!/usr/bin/perl

package Scheduler::PBS;

=head1 Micro job management interface to PBS

Contains functions to keep track of a list of PBS jobs which 
run the Perl script "runner.pl" from the 'jobs' subdirectory.

=cut

use strict;

use Compress::Zlib;

use JSON;

use Digest::MD5 qw(md5_hex);
use POSIX qw(strftime ceil floor);
use Data::Dumper;
use Cwd qw(abs_path);
use File::Temp qw(tempfile tempdir);
use FindBin qw($Bin);

use Configuration::ReadConfig qw($basedir);
use Configuration::AppleSeeds;

use Sys::Hostname;

use Runtime;
use Carp;

use Serialization::Serializable;

=head2 Function submit_jobs

 Parameters:
 $jobs : a hash specifiying the types of jobs in the database
		keys can be either
		
		count => x   : we have x jobs overall
				
		or
		
		o     => n   : we have n jobs requiring o MPI cpus

		or 
		
		a:b   => m   : we have m jobs requiring a nodes with b threads

 !! Override this. !!
 
=cut

sub submit_jobs {
	confess("Virtual function called.");
}

=head2 Function get_user_jobs

 Get all jobs for the current user.
 
 Parameters:
 none

 Returns:
 job information as an ARRAYREF like this:
 [
  {
  	jobid => ...,
  	cpus => ...,
  	dir => ...,
  	state => ... U : unknown : Q : queued, R: Running
  } ,
  ...
 ]
=cut

sub get_user_jobs {
	## qstat is used to find the job ids we are interested in
	my $qstat = get_config_key('qstat_command');
	## qinfo is something like qstat -f (torque) or qstat -j (SGE)
	my $qinfo = get_config_key('qinfo_command');

	my @result = ();

	my $j = `$qstat`;
	my @jobs = split /\n/, $j;

	foreach $j (@jobs) {
		## try to extract numerical job id
		if ( $j =~ m/^\s*([0-9]+)/ ) {
			my $jobid = $1;

			my $qinf = `$qinfo $jobid`;

			my $state = 'U';
			## on openpbs, the state is part of the qsub -f output
			if ( $qinf =~ m/job_state\s*[\=\:]\s*([^\n])/ ) {
				$state = $1;
			}

			my $dir = undef;
			## on Torque/openPBS, get the path from the output file
			if ( $qinf =~ m/Output_Path\s*[\=\:]\s*([^\n]+)/ ) {
				$dir = $1;
				## remove output file name, just return directory
				$dir =~ s/[^\/]+\.pbs\.o[0-9]+$//;
				## of course, on SGE, things are different.
			} elsif ( $qinf =~ m/sge_o_workdir\s*[\=\:]\s*([^\n]+)/ ) {
				$dir = $1;
			}

			my $overall_cpus = 1;
			my $cpus         = 1;
			my $threads      = 1;
			if (
				$qinf =~ m/Resource_List.nodes\s*\=\s*([0-9]+)\:ppn\=([0-9]+)/ )
			{
				$overall_cpus = int($1) * int($2);
				$cpus         = int($1);
				$threads      = int($2);
			} elsif ( $qinf =~ m/threads\s*[\=\:]\s*([0-9]+)/ ) {
				$threads      = $qinf;
				$overall_cpus = $threads * $cpus;
			}

			push @result,
			  {
				'jobid'   => $jobid,
				'state'   => $state,
				'dir'     => $dir,
				'cpus'    => $cpus,
				'threads' => $threads,
				'cpus'    => $overall_cpus,
			  };
		}
	}

	return \@result;
}

=head2 Submit a PBS job that requests a fixed number of threads+nodes

 Parameters:
 
 Returns:
 nothing.

=cut

sub qsub_runner {
	my $nodes         = shift;
	my $threads       = shift;
	my $script_to_run = shift || "$basedir/Runner/runner.pl";
	my $script_args   = shift || "";
	my $env           = shift || {};

	if ( defined( $ENV{ORANGES_USERID} ) && !defined( $env->{ORANGES_USERID} ) ) {
		$env->{ORANGES_USERID} = $ENV{ORANGES_USERID};
	}

	my $jobtempdir = get_config_key('jobtempdir');

	my $tempdir = tempdir( DIR => "$jobtempdir/", CLEANUP => 0 );

	if ( !defined($tempdir) ) {
		confess("cannot create temporary directory");
	}

	$tempdir = abs_path($tempdir);

	# copy the job script over
	qx(cp $script_to_run $tempdir/);

	## jobs use a special log4perl config file, which
	## uses stderr (so stuff shows up in the log output)
	qx(cp $basedir/Runner/log4perl.properties $tempdir/);

	qx(cp -R $basedir/Runner/Configuration $tempdir/)
	  if ( -d "$basedir/Runner/Configuration" );

	open COMMANDFILE, ">$tempdir/run_ws_job.sh";
	my $perl = get_perl();

	## reduce the threads automatically to what we have available
	my $tpn = get_config_key('threads_per_node');
	$threads = ( $threads > $tpn ) ? $tpn : $threads;

	## this is the string that allows us to configure
	## processors/nodes/threads in MPI
	my $mpipart = get_config_key('qsub_mpispec');

	## number of processors = nodes * threads
	my $procs = $nodes * $threads;

	$mpipart =~ s/\<\<PROCS\>\>/$procs/;
	$mpipart =~ s/\<\<THREADS\>\>/$threads/;
	$mpipart =~ s/\<\<NODES\>\>/$nodes/;

	my $exec_script = $script_to_run;
	$exec_script =~ s/(.*[\/\\])?([^\/\\]+)$/$2/;

	my $strenv = "";
	while ( my ( $k, $v ) = each(%$env) ) {
		$strenv .= "setenv $k $v\n";
	}

	## make a PBS shell file
	print COMMANDFILE "#!/bin/csh\n\n" 
	  . $mpipart
	  . "touch ./running_flag\n"
	  . "setenv BASEDIR "
	  . $basedir . "\n"
	  . "setenv APPLES_THREADS "
	  . $threads . "\n"
	  . "setenv APPLES_NODES "
	  . $nodes . "\n"
	  . "setenv APPLES_PROCS "
	  . $procs . "\n"
	  . "setenv ADDITIONAL_CONFIG_SUFFIX "
	  . hostname . "\n"
	  . "setenv PERL5LIB "
	  . ( get_config_key("PERL5LIB") || "" ) . "\n"
	  . $strenv
	  . "$perl $exec_script $script_args\n"
	  . "rm -rf $tempdir/running_flag && touch $tempdir/finished_flag\n";
	close COMMANDFILE;

	qx(chmod 755 $tempdir/run_ws_job.sh);

	my $qs_out;

	## run it through the queue
	my $qsub         = get_config_key('qsub_command');
	my $qsub_command = "cd $tempdir && $qsub ./run_ws_job.sh 2>&1";
	$qs_out = `$qsub_command`;

	debug("Job for parameter set up to execute in $tempdir: $qs_out.");

	open PROGRESS, ">$tempdir/progress.txt";
	print PROGRESS $qs_out;
	print PROGRESS "\nJob submitted.\n";
	close PROGRESS;

	## Ideal answer: Your job [0-9]+ ("run_ws_job.sh") has been submitted
	if ( $qs_out =~ m/([0-9]+) .*run\_ws\_job\.sh\"\) has been submitted/ ) {
		debug("Job submitted as $1");
	} else {
		die("Failed to submit job to PBS, could not obtain PBS ID: $qs_out .");
	}
}

=head2 Check if we are still within the job limits

Dies with an exception if the job limits are exceeded.

 Returns:
 ( jobs_left, cpus_left )

=cut

sub check_limits {
	## Check if we can still create new jobs.
	my $currentjobs = get_user_jobs();

	## count running and waiting jobs
	my $total_jobs = 0;
	my $total_cpus = 0;

	foreach my $job (@$currentjobs) {
		++$total_jobs;
		$total_cpus += $job->{overall_cpus};
	}

	if ( $total_jobs >= get_config_key('njobs_limit') ) {
		confess(
"Maximum number of jobs was exceeded for current user. $total_jobs / "
			  . get_config_key('njobs_limit') );
	}
	if ( $total_cpus >= get_config_key('ncpus_limit') ) {
		confess(
"Maximum number of cpus was exceeded for current user. $total_cpus / "
			  . get_config_key('ncpus_limit') );
	}

	$total_jobs = get_config_key('njobs_limit') - $total_jobs;
	$total_cpus = get_config_key('ncpus_limit') - $total_cpus;

	return ( $total_jobs, $total_cpus );
}

1;
