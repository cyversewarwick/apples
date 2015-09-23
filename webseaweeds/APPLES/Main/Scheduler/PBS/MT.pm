#!/usr/bin/perl

package Scheduler::PBS::MT;

=head1 Micro job management interface to PBS (Multithreading only)

This scheduler will submit jobs that run on a single node, using 
multiple threads.

=cut

use strict;

our @ISA = qw(Scheduler::PBS);
require Scheduler::PBS;

use Configuration::AppleSeeds;

use Serialization::Serializable;
use POSIX;

=head2 Function submit_jobs

 Parameters:
 $jobs : a hash specifiying the types of jobs in the database
 	these keys are used:

		o     => n   : we have n jobs requiring o cpus

=cut

sub submit_jobs {
	my $jobs       = shift;
	my $parameters = {};

	return if $jobs == 0;

	my $result = {};

	my @threadcounts = ();

	my $tpn = get_config_key('threads_per_node');

	for ( my $t = 0 ; $t < $tpn ; ++$t ) {
		push @threadcounts, 0;
	}

	while ( my ( $k, $v ) = each( %{$jobs} ) ) {
		## single number only, identifies the number of
		## processors.
		if ( $k =~ m/^[0-9]+$/ ) {
			my $p = int($k);

			## limit to threads per node
			if ( $p > $tpn ) {
				$p = $tpn;
			}

			$threadcounts[$p]++;
		}
	}

	## submit big jobs first.
	for ( my $t = $tpn - 1 ; $t > 0 ; --$t ) {
		for ( my $njobs = 0 ; $njobs < ceil( $threadcounts[$t] ) ; ++$njobs )
		{
			Scheduler::PBS::check_limits();
			Scheduler::PBS::qsub_runner( 1, $t );
		}
	}
}

1;
