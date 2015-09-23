#!/usr/bin/perl

=head1 Configuration Helper Script

Will read the configuration and fill in config variables

=cut

package Configuration::ReadConfig;

use strict;
use warnings;

=head2 Export stuff

=cut

require Exporter;

our @EXPORT_OK = qw($basedir $database $hostname $port $user $password
  %config $confpath %webforms
  %webform_links);
our %EXPORT_TAGS = (
	ALL => [
		qw($basedir $database $hostname $port $user $password
		  $confpath %webforms
		  %config)
	]
);
our @ISA = qw(Exporter);

use FindBin qw($Bin $Script);
use File::Spec qw(catdir catfile);
use Sys::Hostname;
use Config::General;
use Data::Dumper;

=head2 Initialisation

Find the configuration and fill in default values

=cut

sub init {
	our $confpath = "";

	foreach my $p (@INC) {
		$confpath = File::Spec->catdir( $p, "Configuration" );
		last if -d $confpath;
	}

	unless ( defined($confpath) && -d $confpath ) {
		die "Configuration files not found.";
	}

	# if config reading doesn't work, here's a debug output that might help
	# print "Confpath: $confpath (@INC) \n";

	## This is the configuration read from the host-specific file
	our %config;

	## read config from file
	my $host             = hostname;
	my $conf_read        = 0;
	my @config_locations = (
		File::Spec->catfile( $confpath, "config.txt" ),
		File::Spec->catfile( $confpath, "config_$host.txt" ),
		File::Spec->catfile( ( $Bin, "Configuration" ), "config.txt" ),
		File::Spec->catfile( ( $Bin, "Configuration" ), "config_$host.txt" ),
		File::Spec->catfile( $Bin, "config.txt" ),
		File::Spec->catfile( $Bin, "config_$host.txt" ),
	);

	if (defined ($ENV{ADDITIONAL_CONFIG_PATHS})) {
		my @paths = split/\:/, $ENV{ADDITIONAL_CONFIG_PATHS};
		foreach my $p (@paths) {
			push @config_locations, File::Spec->catfile( $p, "config.txt" );
			push @config_locations, File::Spec->catfile( $p, "config_$host.txt" );
		}
	}

	if (   defined( $ENV{ADDITIONAL_CONFIG_SUFFIX} )
		&& $ENV{ADDITIONAL_CONFIG_SUFFIX} ne ""
		&& $ENV{ADDITIONAL_CONFIG_SUFFIX} ne $host )
	{
		my $suffix = $ENV{ADDITIONAL_CONFIG_SUFFIX};
		push @config_locations,
		  File::Spec->catfile( $confpath, "config_$suffix.txt" );
		push @config_locations,
		  File::Spec->catfile( ( $Bin, "Configuration" ),
			"config_$suffix.txt" );
		push @config_locations,
		  File::Spec->catfile( $Bin, "config_$suffix.txt" );
	}

	$config{parsed_config_files} = [];
	foreach my $config_file_name (@config_locations) {
		if ( -e $config_file_name ) {
			push @{ $config{parsed_config_files} }, $config_file_name;
			my $conf = Config::General->new(
				-ConfigFile => $config_file_name,

			 # TODO wsbc doesn't quite have this yet. we can probably do without
			 #				-ForceArray => 1,
			);

			my %readconf = $conf->getall();
			foreach my $k ( keys %readconf ) {
				if ( $k eq 'PERL5LIB' ) {
					$config{PERL5LIB} .= ":" . $readconf{$k};
				} elsif ( $k eq 'webforms' ) {
					foreach my $kk (%{$readconf{$k}}) {
						$config{$k}->{$kk} = $readconf{$k}->{$kk};
					}
				} else {
					$config{$k} = $readconf{$k};
				}

				$conf_read++;
			}

		}

	}
	if ( !$conf_read ) {
		die "You must provide a host specific configuration file at "
		  . "one of these locations: @config_locations \n";
	}

	## re-execute using correct version of perl.
	if ( defined( $config{perl} ) && $config{perl} ne $^X ) {
		my $pscript = File::Spec->catfile( $Bin, $Script );
		my @args = @ARGV;
		unshift @args, $pscript;
		if ( defined( $config{PERL5LIB} ) ) {
			$ENV{PERL5LIB} = $config{PERL5LIB};
		}
		exec( $config{perl}, @args );
		exit;
	}

=pod install additional include paths (like a local CPAN repository)

=cut

	if ( defined( $config{PERL5LIB} ) ) {
		my @libs = split /\:/, $config{PERL5LIB};
		push @INC, $_ foreach @libs;
	}

=pod set up web forms
=cut

	our %webforms      = ();
	our %webform_links = ();

	if ( ref( $config{webforms} ) eq 'HASH' ) {
		%webforms = %{ $config{webforms} };

		## generate link URLs for the forms service
		foreach my $f ( keys %{ $webforms{form} } ) {
			$webform_links{ $webforms{form}->{$f}->{class} } =
			  "form.pl?formid=$f";
		}
	}

=pod parameters for database access
 TODO: perhaps add table creation using sql script
=cut

	if ( defined( $config{resultdb} ) ) {
		if ( ref( $config{resultdb} ) ne 'HASH' ) {
			die "You need to specify resultDB access : " . Dumper( \%config );
		} else {
			our $database = $config{resultdb}->{database}
			  || 'apples_webpage_users';
			our $hostname = $config{resultdb}->{hostname} || $host;
			our $port     = $config{resultdb}->{port}     || '3306';
			our $user     = $config{resultdb}->{user}     || 'root';
			our $password = $config{resultdb}->{password} || '';
		}
	}

=pod parameters for data storage and locations.
=cut

	our $basedir =
	     $config{basedir}
	  || $ENV{'BASEDIR'}
	  || "$Bin";

=pod defaults for job scheduler
=cut

	$config{jobtempdir} = $config{jobtempdir}
	  || '/Users/Shared/www/webservices/webseaweeds';
	$config{njobs_limit} = $config{njobs_limit}
	  || 1;
	$config{ncpus_limit} = $config{ncpus_limit}
	  || 1;
	$config{job_retry_count} = $config{job_retry_count}
	  || 3;
	$config{threads_per_node} = $config{threads_per_node}
	  || 4;
	$config{qsub_mpispec} = $config{qsub_mpispec}
	  || "#\$ -pe lam 1-<<PROCS>>\n"
	  . "#PBS -l nodes=<<NODES>>:ppn=<<THREADS>>,pvmem=1gb\n";

	## PBS interface
	$config{qstat_command} = $config{qstat_command}
	  || File::Spec->catfile( ( $basedir, "Runner" ), "dummyqstat.pl" );
	$config{qinfo_command} = $config{qinfo_command}
	  || File::Spec->catfile( ( $basedir, "Runner" ), "dummyqinfo.pl" );
	$config{qsub_command} = $config{qsub_command}
	  || File::Spec->catfile( ( $basedir, "Runner" ), "dummyqsub.pl csh" );
	$config{qdel_command} = $config{qdel_command}
	  || File::Spec->catfile( ( $basedir, "Runner" ), "dummyqdel.pl" );


=pod defaults for runner
=cut

	$config{runner_mem_limit} = $config{runner_mem_limit}
		|| 5 * 1024*1024*1024;

=pod Read logging config and set up logging
=cut

	my $log_conf;
	my @log_paths = ( $Bin, $basedir );

	foreach my $cp (@log_paths) {
		$log_conf = File::Spec->catfile( $cp, "log4perl.properties" );
		last if ( -e $log_conf );
		$log_conf =
		  File::Spec->catfile( $cp, "Configuration", "log4perl.properties" );
		last if ( -e $log_conf );
	}

#	print "Logging config: $log_conf\n";

	{
		local $/ = undef;
		local *FILE;

		# We read this from the current directory, so job scripts can
		# have a different logging configuration
		open FILE, "<", $log_conf
		  or do {
			goto DEFAULTS;
		  };
		$log_conf = <FILE>;
		close FILE;
		goto DONE;

		# By default, only output to stderr
	  DEFAULTS:
		$log_conf = <<END;
log4perl.category = DEBUG, Screen

log4perl.appender.Screen = Log::Log4perl::Appender::Screen
log4perl.appender.Screen.layout = PatternLayout
log4perl.appender.Screen.layout.ConversionPattern = %d> [%p] %C %F:%L %m%n
log4perl.appender.ScreenApp.Threshold = DEBUG
END
	  DONE:
	}

	qx(touch $basedir/logs/webseaweeds.log 2>&1);

	eval "use Log::Log4perl qw(get_logger);";
	if ($@) {
		die "Log4perl not found: $@, Config is: " . Dumper( \%config );
	}

	Log::Log4perl::init( \$log_conf );

	## Set up logging for the case when stuff went seriously wrong
	$SIG{__DIE__} = sub {
		return unless defined $^S and $^S == 0; # Ignore errors in eval

		local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
		my $logger = get_logger("");
		$logger->fatal( Dumper( \@_ ) );
		die @_;    # Now terminate really
	};

	if ( defined( $config{print_config} ) && $config{print_config} ) {
		my $logger = get_logger("");
		$logger->debug( "Configuration read: these are the settings: "
			  . Dumper( \%config ) );
	}

	## marker that configuration has been performed
	our $config_done = 1;
}

=head2 Find the Configuration directory

=cut

BEGIN {
	init;
}

1;
