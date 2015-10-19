
=head1 AppleSeeds Package

This package provides basic logging and configuration functionality for
the data interface.

Also, the package wraps creation of CGI sessions and ResultDB access
for webtool use.

=cut

package Configuration::AppleSeeds;

use strict;
no warnings;    # in EPIC, somehow the include paths get messed up.
                # this prevents the warnings about redefinitions

our @ISA = qw(Exporter);
require Exporter;
our @EXPORT    = qw(logger get_config_key find_executable);
our @EXPORT_OK = qw(load_APPLES add_fasta_dbs fasta_prot_or_nuc);

our $perl_command = $^X;
$perl_command = qq{"$perl_command"}
  if $perl_command =~ /\s/;

use Carp;
use Data::Dumper;
use Scalar::Util qw(blessed);

use Configuration::ReadConfig qw(:ALL);
use Log::Log4perl qw(get_logger);
use File::Spec qw (catdir);

=head2 Add a set of FASTA databases 

 Parameters:
 @_ : a list of FASTA file names to be added both as
      BLAST and Sequence databases

=cut

sub add_fasta_dbs {
	my @a = @_;
	our %tools;
	our %BLAST_parameters;
	our %sequencedatabases;
	our %searchdatabases;

	foreach my $db (@a) {
		my $pn = fasta_prot_or_nuc($db);

		$db =~ s/^$config{fastapath}//;
		$db =~ s/^[\/\\]//;

		# Runtime::debug("Adding $db ($pn). \n");

		foreach my $tool ( @{ $tools{$pn} } ) {
			my $searchdb =
			  Sequences::Database::Search::BLAST_Fasta->new( $db, $tool, );

			my $blast_parameters = $BLAST_parameters{$tool}
			  || $BLAST_parameters{''};

			$searchdatabases{ $db . "_" . $tool } = {
				name        => $db,
				description => "Local BLAST on $db",
				db          => $searchdb,
				key         => $db . "_" . $tool,
				parameters  => $blast_parameters,
			};
			$sequencedatabases{$db . "_" . $tool} = {
				name        => $db,
				description => "Local file $db" . "_" . $tool,
				db          => $searchdb,
				key         => $db,
			};
		}
	}
}


=head2 Load the APPLES Library

For scripts that use APPLES, we can use this subroutine to load the APPLES
library path.

This will add the library path using use lib, so this is the way it should
be used:

 BEGIN {
 	use Configuration::AppleSeeds qw(load_APPLES);

 	load_APPLES();
 }

=cut

sub load_APPLES {
	if ( !defined($confpath) ) {
		$confpath = "";
	}
	our $apples_dir = $config{APPLES_directory}
	  || File::Spec->catdir( $basedir, "APPLES", "Main" );
	# print "load_APPLES 105, config\n"  . Dumper(\%config);
	$apples_dir = File::Spec->rel2abs($apples_dir);

	if ( $apples_dir
		&& -d $apples_dir )
	{
		push @INC, $apples_dir;
		eval "use Runtime;";
		if ($@) {
			confess(
				"APPLES modules not found in '$apples_dir' / $confpath : $@."
				  . Dumper( \%config ) );
		}
	} else {
		eval "use Carp qw(confess);";
		confess("APPLES installation could not be found in '$apples_dir'.");
	}

	# load ensembl
	my $ensembl_loaded = 0;
	if ( defined( $config{ensembl} ) && $config{ensembl} ) {
		my @libs = (
			"ensembl",           "ensembl-compara",
			"ensembl-variation", "ensembl-functgenomics",
		);
		foreach (@libs) {
			my $ensemblpath =
			  File::Spec->catdir( $config{ensembl}, $_, "modules" );
			if ( $ensemblpath && -d $ensemblpath ) {
				push @INC, $ensemblpath;
				++$ensembl_loaded;
			}
		}
	}

	if ($ensembl_loaded < 1) {
		Runtime::warn("Could not load EnsEMBL libraries. Config is: " . Dumper (\%config));
	}

	# Load sequence databases
	our %searchdatabases   = ();
	our %sequencedatabases = ();

	our %tools = (
		n => [ 'blastn', 'tblastn', 'tblastx' ],
		p => [ 'blastp', 'blastx' ],
	);

	our %BLAST_parameters = (
		''        => { '-num_threads' => '4', },
		'tblastx' => {
			'-word_size'   => 2,
			-'evalue'      => 1,
			'-seg'         => 'no',
			'-num_threads' => '4',
		},
		'tblastn' => {
			'-word_size'   => 2,
			-'evalue'      => 1,
			'-seg'         => 'no',
			'-num_threads' => '4',
		}
	);

	if ($config{blastdb}) {
		eval("require Sequences::Database::Search::BLAST_Fasta");
		die $@ if $@;
		eval("require Sequences::Database::Search::BLAST_Genbank");
		die $@ if $@;

		if (ref($config{blastdb}) ne 'ARRAY') {
			$config{blastdb} = [ $config{blastdb} ];
		}

		foreach my $gbc (@{$config{blastdb}}) {
			my $searchdb;
			debug ("adding DB " . Dumper($gbc));
			my $db = $gbc->{db} || 'nt';
			my $tool = $gbc->{tool} || 'blastn';
			if ($gbc->{type} =~ m/fasta/i) {
				$searchdb =
				  Sequences::Database::Search::BLAST_Fasta->new( 
				  	-DB => $db, 
				  	-PROGRAM => $tool, );
			} elsif ($gbc->{type} =~ m/genbank/i) {
				$searchdb =
				  Sequences::Database::Search::BLAST_Genbank->new( 
				  	-DB => $db, 
				  	-PROGRAM => $tool, 
				  	-ENTREZQUERY => $gbc->{entrez_query} || ''
				  	);
 			} else {
				die "Invalid BLAST DB declaration " . Dumper($gbc);
			}

			my $blast_parameters = $gbc->{BLAST_parameters}
			  || $BLAST_parameters{$tool}
			  || $BLAST_parameters{''};

			my $key = $gbc->{name} || $db;
			$searchdatabases{ $key . "_" . $tool } = {
				name        => $key,
				description => $gbc->{description} || "$db ($tool)",
				db          => $searchdb,
				key         => $db . "_" . $tool,
				parameters  => $blast_parameters,
			};
		}
	}

	unless ( defined ( $main::NO_APPLES_DBS ) && $main::NO_APPLES_DBS ) {
		eval("require Sequences::Database::Search::BLAST_Fasta");
		die $@ if $@;

		if ( defined( $config{fastapath} ) && -d $config{fastapath} ) {
            my $fn = File::Spec->catfile ($config{fastapath}, '*.fa');
            my @a = glob $fn;
            my $fn2 = File::Spec->catfile ($config{fastapath}, 'nr');
			if (-e $fn2) {
				push @a, $fn2;
			}

			add_fasta_dbs(@a);
		}

		eval "require Sequences::Database::Sequence::Ensembl;";
		if ($@) {
			die $@;
		}
		
		eval "require Sequences::Database::Sequence::Genbank;";
		if ($@) {
			die $@;
		}

		$sequencedatabases{'ensembl'} = {
			name        => 'ensembl',
			description => "Ensembl Online Database",
			db          => Sequences::Database::Sequence::Ensembl->new,
			key         => 'ensembl',
		};

		$sequencedatabases{'ensemblgenomes'} = {
			name        => 'ensembl',
			description => "EnsemblGenomes Online Database",
			db  => Sequences::Database::Sequence::Ensembl->new('ensemblgenomes'),
			key => 'ensemblgenomes',
		};

    	$sequencedatabases{'ensembl_local'} = {
			name        => 'ensembl_local',
			description => "Ensembl Local Databases",
			db  => Sequences::Database::Sequence::Ensembl->new('local'),
			key => 'ensembl_local',
		};

		# $sequencedatabases{'genbank'} = {
		# 	name        => 'genbank',
		# 	description => "Genbank Nucleotide Online Database",
		# 	db          => Sequences::Database::Sequence::Genbank->new,
		# 	key         => 'genbank',
		# };
	}

	$Runtime::sequencedatabases = \%sequencedatabases;
	$Runtime::searchdatabases   = \%searchdatabases;

	unless ( defined( $config{usertempdir} )
		&& -d $config{usertempdir}
		&& -w $config{usertempdir} )
	{
		confess(
"Cannot write user-specific tempdir, please specify in config file... ( "
			  . "$config{usertempdir} )" );
	}

	## set $Runtime::jobservice according to the job services config
	eval {
		## use the first one
		if ( ref( $config{jobservice_config} ) eq 'ARRAY' ) {
			$config{jobservice_config} = $config{jobservice_config}->[0];
		}
	};
	if (  !defined( $config{jobservice_config} )
		|| ref( $config{jobservice_config} ) ne 'HASH' )
	{
		confess(
"Job service configuration was not found in config files : @{$config{parsed_config_files}} "
			  . Dumper( \%config ) );
	}

	my %jobservice_config = %{ $config{jobservice_config} };

	$Runtime::jobservice = undef;
	eval "require Runtime::Job_Service";
	if ($@) {
		die "Runtime::Job_Service not found, config is: " . Dumper (\%config);
	}

	eval {
		if ( $jobservice_config{type} eq 'remote' )
		{
			eval "require Runtime::Remote_Job_Service";
			if ($@) {
				die "Runtime::Remote_Job_Service not found, config is: " . Dumper (\%config);
			}
			$Runtime::jobservice = Runtime::Remote_Job_Service->new(
				$jobservice_config{url},
				$jobservice_config{user},
				$jobservice_config{password},
				$jobservice_config{http_user},
				$jobservice_config{http_password},
			);
			## see if it authenticates without error
			## This will be done automatically when the service is used,
			## we don't do it here to make the scripts start up faster.
			# $Runtime::jobservice->_authenticate;
		} elsif ( $jobservice_config{type} eq 'web' ) {
			eval "require Runtime::Web_Job_Service";
			if ($@) {
				die "Runtime::Web_Job_Service not found, config is: " . Dumper (\%config);
			}

			my $userid = $jobservice_config{userid}
			  || $ENV{ORANGES_USERID};
			$Runtime::jobservice = Runtime::Web_Job_Service->new($userid);

		} elsif ( $jobservice_config{type} eq 'local' ) {
			$Runtime::jobservice =
			  Runtime::Job_Service->new( $jobservice_config{datadir} );
		}
	};
	if ($@) {
		warn("Could not create job service, using local mode. $@");
	}

	if ( !defined($Runtime::jobservice) && $jobservice_config{type} ne 'none' )
	{
		warn("No job service was set up.");
	}
}

=head2 Check if a fasta file is protein or nucleotide database

 Parameters:
 $name : file name

 Returns : 'n' or 'p'

=cut

sub fasta_prot_or_nuc {
	my $name = shift;

	open FASTAFILE, "<", $name
	  or die "File $name not found";

	my $samples = 10;    ## we sample 10 lines.
	my $ret     = 'n';
	while (<FASTAFILE>) {
		my $r = $_;
		unless ( $r =~ m/^\>/ ) {    # skip comments
			                         # these characters only exist in proteins
			if ( $r =~ m/[EFILOPQZ*]/ ) {
				$ret = 'p';
			}
		}
		--$samples;
		last if $samples <= 0 || $ret eq 'p';
	}

	close FASTAFILE;

	return $ret;
}

=head2 Function logger
Returns a log4perl logger
=cut

sub logger {
	return get_logger;
}

=head2 Get a key from the ORANGES configuration.

 Parameters:
 $key : key string for the configuration

 Returns:
 the value of this config key (or undef)

=cut

sub get_config_key {
	my $key = shift;
	return $config{$key};
}

=head2 Platform-independent function to get executable names.

Uses File::Spec to construct a file name, and appends ".exe" when appropriate.

The paths to executables can be set in the configuration. You can either
specify individual keys like

 path_to_executable = /usr/local/bin/executable

or add extra directories to the path, by setting the config setting

 PATH=/opt/mybin:/usr/local/mypackage

in which we will then look for the executables specified.

=cut

sub find_executable {
    my $p = File::Spec->catfile (@_);
    my $found = undef;

    if (get_config_key ("path_to_$p")) {
        $found = get_config_key ("path_to_$p");
    } else {
        if ($^O eq 'MSWin32') {
            $p.= ".exe";
        }

        my $extrapath = get_config_key("PATH") || "";
        my $path = $extrapath . ":" . $ENV{PATH};

        my @paths = split /\:/, $path;

        foreach my $pp (@paths) {
            my $try = File::Spec->catfile($pp, $p);
            if (-e $try && -x $try) {
                $found = $try;
                last;
            }
        }
    }

    if (!defined ($found) || !-x $found) {
        die "Could not find executable $p / $found";
    }

    return $found;
}

1;
