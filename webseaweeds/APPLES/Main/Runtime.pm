#!/usr/bin/perl

=head1 APPLES Runtime Class

=cut

package Runtime;

use strict;
use feature qw(state);
use feature ":5.10";

use Fcntl;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(warn info debug error user_info config cache get_perl
  get_sequence_database_ids get_sequence_database get_sequence_database_info
  get_search_database_ids get_search_database get_search_database_info
  jobservice APPLES_FULLVERSION APPLES_SERIALIZATIONVERSION APPLES_BRANCH 
  trim uniq read_file write_file eval_job eval_R
  run_perl dump_json
);

use constant {
	APPLES_BRANCH => "grannysmith",
};

use constant {
			   UI_ERROR   => 0,
			   UI_WARNING => 1,
			   UI_INFO    => 2,
			   UI_DEBUG   => 3,
};

use Carp;
use Config::General;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use File::Temp qw (tempdir tempfile);
use Storable;
use MooseX::Declare;
use JSON;
use Cwd qw(getcwd);
use File::Spec qw(catfile);

use Configuration::ReadConfig qw(%config $basedir);

use Log::Log4perl qw(get_logger);

require Runtime::Cache;

=head2 Exception that is "thrown" when Computations were
queued rather than executed

=cut

class Runtime::Job_Queued_Exception {
	has 'job' => (
				   is  => "rw",
				   isa => "Runtime::Computation",
	);
};

## The job submission service
##
our $jobservice = undef;

=head2 Wrapper for old user_info to use Log4Perl

Parameters:
 $level : the logging level (see above)
 $text : a message

Returns:
 Nothing.

=cut

sub user_info {
	my $level = shift;
	my $text  = shift;
	state $logger = get_logger;

	local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;

	given ($level) {
		when (UI_ERROR) {
			$logger->error($text);
		}
		when (UI_WARNING) {
			$logger->warn($text);
		}
		when (UI_INFO) {
			$logger->info($text);
		}
		when (UI_DEBUG) {
			$logger->debug($text);
		}
		default {
			die "Undefined debug level: $text";
		}
	}
}    # user_info #

=head2 Wrapper for outputting info text

Parameters:
 $text : a message

Returns:
 Nothing.

=cut

sub info {
	my $text = shift;

	#local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
	#user_info( UI_INFO, $text );
}

=head2 Wrapper for outputting error text

Parameters:
 $text : a message

Returns:
 Nothing.

=cut

sub error {
	my $text = shift;
	local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
	user_info( UI_ERROR, $text );
}

=head2 Wrapper for outputting info text

Parameters:
 $text : a message

Returns:
 Nothing.

=cut

sub warn {
	my $text = shift;
	#local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
	#user_info( UI_WARNING, $text );
}

=head2 Wrapper for outputting info text

Parameters:
 $text : a message

Returns:
 Nothing.

=cut

sub debug {
	my $text = shift;
	#local $Log::Log4perl::caller_depth = $Log::Log4perl::caller_depth + 1;
	#user_info( UI_DEBUG, $text );
}

=head2 Get the result cache

 Returns:
 the cache object associated with this instance

=cut

sub cache {
	state $cache = Runtime::Cache->new;

	if ( !defined($cache) ) {
		confess("Runtime object used before initialization.");
	}

	return $cache;
}

=head2 Get the APPLES Configuration

 Returns:
 a HASHREF
=cut

sub config {
	return \%config;
}

=head2 get the path to the current perl executable

=cut

sub get_perl {
	my $perl_command = $^X;
	$perl_command = qq{"$perl_command"}
	  if $perl_command =~ /\s/;
	return $perl_command;
}

=head2 Get a sequence database

 Parameters:
 $dbname : the name of the database

 Returns:
 A Sequences::Database::Sequence reference

=cut

sub get_sequence_database {
	## this gets initialized in AppleSeeds::load_APPLES
	our $sequencedatabases;
	my $dbname = _get_matching_key( $sequencedatabases, @_ ) || "@_";

	if ( defined($dbname) && exists( $sequencedatabases->{$dbname} ) ) {
		return $sequencedatabases->{$dbname}->{db};
	}
	confess("Unknown sequence database : $dbname");
}

=head2 Get information about a sequence database

 Parameters:
 $dbname : the name of the database

 Returns:
 { name => <verbose name>, description => <more verbose stuff>}

=cut

sub get_sequence_database_info {
	## this gets initialized in AppleSeeds::load_APPLES
	our $sequencedatabases;
	my $dbname = _get_matching_key( $sequencedatabases, @_ )  || "@_";

	if ( defined($dbname) && exists( $sequencedatabases->{$dbname} ) ) {
		return {
			name => $sequencedatabases->{$dbname}->{name},
			key  => $sequencedatabases->{$dbname}->{key},
			description => $sequencedatabases->{$dbname}->{description},
		};
	}
	confess("Unknown sequence database : $dbname");
}

=head2 Get matching keys for a prefix

 Parameters:
 $hash : a HASHREF
 $string, ... : a list of strings to be contained in the name

 The strings will be split by spaces (split /\s/)

 Returns:
 an name of a key in the hash that matches the prefix

 undef if no such key exists

 dies if key name is ambiguous

=cut

sub _get_matching_key {
	my $hash = shift;

	my @arr = keys %$hash;
	my @px  = ();

	my $wholename = "";
	my $prefix = shift;
	while ( defined($prefix) ) {
		$wholename.= $prefix;
		my @pxs = split /\s/, $prefix;

		push @px, $_ foreach @pxs;
		$prefix = shift;
	}

	foreach $prefix (@px) {
		@arr = grep { $_ =~ m/$prefix/i } @arr;
	}

	if ( scalar @arr > 1 ) {
		foreach my $v (@arr) {
			if ($v eq $wholename) {
				## see if we have an exact match
				return $wholename;
			}
		}

		confess("Name is not unique: @px (matches : @arr)");
	} elsif ( scalar @arr == 1 ) {
		return $arr[0];
	}

	return undef;
}

=head2 Get sequence database id keys to use with the two subs shown above

 Returns:
 ARRAYREF[String]

=cut

sub get_sequence_database_ids {

	## this gets initialized in AppleSeeds::load_APPLES
	our $sequencedatabases;

	my @k = keys %{$sequencedatabases};
	return \@k;
}

=head2 Get a search database

 Parameters:
 $dbname : the name of the database

 Returns:
 A searchs::Database::search reference

=cut

sub get_search_database {
	## this gets initialized in AppleSeeds::load_APPLES
	our $searchdatabases;
	my $dbname = _get_matching_key( $searchdatabases, @_ ) || "@_";

	if ( defined($dbname) && exists( $searchdatabases->{$dbname} ) ) {
		return $searchdatabases->{$dbname}->{db};
	}

	my @sdbs = keys %$searchdatabases;

	confess("Unknown search database : $dbname (I have @sdbs)");
}

=head2 Get information about a search database

 Parameters:
 $dbname : the name of the database

 Returns:
 { name => <verbose name>, description => <more verbose stuff>}

=cut

sub get_search_database_info {
	## this gets initialized in AppleSeeds::load_APPLES
	our $searchdatabases;
	my $dbname = _get_matching_key( $searchdatabases, @_ ) || "@_";

	if ( defined($dbname) && exists( $searchdatabases->{$dbname} ) ) {
		return {
				 name        => $searchdatabases->{$dbname}->{name},
				 key         => $searchdatabases->{$dbname}->{key},
				 description => $searchdatabases->{$dbname}->{description},
				 parameters  => $searchdatabases->{$dbname}->{parameters} || {},
		};
	}
	confess("Unknown search database : $dbname");
}

=head2 Get search database id keys to use with the two subs shown above

 Returns:
 ARRAYREF[String]

=cut

sub get_search_database_ids {

	## this gets initialized in AppleSeeds::load_APPLES
	our $searchdatabases;

	my @k = keys %{$searchdatabases};
	return \@k;
}

=head2 Get the job service

 Returns:
 a Runtime::Job_Service or derived class object to submit
 jobs to and store results with

=cut

sub jobservice {
	our $jobservice;
	return $Runtime::jobservice;
}

=head2 Perl helper: Trim a string by removing leading and trailing whitespace

 Parameters:
 $v : a string

 Returns:
 trimmed version of $v

=cut

sub trim {
	my $v = shift;
	$v=~ s/^\s+//;
	$v=~ s/\s+$//;
	return $v;
}

=head2 Perl helper: Remove duplicates from an array

 Parameters:
 @_ : an array

 Returns:
 Copy of @_ without duplicates

=cut

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

=head2 Read a file into a buffer

See http://www.perl.com/pub/2003/11/21/slurp.html

 Parameters:
 $file_name : the filename

Returns:
 buffer containing the file data in scalar context
 array of lines in array context

=cut

sub read_file {
    my( $file_name ) = @_ ;

    my $buf = "";
    my $buf_ref = \$buf ;

    local( *FH ) ;
    open( FH, "<", $file_name ) or
        die "Can't open $file_name: $!" ;

	binmode FH, ":raw";

    my $size_left = -s FH ;

    while( $size_left > 0 ) {

        my $read_cnt = read( FH, ${$buf_ref},
            $size_left, length ${$buf_ref} ) ;

        unless( $read_cnt ) {

            die "read error in file $file_name: $!" ;
            last ;
        }

            $size_left -= $read_cnt ;
    }
	close FH;

# handle void context (return scalar by buffer reference)
    return unless defined wantarray ;

# handle list context
    return split m|$/|, ${$buf_ref} if wantarray ;

# handle scalar context
    return ${$buf_ref} ;
}

=head2 Write buffer to a file

 Parameters:
 $file_name : the file name
 [optional] $args : hashref arguments
 the rest of @_ : stuff to write to the file

=cut

sub write_file {

    my $file_name = shift ;

    my $args = ( ref $_[0] eq 'HASH' ) ? shift : {} ;
    my $buf = join '', @_ ;

    my $mode = O_WRONLY | O_CREAT;
    $mode |= O_BINARY if $args->{'binmode'} ;
    $mode |= O_APPEND if $args->{'append'} ;

    local( *FH ) ;
    sysopen( FH, $file_name, $mode ) or
        die "Can't open $file_name: $!" ;

    my $size_left = length( $buf ) ;
    my $offset = 0 ;

    while( $size_left > 0 ) {

        my $write_cnt = syswrite( FH, $buf,
                $size_left, $offset ) ;

        unless( $write_cnt ) {

            die "write error in file $file_name: $!" ;
            last ;
        }

        $size_left -= $write_cnt ;
        $offset += $write_cnt ;
    }
	close FH;

    return ;
}

=head2 Run Perl for a given script

This takes care of the perl library path, and runs the script in the perl interpreter we 
are currently using.

 Parameters:
 $filename : name of the file to execute

=cut

sub run_perl {
	my $filename = shift;
	
	use IPC::Open3;
	use Symbol qw(gensym);
	use IO::File;
	use FindBin qw($Bin);
	my $libpath =  $config{"basedir"} . ":" . $Bin;
	if (defined ($config{PERL5LIB})) {
		$libpath.= ":" . $config{PERL5LIB};
	}
	if (defined ($ENV{PERL5LIB})) {
		$libpath.= ":" . $ENV{PERL5LIB};
	}
	local *CATCHERR = IO::File->new_tmpfile;
	my $pid = open3(gensym, \*CATCHOUT, ">&CATCHERR", "PERL5LIB=$libpath $^X $filename");
	while( <CATCHOUT> ) {
		print $_;
	}
	waitpid($pid, 0);
	my $err = $?;
	seek CATCHERR, 0, 0;
	while( <CATCHERR> ) {
		print STDERR $_;
	}
	if ($err) {
		die "Child process exited with $err (script: $filename)"
	}
}

=head2 Eval by running a separate perl instance.

This is a bit of a workaround for the problem that Perl
doesn't give memory back to the system. Also, it can
(but should not) be used to work around leaky code.

 Parameters:
 $job : a string, or a Runtime::Computation
 $parallel : (default : 0) run in parallel

 Returns:
 undef

=cut

sub eval_job {
	my $job = shift;
	my $parallel = shift || 0;

	my ($fh, $filename) = tempfile( DIR => getcwd );
	my $perl =
	print $fh <<END;
#!$^X

use strict;

# This is the loader for APPLES. Go fruit.
BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
	1;
}

use Runtime;

use POSIX;
use Data::Dumper;
use JSON;
use Serialization::Serializable;

END

	if (UNIVERSAL::isa ($job, "Runtime::Computation")) {
		my $json = to_json (Serialization::Serializable::to_hash ($job), {allow_blessed => 1});
		print $fh 'my $jobvar = <<' . "'END';\n";
		print $fh "$json";
		print $fh "\nEND\n";
		print $fh <<END;
my \$job = Serialization::Serializable::from_hash (from_json ( \$jobvar ));
\$job->run();
END
	} else {
		print $fh $job;
	}
	close $fh;
	
	run_perl($filename);

	unlink $filename;
}

our $R = undef;

=head2 Run R code, only allocating a single Statistics::R instance

...since apparently, creating one every time kinda fills up memory

 Parameters:
 $rcode : the R code as a string
 $vars  : a hashref containing variables to assign in R

=cut

sub eval_R {
	my $rcode = shift;
	my $vars = shift || {};

	eval "require Statistics::R;";
	if ($@) {
		die "Cannot run $R script : " . $@;
	}
	our $R;
	if (!defined $R) {
		$R = Statistics::R->new();
	} else {
		$R->run ("rm(list = ls(all = TRUE))");
	}
	$R->run ($rcode);
}

=head2 Return full APPLES version

 Returns:
 Full version which is read from the file version.txt in the 
 APPLES base directory (this is updated via git filters).

 If the file cannot be found, we return

 	unknown.s2000

 (the first serialization version that was committed was 2000)

=cut
our $APPLES_FV = '';
sub APPLES_FULLVERSION {
	our $APPLES_FV;
	return $APPLES_FV if $APPLES_FV ne '';
	
	our $basedir;
	my $versionfile = File::Spec->catfile($basedir, 'version.txt');
	
	my $version;
	{ local $/ = undef; local *FILE; open FILE, "<$versionfile"; $version = <FILE>; close FILE }
	$version =~ s/[\s\n\r\t]//g;
	$APPLES_FV = $version || 'unknown.s2000';
	return $APPLES_FV;
}

=head2 Return the serialization version for Serializable
	
	Returns:
	an APPLES version tag must have a substring like \.s[0-9]+ which 
	gives the version number for Serializable

	The last version from Subversion was 1610. We start at 2000, which is used
	if the version cannot be determined.

=cut
our $APPLES_SV = 0;
sub APPLES_SERIALIZATIONVERSION {
	our $APPLES_SV;
	return $APPLES_SV if $APPLES_SV > 0;
	my $ver = APPLES_FULLVERSION;
	Runtime::warn ("Unknown APPLES version $ver") if $ver =~ m/^unknown/;
	if ($ver =~ m/\.s([0-9]+)/) {
		$APPLES_SV = $1;
	} else {
		Runtime::warn ("Unable to get APPLES serialization version from $ver");
		$APPLES_SV = 2000;
	}
	return $APPLES_SV;
}

=head2 Dump a serializable object to stdout as JSON

 Parameters: 
 $obj : the object to dump

=cut

sub dump_json {
	eval "use Serialization::Serializable";
	my $obj = shift;
	print to_json (
		Serialization::Serializable::to_hash ($obj), 
			{allow_blessed => 1, pretty=>1});
}

1;
