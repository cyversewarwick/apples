#!/dev/null

package Runtime::Web_Session_Helpers;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_cgi_and_session print_header get_dbh
 check_access get_resultdb);

=head1 The ORANGES Web tool Runtime

Everything you need to set up an ORANGES CGI script

=cut

use CGI;
use CGI::Session;

use DBI;
use JSON;

use Configuration::ReadConfig qw(:ALL);
use Configuration::AppleSeeds;

# Initialize CGI session and such.
$CGI::POST_MAX = 1024 * 10000;

our $cgi     = undef;
our $session = undef;

=head2 Function get_cgi_and_session
Return a CGI object and session handle.
=cut

sub get_cgi_and_session() {
	our ( $cgi, $session );
	unless ( $cgi && $session ) {
		$cgi = CGI->new;
		my $sessionid = $cgi->param("CGISESSID") || undef;
		$session =
		  new CGI::Session( undef, $sessionid, { Directory => '/tmp' } );
		$session->expire('+1h');
	}
	return ( $cgi, $session );
}


our $header_was_printed = 0;

=head2 Function print_header
Print content type and cookie header when used as CGI service
=cut

sub print_header {
	return if $header_was_printed;

	# make sure we have a session created
	get_cgi_and_session();

	my $type = shift;

	# we generate JSON output by default.
	$type = 'application/json' if !defined($type);
	print $cgi->header( -type   => $type,
						-cookie => $cgi->cookie( CGISESSID => $session->id ) );
	$header_was_printed = 1;
}

our $dbh = undef;

=head2 Function get_dbh
Open a database connection and return the handle.
=cut

sub get_dbh() {
	unless ($dbh) {
		use Data::Dumper;
		use Configuration::ReadConfig qw(%config);
		my $dsn = "DBI:mysql:database=$database;host=$hostname;port=$port";
		$dbh = DBI->connect( $dsn, $user, $password )
			or die "Could not connect to database: $@ " . Dumper (\%config);
		$dbh->{RaiseError} = 1;
	}
	return $dbh;
}


=head2 Function check_access

Check if the current session allows database access.

If we cannot access the database, we print a JSON
error message and exit.

 Parameters:
 $quit_if_na : 0/undef => (optional) exit and print JSON error
               1 => don't exit, don't print

 Returns:
 a Runtime::ResultDB (if access was granted)

=cut

sub check_access {
	my $no_exit_if_denied = shift || 0;
	get_cgi_and_session();
	if ( defined( $session->param('info') ) ) {
		my $info = $session->param('info');
		if ( $info->{'access'} ne 'granted' ) {
		  DENIED:
			unless ($no_exit_if_denied) {
				print_header();
				my $result = { 'access' => "denied", };
				$session->delete();
				print encode_json($result);
				exit;
			}
			return undef;
		}
	} else {
		goto DENIED;
	}
	return get_resultdb();
}

=head2 Function get_resultdb

Create and return a result database object

 Returns:
 a Runtime::ResultDB

=cut

sub get_resultdb {
	eval "use Runtime::ResultDB;";
	if ($@) {
		die "Unable to load the ResultDB module: $@.";
	}

	my $userid = shift || $session->param("userid");
	my ( $cgi, $session ) = get_cgi_and_session;
	## we load this module here so
	## we don't get import errors when
	## using this in APPLES only
	## if anyone calls check_access without a
	## user id, of course there needs to be an error
	eval "require Runtime::ResultDB"
	  or die "Cannot find ResultDB module: $@";
	our $db = Runtime::ResultDB->new( get_dbh, $userid );
	return $db;
}

1;
