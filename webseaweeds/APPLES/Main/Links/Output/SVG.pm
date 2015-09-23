#!/usr/bin/perl

=head1 Class to render a links network as a dot file

=cut

package Links::Output::SVG;

use strict;

use Runtime;
use Configuration::AppleSeeds;

use File::Temp qw(tempfile);
use IPC::Open3;
use Symbol qw(gensym);
use IO::File;

use Links::Output::Dot;

use Scalar::Util qw(blessed);

=head2 Render an SVG of a Links_Database via dot 

 Parameters:
 $ldb : the database
 $wattr : (optional) the name of the weight attribute
 			default : "weight"
 
 Returns:
 the SVG code as a string

=cut

sub render {
	my $ldb = shift;
	my $wattr = shift || "weight";

	my ( $fh, $filename ) = tempfile();
	print $fh Links::Output::Dot::render($ldb, $wattr);
	close $fh;

	my $dot = find_executable("dot");

	local *CATCHERR = IO::File->new_tmpfile;
	my $pid = open3(gensym, \*CATCHOUT, ">&CATCHERR", "$dot -Tsvg $filename");
	my $ll = "";
	while( <CATCHOUT> ) {
		$ll.= $_;
	}
	waitpid($pid, 0);
	my $err = $?;
	my $errtext = "";
	seek CATCHERR, 0, 0;
	while( <CATCHERR> ) {
		$errtext.= $_;
	}

	if ($err != 0) {
		die "$err : $errtext";
	} 

	return $ll;
}

1;
