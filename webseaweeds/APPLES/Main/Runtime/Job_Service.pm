#!/usr/bin/perl

=head1 Class Runtime::Job_Service

Dummy job service implementation. Runs everything sequentially on the
local system.

=cut

package Runtime::Job_Service;

use strict;

use Runtime;
use Carp;
use Time::HiRes qw (stat gettimeofday);

use JSON;
use Cwd;
use File::Spec;

=head2 Constructor

 Parameters:
 $class : Runtime::Job_Service
 $outputdir : the output directory ( defaults to cwd )

 Returns: 
 A new Job_Service
 
=cut

sub new {
	my $class = shift;
	my $outputdir = shift || getcwd;

	if ( !-d $outputdir || !-w $outputdir ) {
		warn("Output directory $outputdir is not writable, using cwd.");
		$outputdir = getcwd;
	}

	if ( !-d $outputdir || !-w $outputdir ) {
		confess("Cannot access the directory for storing results: $outputdir.");
	}

	return bless { outputdir => $outputdir, }, $class;
}

=head2 Run a computation

 Parameters:
 $self : a self object
 $comp : a Runtime::Computation
 
 Returns : 
 the result of the computation's run
 function.
 
=cut

sub run_computation {
	my $self = shift;
	my $comp = shift;
	return $comp->run;
}

=head2 Get a data item from the result database

Jobs can retrieve data items from a result database. These data items
are associated with a dataset ID.

In this implementation, we will use a flat file system database for this.

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
	my $version = shift;

	( $version, $type ) =
	  $self->_mostrecentversion( $dataset, $source, $version );
	my $data = undef;

	{
		open FILE, "<",
		  $self->_makefilename( $dataset, $source, $type, $version )
		  or die
		  "Data item $dataset / $source / $version / $type could not be found.";
		local $/ = undef;
		$data = <FILE>;
		close FILE;
	}

	if ( $type =~ m/json/i ) {
		$data = from_json($data);
	}

	return {
		type => $type,
		data => $data,
	};
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

	my $versionno = 1;
	my $version;
	my $filename;

	# make a unique filename
	do {

		my ( $seconds, $ustime ) = gettimeofday;
		$version = "${versionno}_$ustime";
		$filename = $self->_makefilename( $dataset, $source, $type, $version );
		++$versionno;
	} while ( -e $filename );

	debug("Writing result file : ${filename}");
	open OUT, ">", $filename;

	## is it ascii?
	unless ( $type =~ m/JSON/i
		|| $type =~ m/XML/
		|| $type =~ m/text/i )
	{
		binmode OUT;
	}

	## make hashes into JSON
	if ( $type =~ m/JSON/i && Serialization::Serializable::is_hash($item) )
	{
		$item = to_json( $item, { allow_blessed => 1 } );
	}

	print OUT $item;

	close OUT;
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
	my $type    = '';

	$dataset =~ s/[^A-Za-z0-9\._]/_/g;
	$source  =~ s/[^A-Za-z0-9\._]/_/g;

	my $datadir = File::Spec->catdir( $self->{outputdir}, $dataset );

	## dataset doesn't exist -> no versions exist either
	if ( !-d $datadir ) {
		return [];
	}

	my $name  = $self->_makefilename( $dataset, $source, $type, '*' );
	my @files = glob $name;
	my $vre   = $self->_makefilename( $dataset, $source, $type, '(.*)' );
	$vre =~ s/^$datadir//;
	map {
		$_ =~ s/^$datadir//;
		$_ =~ s/$vre$/$1/;
		$_ =~ s/\..+$//;
	} @files;

	return \@files;
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

	my $type = '';

	my $versions = $self->get_dataitem_versions( $dataset, $source );

	$dataset =~ s/[^A-Za-z0-9\._]/_/g;
	$source  =~ s/[^A-Za-z0-9\._]/_/g;

	my $datadir = File::Spec->catdir( $self->{outputdir}, $dataset );
	my $mostrecent;

	( $mostrecent, $type ) = $self->_mostrecentversion( $dataset, $source );

	if ( scalar @$rversions == 0 ) {
		@{$rversions} = grep { $_ ne $mostrecent } @$versions;
	}

	foreach my $ver (@$rversions) {
		my $fn = $self->_makefilename( $dataset, $source, $type, $ver );
		debug("Element $fn : purging  $ver / $type (MR is $mostrecent)");
		unlink($fn);
	}
}

=head2 Make a file name
  
 Parameters:
 $self    : a self object
 $dataset : a dataset name (folder)
 $source  : the dataitem source name
 $type    : type of the item  

 Returns: 
 The full filename of this data item
 
=cut

sub _makefilename {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;
	my $type    = shift;
	my $version = shift;

	$dataset =~ s/[^A-Za-z0-9\._]/_/g;
	$source  =~ s/[^A-Za-z0-9\._]/_/g;

	my $datadir = File::Spec->catdir( $self->{outputdir}, $dataset );
	if ( !-d $datadir ) {
		mkdir $datadir;
	}

	my $datafilename = File::Spec->catfile( $datadir, $source );
	my $filetype = "";

	if ( $type =~ m/[\/\.]?([A-Za-z0-9]+)$/ ) {
		$filetype = ".$1";
	} elsif ( $type =~ m/json/i ) {
		$filetype = ".json";
	} elsif ( $type =~ m/xml/i ) {
		$filetype = ".xml";
	}

	return "${datafilename}_v${version}${filetype}";
}

=head2 Get the most recent version/type of a data item

 Parameters:
 $self    : a self object
 $dataset : a dataset name (folder)
 $source  : the dataitem source name
 $fver    : (optional) the version to retrieve (undef to get the most recent)
 
 Returns:
 The version string of the most recent version object and its type as
 ($version, $type)

=cut

sub _mostrecentversion {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;
	my $fver    = shift;

	my $mrtype     = undef;
	my $mrusver    = '';
	my $mostrecent = undef;
	my $mrtime     = undef;

	my $versions = $self->get_dataitem_versions( $dataset, $source );

	my @vertypes;
	foreach my $ver (@$versions) {
		my $fn = $self->_makefilename( $dataset, $source, '', "$ver.*" );
		my @files = glob $fn;
		map {
			$_ =~ s/.*($ver\..+?)$/$1/;
			push @vertypes, $_;
		} @files;
	}

	foreach my $vertype (@vertypes) {
		my ( $ver, $type ) = split /\./, $vertype;
		my $fn = $self->_makefilename( $dataset, $source, $type, $ver );

		my $usver = -1;
		if ( $ver =~ m/([0-9]+)\_(.*)/ ) {
			$usver = $2;
		}

		my $tme = ( stat($fn) )[10];
		if (   !defined($mrtime)
			|| $tme > $mrtime
			|| ( $tme == $mrtime && ( $usver cmp $mrusver ) > 0 ) )
		{
			$mrtype     = $type;
			$mrtime     = $tme;
			$mostrecent = $ver;
		}
		## check if we have the required version
		if ( ( defined($fver) ) && $ver =~ m/^$fver(\..+)?/ ) {
			$mostrecent = $ver;
			$mrtype     = $type;
			$mrtime     = $tme;
			last;
		}
	}

	debug(
"Most recent version for $dataset / $source is $mostrecent.$mrtype (time $mrtime)"
	);

	return ( $mostrecent, $mrtype );
}
1;
