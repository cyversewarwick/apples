package Runtime::ResultDB;

=head1 Results database

Dataset-based simple JSON database.

=cut

use strict;

use Compress::Zlib;
use XML::Simple;
use JSON;
use File::Temp qw(tempfile);
use Digest::MD5 qw(md5_hex);
use POSIX qw(strftime);
use Data::Dumper;
use Scalar::Util qw(blessed);

BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
}

use Configuration::AppleSeeds;

use Runtime;
use Carp;

use DateTime;
use DateTime::Format::DBI;

use Serialization::Serializable;

=head2 ResultDB constructor

 Parameters:
 $class    : Runtime::ResultDB
 $dbh      : a DBI handle
 $user_id  : the user id to use

=cut

sub new {
	my ( $class, $dbh, $user_id ) = @_;
	my $self = {
		dbh      => $dbh,
		user_id  => $user_id,
	};

	bless $self, $class;
	return $self;
}

=head2 Check if we can access a data item

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_di  :  the ID of the data item

 Returns:
 String "access : <val>" with the value of the access value
 in the corresponding dataset.

 undef if the dataset wasn't found.

=cut

sub check_access_to_di($$) {
	my ( $self, $id_di ) = @_;
	my $sth = $self->{'dbh'}->prepare(

		"SELECT `access` FROM `datasets`, `dataitems` WHERE "

		  # param 1: id_user
		  . " `id_user` = ?" . " AND `datasets`.`id_ds` = `dataitems`.`id_ds`" .

		  # param 2: id_di
		  " AND `dataitems`.`id_di` = ? "
	);
	$sth->execute( $self->{'user_id'}, $id_di );
	if ( my $ref = $sth->fetchrow_hashref() ) {
		return "access: " . $ref->{'access'};
	}
	return 0;
}

=head2 Check if we can access a dataset

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the ID of the dataset

 Returns:
 String "access : <val>" with the value of the access value
 in the dataset.

 undef if the dataset wasn't found.

=cut

sub check_access_to_ds($$) {
	my ( $self, $id_ds ) = @_;
	my $sth = $self->{'dbh'}->prepare(

		# param 1: id_user
		"SELECT `access` FROM `datasets` WHERE `id_user` = ?" .

		  # param 2: id_ds
		  " AND `id_ds` = ? "
	);
	$sth->execute( $self->{'user_id'}, $id_ds );
	if ( my $ref = $sth->fetchrow_hashref() ) {
		return "access: " . $ref->{'access'};
	}
	return 0;
}

=head2 Retrieve datasets

 Parameters:
 $self   :  an Runtime::ResultDB
 $dsname :  the name of the dataset,
            or an SQL match expression

            if $dsname contains the percent character,
            it will be used in an SQL LIKE expression

            otherwise, we compare for equality

 Returns:
 an ARRAYREF of matching datasets, for each dataset, we return a
 HASHREF with all the values associated with the dataset, and additional
 values
 {
 	identifier => <ident>,
 	id_user => <creating user's id>,
 	access => <access>
 	description => <a description of the dataset>
 }

=cut

sub retrieve_datasets($$) {
	my ( $self, $dsname ) = @_;

	my $sth;

	# check if we have wildcards
	if ( $dsname =~ m/\%/ ) {
		$sth =
		  $self->{'dbh'}->prepare(
			"SELECT * FROM `datasets` WHERE id_user=? AND `identifier` LIKE ?");
	} else {
		$sth =
		  $self->{'dbh'}->prepare(
			"SELECT * FROM `datasets` WHERE `id_user`=? AND `identifier`=?");
	}

	my @results = ();
	$sth->execute( $self->{'user_id'}, $dsname );
	while ( my $ref = $sth->fetchrow_hashref() ) {
		my %ds = %$ref;
		push @results, \%ds;
	}
	return \@results;
}

=head2 Retrieve dataset by ID

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  a dataset id

 Returns:
 the dataset's row from the SQL database,
 or undef if no such dataset exists.

=cut

sub retrieve_dataset_by_id($$) {
	my ( $self, $id_ds ) = @_;

	my $sth =
	  $self->{'dbh'}
	  ->prepare("SELECT * FROM `datasets` WHERE id_user=? AND `id_ds`=?");

	my @results = ();
	$sth->execute( $self->{'user_id'}, $id_ds );
	if ( my $ref = $sth->fetchrow_hashref() ) {
		my %ds = %$ref;
		return \%ds;
	}
	return undef;
}

=head2 Retrieve dataset item IDs

This function retrieves all item ids for a dataset.

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the dataset id

 Returns:
 an ARRAYREF of HASHREF's :
 [
  	{
  		id_di => <the unique dataitem id>,
  		source => <the source of the data item>,
  		type => <the type of data stored>,
  		created => <the date when the item was created>,
  	}
 ]

=cut

sub retrieve_dataitem_ids($$) {
	my ( $self, $id_ds ) = @_;

	unless ( $self->check_access_to_ds($id_ds) =~ m/^access.*/ ) {
		die("access denied to dataset $id_ds for current user");
	}

	my $sth =
	  $self->{'dbh'}->prepare(
"SELECT `id_di`, `type`, `created`, `source` FROM `dataitems` WHERE `id_ds` = ?"
	  );

	my @results = ();
	$sth->execute($id_ds)
	  or die("Failed to retrieve dataitem ids for $id_ds");
	while ( my $ref = $sth->fetchrow_hashref() ) {
		push @results, $ref;
	}
	return \@results;
}

=head2 Retrieve dataset item IDs, filter by source

This function retrieves all item ids for a dataset
which match a certain source.

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the dataset id

 Returns:
 an ARRAYREF of HASHREF's :
 [
  	{
  		id_di => <the unique dataitem id>,
  		source => <the source of the data item>,
  		type => <the type of data stored>,
  		created => <the date when the item was created>,
  	}
 ]

=cut

sub retrieve_dataitem_ids_by_source($$$) {
	my ( $self, $id_ds, $source ) = @_;

	unless ( $self->check_access_to_ds($id_ds) =~ m/^access.*/ ) {
		die("access denied to dataset $id_ds for current user");
	}

	my $sth =
	  $self->{'dbh'}
	  ->prepare( "SELECT `id_di`, `type`, `created`, `source` FROM `dataitems`"
		  . " WHERE `id_ds` = ? AND `source` = ?" );

	my @results = ();
	$sth->execute( $id_ds, $source )
	  or die("Failed to retrieve dataitem ids for $id_ds/$source");
	while ( my $ref = $sth->fetchrow_hashref() ) {
		push @results, $ref;
	}
	return \@results;
}


=head2 Retrieve dataset item IDs, filter by source

This function retrieves all item ids for a dataset
which are LIKE a certain source.

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the dataset id

 Returns:
 an ARRAYREF of HASHREF's :
 [
  	{
  		id_di => <the unique dataitem id>,
  		source => <the source of the data item>,
  		type => <the type of data stored>,
  		created => <the date when the item was created>,
  	}
 ]

=cut

sub retrieve_dataitem_ids_by_source_similar($$$) {
	my ( $self, $id_ds, $source ) = @_;

	unless ( $self->check_access_to_ds($id_ds) =~ m/^access.*/ ) {
		die("access denied to dataset $id_ds for current user");
	}

	my $sth =
	  $self->{'dbh'}
	  ->prepare( "SELECT `id_di`, `type`, `created`, `source` FROM `dataitems`"
		  . " WHERE `id_ds` = ? AND `source` LIKE ?" );

	my @results = ();
	$sth->execute( $id_ds, $source )
	  or die("Failed to retrieve dataitem ids for $id_ds/$source");
	while ( my $ref = $sth->fetchrow_hashref() ) {
		push @results, $ref;
	}
	return \@results;
}


=head2 Create a new dataset

 Parameters:
 $self        :  an Runtime::ResultDB
 $name        :  the name for the dataset
 $description :  a textual description of the dataset

 Returns:
 $id_ds  :  the dataset id

=cut

sub create_dataset($$$) {
	my ( $self, $name, $description ) = @_;

	my $sth =
	  $self->{'dbh'}->prepare(
"INSERT INTO `datasets` (`identifier`, `description`, `id_user`, `access`) VALUE (?,?,?,?)"
	  );
	$sth->execute( $name, $description, $self->{'user_id'}, 0 )
	  or die("Could not create dataset $name");

	my $id_ds =
	  $self->{'dbh'}->last_insert_id( undef, undef, qw(datasets id_ds) );
	return $id_ds;
}

=head2 Delete a dataset

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the dataset id

 Returns:
 Nothing.

=cut

sub delete_dataset($$) {
	my ( $self, $id_ds ) = @_;
	if ( $self->check_access_to_ds($id_ds) =~ m/^access.*/ ) {
		my $sth =
		  $self->{'dbh'}->prepare(
			"DELETE FROM `datasets` where `id_ds` = ? AND `id_user` = ?");
		$sth->execute( $id_ds, $self->{'user_id'} )
		  or die("Could not delete dataset $id_ds");
	} else {
		die("cannot delete $id_ds : access denied");
	}

}

=head2 Retrieve the data for a given data item

The data associated with each data item is stored in a
file (to cater for large data objects).

This sub retrieves the data from the file system.

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_di  :  the data item id

 Returns:
 {
		'type' => data type,
		'size' => the length of the data,
		'data' => the data retrieved from the file
		          (for JSON/XML, the data will be decoded
		           into a hash)
 }

=cut

sub retrieve_dataitem_data($$) {
	my ( $self, $id_di ) = @_;

	unless ( $self->check_access_to_di($id_di) =~ m/^access.*/ ) {
		die("access denied to dataset $id_di for current user");
	}

	debug("getting $id_di from db");
	my $sth =
	  $self->{'dbh'}->prepare("SELECT * FROM `dataitems` WHERE `id_di` = ?");

	$sth->execute($id_di);

	my $result = {};
	my $type;
	my $size;
	if ( my $ref = $sth->fetchrow_hashref() ) {

		$type = $ref->{'type'};

		my $sth0 = $self->{'dbh'}->prepare("SELECT data FROM blobstore WHERE id_blob = ? ORDER BY chunk")
		  ;
		my $numrows = $sth0->execute ($ref->{'data'});
		if ($numrows < 1) {
			die "Failed to get data item $id_di (blob id $ref->{data} ) : $numrows rows found.";
		}

		my $data = '';
		while ( my $ref2     = $sth0->fetchrow_hashref ) {
			$data .= $$ref2{'data'};
		}
		
		$sth0->finish;

		if ( $type =~ m/^XML/ ) {
			if ( $type =~ m/\.gz$/ ) {
				$data = Compress::Zlib::memGunzip($data);
			}
			$result = XML::Simple::XMLin($data);
			$type   = 'XML';
		} elsif ( $type =~ m/^JSON/ ) {
			if ( $type =~ m/\.gz$/ ) {
				$data = Compress::Zlib::memGunzip($data);
			}
			eval { $result = decode_json($data); };
			if ($@) {
				die(
					"Could not decode (supposedly) JSON text $data in item "
					  . Dumper($ref) );
			}
			$type = 'JSON';
		} else {
			$result->{'data'} = $data;
		}
	} else {
		die("data item $id_di cannot be found");
	}
	my $retval = {
		'type' => $type,
		'size' => $size,
		'data' => $result,
	};

	return $retval;
}

=head2 Store a new data item

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_ds  :  the dataset id
 $type   :  the type of the data item
 $source :  an identifier of the source of the data item
 $data   :  the contents of this element will be stored
            on the disk

 Returns:
 $id_di : the unique ID for the data item

=cut

sub store_dataitem($$$$$) {
	my ( $self, $id_ds, $type, $source, $data ) = @_;
	
	my $chunk_length = 256*1024;

	unless ( $self->check_access_to_ds($id_ds) =~ m/^access.*/ ) {
		die("access denied to dataset $id_ds for current user");
	}

	my $sth0 =
	  $self->{'dbh'}->prepare("INSERT INTO blobstore (`data`, `chunk`) VALUE (?, ?)");

	my 	($fh, $filename) = tempfile();
	binmode ($fh);
	syswrite ($fh, $data, length ($data));
	close $fh;
		
	open TEMPFILE, "<", $filename;
	binmode (TEMPFILE);
	my $chunk;
	
	my $bytes_read = sysread (TEMPFILE, $chunk, $chunk_length);
	
	$sth0->execute($chunk, 0)
	  or die
	  "Could not insert data for dataset $id_ds / $source into database.";
	
	$sth0->finish();

	$filename =
	  $self->{'dbh'}->last_insert_id( undef, undef, qw(blobstore id_blob) );

	my $sth0a =
	  $self->{'dbh'}->prepare("INSERT INTO blobstore (id_blob, data, chunk) VALUE ( ? , ? , ? )")
		or die "Could not insert data for dataset $id_ds / $source into database.";
	
	$bytes_read = sysread (TEMPFILE, $chunk, $chunk_length);
	my $chunkid = 1;
	while ($bytes_read > 0) {
		$sth0a->execute ( $filename, $chunk, $chunkid )
			or die
	  	"Could not insert data for dataset $id_ds / $source into database.";
	  	$chunkid++;
		$bytes_read = sysread (TEMPFILE, $chunk, $chunk_length);
	}
	
	close TEMPFILE;
	$sth0a->finish();

	my $sth1 =
	  $self->{'dbh'}->prepare(
"INSERT INTO dataitems (`id_ds`, `created`, `type`, `source`, `data`) VALUE (?, NOW(), ? , ? , ?)"
	  );
	$sth1->execute( $id_ds, $type, $source, $filename )
	  or die("Failed to create a new data item.");
	my $id_di =
	  $self->{'dbh'}->last_insert_id( undef, undef, qw(dataitems id_di) );

	return $id_di;
}

=head2 Delete a data item

 Parameters:
 $self   :  an Runtime::ResultDB
 $id_di  :  the data item id

 Returns:
 Nothing

=cut

sub delete_dataitem($$) {
	my ( $self, $id_di ) = @_;
	my $access = $self->check_access_to_di($id_di);
	if ( $access =~ m/^access.*/ ) {
		my $sth =
		  $self->{'dbh'}->prepare("DELETE FROM `dataitems` WHERE `id_di` = ?");
		$sth->execute($id_di)
		  or die("cannot delete data item $id_di");
	} else {
		die("cannot delete $id_di : access denied ($access)");
	}
}

=head2 Get a dataset hash, create if doesn't exists

 Parameters:
 $self   :  an Runtime::ResultDB
 $ds     :  dataset name string

 Returns:
 the result of L<retrieve_dataset_by_id> for the
 (possibly newly created) dataset

=cut

sub get_or_create_dataset($$) {
	my ( $self, $ds ) = @_;

	my $result = undef;

	my $qdatasets = $self->retrieve_datasets($ds);

	$result = $qdatasets->[0]
	  if ( scalar @$qdatasets ) > 0;

	warn(   "Multiple datasets matching $ds were found. Let's hope "
		  . Dumper($result)
		  . " is the one we wanted..." )
	  if ( scalar @$qdatasets ) > 1;

	if ( !defined($result) ) {
		debug("Creating new dataset $ds.");

		my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
		my $id_ds =
		  $self->create_dataset( $ds, "Dataset created at " . $now_string );
		debug("Inserted dataset has id $id_ds");
		$result = $self->retrieve_dataset_by_id($id_ds);
		if ( !defined($result) ) {
			die( "Failed to create a new dataset $ds, having "
				  . Dumper($result) );
		}
	}

	return $result;
}

=head2 This is a convenience function for getting a
data item using time-stamped versioning

If the source points to a serializable class, we will
deserialize the data item.

Otherwise, we will return the item just as returned by ResultDB.

When different versions exist, only the most recent item for this source
is returned.

 Parameters:
 $self : a self object
 $dataset : a dataset name
 $source : the source id of the data item
 $id_di  : (optional) the id of the data item (this would be sufficient to
           query uniquely, but we also check source and dataset
           here)

 Returns:
 a single item or undef if this item does not exist yet

=cut

sub get_versioned_dataitem {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;
	my $id_di   = shift;

	my $id_ds = $self->get_or_create_dataset($dataset)->{id_ds};
	my $dis = $self->retrieve_dataitem_ids_by_source( $id_ds, $source );

	my $data = undef;

	my $last_created = undef;

	my $db_parser = DateTime::Format::DBI->new( $self->{'dbh'} );

	if ( defined($id_di) && $id_di ne "" ) {
		## if we know the id_di, check it's in the set of dataitems
		## with the correct source and dataset
		my $found = 0;
		foreach my $di (@$dis) {
			if ( int( $di->{id_di} ) == int($id_di) ) {
				$found = 1;
				last;
			}
		}
		if ( !$found ) {
			die("No such dataset/dataitem/source found.");
		}
	} else {
		## don't know a data item id? find the most recent one
		foreach my $di (@$dis) {
			my $creation_date = $db_parser->parse_datetime( $di->{created} );

			## use this item either if it's newer or if we haven't found one yet
			if ( !defined($last_created) ) {
				$last_created = $creation_date;
				$id_di        = $di->{id_di};
			} elsif ( DateTime->compare( $creation_date, $last_created ) > 0 ) {
				$last_created = $creation_date;
				$id_di        = $di->{id_di};
			}
		}
	}

	if ( defined($id_di) && $id_di ne "" ) {
		$data = $self->retrieve_dataitem_data($id_di);
		$data->{id_di} = $id_di;
	}
	return $data;
}

=head2 List the versions of a versioned data item

 Parameters:
 $self : a self object
 $dataset : a dataset name
 $source : the source id of the data item

 Returns:
 a list of data item ids in an ARRAYREF

=cut

sub list_versions {
	my $self    = shift;
	my $dataset = shift;
	my $source  = shift;

	my $id_ds = $self->get_or_create_dataset($dataset)->{id_ds};
	my $dis = $self->retrieve_dataitem_ids_by_source( $id_ds, $source );

	my @res = map { int( $_->{id_di} ) } @$dis;

	return \@res;
}

=head2 Perform DB maintenance, i.e. remove orphaned blobs and dataitems

 Parameters:
 $self : a self object

=cut

sub db_maintenance {
	my $self = shift;

	my $sth0 =
	  $self->{'dbh'}->prepare("delete from `dataitems` where  `dataitems`.`id_ds` not in (select `datasets`.`id_ds` from `datasets`)")
		or die "Could not delete orphaned items (1).";
	;
	$sth0->execute ()
		or die "Could not delete orphaned items (2).";
	
	my $sth1 =
	  $self->{'dbh'}->prepare("delete `blobstore` from `blobstore` left outer join `dataitems` on `blobstore`.`id_blob` = `dataitems`.`data` where id_di is null;")
		or die "Could not delete orphaned blobs (1).";
	$sth1->execute ()
		or die "Could not delete orphaned blobs (2).";
}

1;
