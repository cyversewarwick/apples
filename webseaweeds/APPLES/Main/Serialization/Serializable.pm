
=pod

=head1 Class Serializable

Serializable - Class for outputting various data structures to XML for re-use in webservices and other places.

=head1 SYNOPSIS
If a class wishes to be able to export itself and all the data it holds to an XML file, then it needs to extend Serializable.

=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<to_hash>

Takes all class attributes and puts them into hash ready to be written to XML file. N.B. all members starting 
with "ns__" will be skipped and not exported to XML file. 

=item C<from_hash>

NOT IMPLEMENTED. 
Reverse of to_hash(), reads data from the hash and constructs the class. Maybe should be implemented by each class?

=back
=head1 LICENSE

This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package and distributed 
under Academic Non-Commercial Use Licence.

=head1 COPYRIGHT

(c) Copyright University of Warwick 2009-2010

=head1 AUTHOR

=cut

package Serialization::Serializable;
use strict;
use feature qw(state);

use Moose;

use Scalar::Util qw(blessed reftype);
use Data::Dumper;
use Devel::Cycle;
use JSON;
use Data::UUID;
use MIME::Base64;

use Runtime;
use Carp;

=head2 Check if a given object is a HASHREF

 Parameters:
 $hash_obj : an object
 
 Returns : true if $hash_obj is a hashref
 
=cut

sub is_hash {
	my $hash_obj = shift;
	my $is_a_hash = ( ref($hash_obj) && ref($hash_obj) eq 'HASH' )
	  || ( reftype($hash_obj) && reftype($hash_obj) eq 'HASH' );
	return $is_a_hash;
}

=pod
    
=head2 from_hash()

Reverse of to_hash(), reads data from the hash and constructs the class. 

Override if your class needs special code to restore non-serialized fields.

This code checks if the version information matches the versions present in the
current code. 

Reading old version objects into newer version code is forbidden
and will cause the code to croak(). 

Reading new version objects into old version code will not cause a warning, but
this can be enabled by uncommenting the lines in the code below.

Serializable objects can provide a sub 

 sub serialization_restore {
 	## this will be called after an object is restored
 }

to execute after de-serialization.

Blessing may fail as well if the necessary packages have not been loaded.

=cut

sub from_hash {
	my $hash_obj = shift;
	my $ignore_versionid = shift || 0;

	if ( !defined($hash_obj) ) {
		return undef;
	}

	if ( _weakened_memory_cycle_exists($hash_obj) ) {
		confess( "Cannot restore cyclic object: " . Dumper($hash_obj) );
	}

	my $is_a_hash = ( ref($hash_obj) && ref($hash_obj) eq 'HASH' )
	  || ( reftype($hash_obj) && reftype($hash_obj) eq 'HASH' );

	if ( $is_a_hash
		&& defined( $hash_obj->{'SERIAL_VERSIONID'} ) )
	{
		my $object = {};
		while ( my ( $key, $value ) = each %$hash_obj ) {
			## member is also Serializable.
			if ( ref($value) eq "HASH"
				and defined( $value->{'SERIAL_VERSIONID'} ) )
			{
				$object->{$key} =
				  Serialization::Serializable::from_hash($value);
			} elsif ( ref($value) eq "HASH" )
			{    ## traverse into hashes and arrays
				$object->{$key} = {};
				while ( my ( $key2, $value2 ) = each %$value ) {
					$object->{$key}->{$key2} =
					  Serialization::Serializable::from_hash( $value2, 1 );
				}
			} elsif ( ref($value) eq "ARRAY" ) {
				$object->{$key} = [];
				foreach my $value2 (@$value) {
					push @{ $object->{$key} },
					  Serialization::Serializable::from_hash( $value2, 1 );
				}
			} else {

				# otherwise: this is an attribute. Just copy it.
				$object->{$key} = $value;
			}
		}
		my $class = $hash_obj->{'SERIAL_VERSIONID'};
		my $version;

		if ( $class =~ m/(.*)\:([0-9]*)$/ ) {
			$class   = $1;
			$version = $2;
		}

		eval "require $class;";
		## we don't check for errors here, since this rules out using locally declared classes

		bless $object, $class;
		if ( $version > 1 && $object->version() > $version ) {

			# we warn when restoring newer objects from older versions
			# warn( "Restoring newer $class object from older version $version" . " to "
			# 	  . $object->version() );
		} elsif ( $object->version() < $version ) {

			# we warn if we restore older objects from newer versions
			# warn( "Restored older $class object from newer version $version" . " to "
			# 	  . $object->version() );
		}

		if ( UNIVERSAL::can( $object, 'serialization_restore' ) ) {
			$object->serialization_restore();
		}

		return $object;
	} else {
		unless ($ignore_versionid) {
			confess("Serializable object cannot be restored from hash"
				  . " because version and object type are missing : "
				  . Dumper($hash_obj) . "\n"
				  . ref($hash_obj) . " : "
				  . $hash_obj->{'SERIAL_VERSIONID'} );
		} else {    # when restoring recursively, see what we have
			my $object;
			if ( ref($hash_obj) eq "HASH" ) { ## traverse into hashes and arrays
				$object = {};
				while ( my ( $key, $value ) = each %$hash_obj ) {
					$object->{$key} =
					  Serialization::Serializable::from_hash( $value, 1 );
				}
			} elsif ( ref($hash_obj) eq "ARRAY" ) {
				$object = [];
				foreach my $value (@$hash_obj) {
					push @{$object},
					  Serialization::Serializable::from_hash( $value, 1 );
				}
			} else {
				$object = $hash_obj;
			}
			return $object;
		}
	}
	return undef;
}

=pod
    
=head2 to_hash()

Takes all class attributes and puts them into hash ready to be written to XML/JSON file. 
N.B. all members starting 
with "ns__" will be skipped.

Returns a hash ready to be written by XML_Writer or JSON.

=cut

sub to_hash {
	my $self_obj = shift;
	my $sub_call = shift || 1;
	my $object   = {};

	if ( !defined($self_obj)
		|| UNIVERSAL::isa( $self_obj, "Serialization::Non_Serializable" ) )
	{
		return undef;
	}

	if (   is_hash($self_obj)
		&& UNIVERSAL::isa( $self_obj, 'Serialization::Serializable' ) )
	{
		## hack to make sure we use overridden to_hash functions
		if ( $sub_call == 1 ) {
			$object = $self_obj->to_hash(2);
		} else {
			while ( my ( $key, $value ) = each %$self_obj ) {
				unless ( $key =~ m/^ns__/ )
				{    # filter bits which should not be serialized
					$object->{$key} = to_hash( $value, 1 );
				}
			}
		}
		$object->{'SERIAL_VERSIONID'} =
		  blessed($self_obj) . ":" . $self_obj->version();
	} else {
		if ( ref($self_obj) eq "ARRAY" ) {
			$object = [];
			foreach (@$self_obj) {
				push @$object, to_hash( $_, 1 );
			}
		} elsif ( is_hash($self_obj) ) {
			while ( my ( $key, $value ) = each %$self_obj ) {
				$object->{$key} = to_hash( $value, 1 );
			}
		} else {
			$object = "$self_obj";
		}
	}

	return $object;
}

=head2 Return a unique string for this object

 Parameters:
 $self : a self object
 
 Returns:
 a string that is a unique representation of this object (i.e. can be 
 used for comparison)

 Note that you should override this method with something that is more 
 sensible than the code given here, and -- most importantly -- is carried
 through serialization. 
 
=cut

sub unique_id {
	my $self = shift;
	if ( !defined( $self->{ns__serializable_uniqueid} ) ) {
		my $ug = Data::UUID->new;
		return $ug->create_str();
	}

	return $self->{ns__serializable_uniqueid};
}

=head2 Compare $self to another serializeable object

Serializable objects can be compared by comparing their
canonical JSON representation.

 Parameters:
 $self : a self object
 $other : another serializable object
 
 Returns:
 1 if the objects are equal
 0 otherwise

=cut

sub equals {
	my $self  = shift;
	my $other = shift;

	if (
		UNIVERSAL::isa( $self,  "Serialization::Serializable" ),
		UNIVERSAL::isa( $other, "Serialization::Serializable" )
	  )
	{
		my $js1 = $self->unique_id;
		my $js2 = $self->unique_id;

		return $js1 eq $js2 ? 1 : 0;
	} else {
		confess("Cannot compare non-serializable objects");
	}
}

=head2 Create a cloned copy of $self

This works by means of serializing to a string, and then deserializing again.

 Parameters:
 $self : a self object
 
 Returns:
 a copy of $self

=cut

sub cloneme {
	my $self = shift;
	my $str = to_json (to_hash ( $self ), {allow_blessed=>1});
	return from_hash (from_json($str));
}

=head2 version()

return a version for identifying when fields have been added or changed. Every datatype implementing
Serializable needs to increment this whenever its specification changes. Alternatively, this can
be coupled with subversion version numbers as shown below. 

Must return an integer.

=cut

sub version {
	return APPLES_SERIALIZATIONVERSION;
}

=head2 Helper function to check for cycles

=cut

sub _weakened_memory_cycle_exists {
	my $obj = shift;

	my $found = 0;
	find_weakened_cycle( $obj, sub { $found = 1; } );
	return $found;
}

1;
