#!/usr/bin/perl

=head1 Class Sequences::Database::Sequence;

Generic sequence database class interface

=cut

package Sequences::Database::Sequence;

use strict;

use Runtime;
use Configuration::AppleSeeds;
use Carp;
use Scalar::Util qw(blessed);
use Data::Dumper;

our @ISA = qw(Serialization::Serializable);
require Serialization::Serializable;

use Serialization::Serializable_Array;

use Sequences::Database::Relative_Location;

use Runtime;
use Carp;

=head2 Get a sequence by its accession from the database

Parameters:
$self : A Database
$location : a Relative_Location

Returns :
Nothing, but dies if location is invalid

=cut

sub validate_location {
	my $self = shift;
	my $loc  = shift;

	## override this to also check if identifier exists
	$loc->validate();

	return;
}

=head2 Get a sequence by its accession from the database

Parameters:
$self : A Database
$acc  : the accession string

Returns :
The sequence as a Genomic_Sequence

=cut

sub get_sequence {
	my $self = shift;
	my $acc  = shift;

	my $location =
	  Sequences::Database::Relative_Location->new( identifier => $acc );
    
	return $self->get_sequence_by_location($location);
}

=head2 Get sequences by sequence location

Location objects specify locations of sequences relative to
one or more accessions.

 Parameters:
 $self : a Database
 $location : a Relative_Location

 Returns:
 An ARRAYREF [Genomic_Sequence]

=cut

sub get_sequence_by_location {
	my $self   = shift;
	my $loc    = shift;
	my @result = ();

	unless ( UNIVERSAL::isa( $loc, "Sequences::Database::Relative_Location" ) )
	{
		confess("Invalid location record.");
	}

	my $cache = undef;
    
	unless ( exists( $self->{no_cache} ) )
	{
		$cache = cache;
	}

	my $cache_location = bless {
		location => $loc,
		db       => $self,
	  },

	  'Serialization::Serializable';

	my $r = ( defined($cache) ) ? $cache->cache_get($cache_location) : undef;

	if ( UNIVERSAL::isa( $r, 'Serialization::Serializable_Array' ) ) {
		$r = $r->data;
	}

	if ( !defined($r) ) {
		info( "Getting sequence " . $loc->identifier . " from DB." );
		# print "\nSequence.pm line 117. get_sequence_by_location, call _get_sequence_by_location, Getting sequence " . $loc->identifier . " from DB.\n";
		$r = $self->_get_sequence_by_location( $loc, );
		# print "\nSequence.pm line 119. get_sequence_by_location, returned from _get_sequence_by_location.\n";
		if (UNIVERSAL::can ($self, "get_meta_data")) {
			$r->{metadata} = $self->get_meta_data ($loc->identifier);
		}

		if ( defined($cache) && defined($r) ) {
			if ( ref($r) eq 'ARRAY' ) {
				my $s = Serialization::Serializable_Array->new;
				$s->data($r);
				$cache->cache_put( $cache_location, $s );
			} else {
				$cache->cache_put( $cache_location, $r );
			}
		}
	} else {
		info( "Sequence " . $loc->identifier . " was in the cache." );
	}

	if ( defined($r) ) {
		if ( ref($r) eq 'ARRAY' ) {
			foreach my $rr (@$r) {
				push @result, $rr;
			}
		} else {
			push @result, $r;
		}
	}
	# print "\nSequence.pm line 146. get_sequence_by_location, exiting.\n";
	return \@result;
}

=head2 Get a single sequence by sequence location

Location objects specify locations of sequences relative to
one or more accessions.

 Parameters:
 $self : a Database
 $location : a Location

 Returns:
 An ARRAYREF [Genomic_Sequence]

=cut

sub _get_sequence_by_location {
	confess("Abstract method called");
}

=head2 Get a list of all accessions in this database

 Parameters:
 $self : a self object
 $species : a species/object filter

 Returns:
 an ARRAYREF containing all the accessions

=cut

sub get_all_accessions {
	my $self = shift;
	die "Method not implemented by sequence database class " . blessed ($self);
}

=head2 Overloaded from Serializable

=cut

sub to_hash {
	my $self = shift;

	return { 'name' => $self->name, };
}

=head2 Overloaded from Serializable

=cut

sub serialization_restore {
	confess("Sequence_Database cannot be deserialized.");
}

1;
