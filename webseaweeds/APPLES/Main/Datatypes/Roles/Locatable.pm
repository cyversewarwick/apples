#!/usr/bin/perl

=head1 Helper role to extract sequence locations
from other objects

=cut

package Datatypes::Roles::Locatable;

use strict;

use feature ':5.10';

use Runtime;

use Moose::Role;
use MooseX::Declare;

use JSON;

use Serialization::Serializable;
use Datatypes::Sequence::Location;
use Sequences::Database::Relative_Location;

=head2 Return the location of $self

 Parameters:
 $self : a self object

 Returns:

 a Datatypes::Sequence::Location as accurately
 as possible for the current object, or undef
 if the object cannot be located

=cut

sub get_sequence_location {
	my $self = shift;

	if ( UNIVERSAL::isa( $self, "Datatypes::Sequence::Location" ) ) {
		return $self;
	}

	if ( UNIVERSAL::isa( $self, "Sequences::Genomic_Sequence" ) ) {
		my $dbi      = undef;
		my $location = undef;
		my $id       = undef;
		if ( Serialization::Serializable::is_hash( $self->{source} ) ) {
			$dbi      = _find_db_identifier( $self->{source}->{db} );
			$location = $self->{source}->{location};
		} else {
			my $sid = $self->id() || $self->accession_number;
			$location =
			  Sequences::Database::Relative_Location->new( identifier => $sid,
			  );
			$dbi = _find_db_identifier($location);
		}

		my $l = Datatypes::Sequence::Location->new(
			db       => $dbi,
			location => $location
		);
		$l->validate();
		return $l;
	}

	if ( UNIVERSAL::can( $self, "_get_location" ) ) {
		return $self->get_location();
	}

	die "Sequence location cannot be obtained for $self.";
}

=head2 Find a db identifier for a database object

 Parameters:
 $db : isa "Sequences::Database::Sequence"
 $lookup : try to validate locations

 Returns:
 a string identifier for get_sequence_database, or undef

=cut

sub _find_db_identifier {
	my $db = shift;
#	my $lookup = shift || 1;
## feature disabled, this takes too much time.
	my $lookup = 0;

	if ( UNIVERSAL::isa( $db, "Sequences::Database::Relative_Location" ) ) {

		if (!$lookup) {
			return _guess_db_by_identifier($db->identifier);
		}

		## we see if we can find a sequence database that contains this id
		my $ids   = get_sequence_database_ids();
		my $found = 0;
		my $dbi;
		foreach $dbi (@$ids) {
			my $db2 = get_sequence_database($dbi);

			eval { $db2->validate_location($db); };

			if ( !$@ ) {
				$found = 1;
				last;
			}
		}

		if ( !$found ) {
			return undef;
		} else {
			return $dbi;
		}
	} elsif ( UNIVERSAL::isa( $db, "Sequences::Database::Sequence" ) ) {
		## for sequence databases, check if we find one that
		## serializes canonically to the same value
		my $ids = get_sequence_database_ids();

		my $js1 = to_json(
			Serialization::Serializable::to_hash($db),
			{ canonical => 1, allow_blessed => 1 }
		);

		foreach my $id (@$ids) {
			my $db2 = get_sequence_database($id);

			my $js2 = to_json(
				Serialization::Serializable::to_hash($db2),
				{ canonical => 1, allow_blessed => 1 }
			);

			if ( $js2 eq $js1 ) {
				return $id;
			}
		}

	} else {
		## is it a string ?
		my $ddb = undef;
		eval { $ddb = get_sequence_database($db); };
		if ( defined($ddb) ) {
			return $db;
		} else {
			return undef;
		}
	}

	return undef;
}

=head2 Guess the sequence database by identifier name

 Parameters:
 $identifier : the identifier

 Returns:
 a sequence DB identifier (genbank, ensembl or ensemblgenomes)

=cut

sub _guess_db_by_identifier {
	my $identifier = shift;

	## stuff that starts with ENS is most likely ensembl
	if ($identifier =~ m/^ENS/) {
		## TODO: distinguish ensembl and ensemblgenomes
		return "ensembl";
	} else {
		return "genbank";
	}
}

1;
