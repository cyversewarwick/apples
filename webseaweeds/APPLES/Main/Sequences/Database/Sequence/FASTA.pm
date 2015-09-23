#!/usr/bin/perl

=head1 Class Sequence_Database_FASTA;

Implementation of Sequence_Database interface for FASTA files.

=cut

package Sequences::Database::Sequence::FASTA;

use strict;

use Sequences::Genomic_Sequence;

use Sequences::Database::Sequence;
our @ISA = qw(Sequences::Database::Sequence);

use Runtime;
use Carp;

use Bio::Index::Fasta;
use Bio::DB::SeqI;

=head2 Constructor

 Parameters:
 
 $class : Sequence_Database_FASTA
 $fastaname : name of the underlying FASTA file
 
=cut

sub new {
	my ( $class, $fastaname ) = @_;
	die "File $fastaname not found." unless -e $fastaname;

	# make an index
	unless ( -e "$fastaname.fastaindex" ) {
		my $inx = Bio::Index::Fasta->new(
			-filename   => "$fastaname.fastaindex",
			-write_flag => 1
		);

		$inx->make_index($fastaname);
	}

	my $inx = Bio::Index::Fasta->new( -filename => "$fastaname.fastaindex" );

	my $self = {
		fastaname => $fastaname,
		sindex    => $inx,
		no_cache  => 1,
	};
	bless $self, $class;
	return $self;
}

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
	die "FASTA databases only support simple locations." 
		unless $loc->is_simple();
	
	my $acc = $loc->identifier ();
	my $seq = $self->{sindex}->fetch($acc);

	if ( !defined $seq ) {
		die "$acc was not found in " . $self->name;
	}
	return;
}

=head2 Get a  single sequence by sequence location

Location objects specify locations of sequences relative to 
one or more accessions. Notice that this implementation 
will ignore upstream/downstream arguments, as these don't
make sense with FASTA inputs. 

 Parameters:
 $self : a Database
 $location : a Location
 
 Returns:
 An ARRAYREF [Genomic_Sequence]

=cut

sub _get_sequence_by_location {
	my $self     = shift;
	my $location = shift;

	my $acc = $location->identifier;
	if ( !$location->is_simple() ) {
		warn(
			"FASTA does not support retrieval for complex locations."
		);
	}

	my $seq = $self->{sindex}->fetch($acc);

	if ( !defined $seq ) {
		die "$acc was not found in " . $self->name;
	}

	return Sequences::Genomic_Sequence->new_from_seq(
		$seq,    # Returns Bio::Seq object
		-five_prime_pos  => 0,    # TODO technically we could try to extract
		-three_prime_pos => $seq->length,    # this from the FASTA description
		-strand          => 1,
		-source => { db => $self, location => $location },
		-species => $self->name,
	);

}

=head2 Serialization override

=cut

sub to_hash {
	my $self = shift;
	return { 'fastaname' => $self->{fastaname}, };
}

=head2 Serialization override

=cut

sub serialization_restore {
	my $oldself = shift;

	my $self =
	  Sequences::Database::Sequence::FASTA->new( $oldself->{fastaname} );
	%$oldself = %$self;
}

1;
