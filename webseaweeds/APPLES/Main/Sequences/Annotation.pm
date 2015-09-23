#!/usr/bin/perl

package Sequences::Annotation;

=head1 Sequence annotation class

A little counterintuitively, this is not the standard thing associated with sequence
annotation, rather, it is just a marker associated with a set of genomic coordinates
on a Genomic_Sequence object.

=cut

use strict;

use Serialization::Serializable;
our @ISA = qw(Serialization::Serializable);

use Hash::Merge;

=head2 Constructor

Construct a new annotation.

 Parameters:
 $class :  Sequences::Annotation
 $id  :  a unique identifier
 $type :  a type of annotation
 $description : a textual description
 $sourceid : an id of the source of the annotation
 $five_prime_pos :  five prime end 
 $three_prime_pos :  three prime end
 $strand : strand (-1 or 1)
 $data : some data to associate with the annotation

 Returns: 
 A new Annotation
=cut

sub new {
	my (
		$class,           $id,       $type,
		$description,     $sourceid, $five_prime_pos,
		$three_prime_pos, $strand,   $data
	) = @_;


	if ( $strand =~ m/([+\-])?[0-9]+/) {
		if ( $strand < 0 ) {
			$strand = 'negative';
		} else {
			$strand = 'positive';
		}		
	} elsif ($strand =~ m/^n/i) {
		$strand = 'negative';			
	} elsif ($strand =~ m/^p/i) {
		$strand = 'positive';			
	} else {
		confess ("Invalid strand information: $strand");
	}

	my $ann = {
		id              => $id,
		type            => $type,
		description     => $description,
		sourceid        => $sourceid,
		five_prime_pos  => $five_prime_pos,
		three_prime_pos => $three_prime_pos,
		strand          => $strand,
		ns__data        => $data,
	};

	bless $ann, $class;
	return $ann;
}

=head2 Return an identifier which says uniquely what this is, and where it is.

This is used for merging, two sequence annotations with the same ident will be 
merged into one.

 Parameters:
 $self : a self object

 Returns:
 a unique identifier string for this ann.

=cut

sub ident {
	my $self = shift;
	return "$self->{type} : $self->{sourceid} ($self->{description}), $self->{five_prime_pos} - $self->{three_prime_pos} / $self->{strand}";
}

=head2 From two make one, store in $self

 This only works if the two ident values are the same (otherwise, there is no )

 Parameters:
 $self : a self object
 $rhs : another Annotation object

=cut

sub merge {
	my $self = shift;
	my $rhs = shift;
	my $merge = Hash::Merge->new('LEFT_PRECEDENT');
	
	my $id = $self->{id};
	my $type= $self->{type};
	my $sourceid = $self->{sourceid};
	my $five_prime_pos = $self->{five_prime_pos};
	my $three_prime_pos = $self->{three_prime_pos};
	my $strand = $self->{strand};
	
	$self = $merge->merge ($self, $rhs);

	$self->{id} = $id;
	$self->{type} = $type;
	$self->{sourceid} = $sourceid;
	$self->{five_prime_pos} = $five_prime_pos;
	$self->{three_prime_pos} = $three_prime_pos;
	$self->{strand} = $strand;	
}

1;
