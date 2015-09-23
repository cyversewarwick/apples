#!/usr/bin/perl

=head1 Serializable Array Helper

Helper class to serialize arrays

=cut

package Serialization::Serializable_Array;

use strict;

use Runtime;
use Carp;

our @ISA = qw(Serialization::Serializable);
require Serialization::Serializable;

=head2 Constructor

 Parameters: 
 $class : Serialization::Serializable_Array
 
 Returns 
 a new Serialization::Serializable_Array object

=cut

sub new {
	my ( $class, @data ) = @_;

	my $self = bless { data => \@data, }, $class;

	return $self;
}

=head2 Getter/Setter for data member

 Parameters:
 $self : $self object
 $val : new value (optional)

=cut

sub data {
	my $self = shift;
	my $val  = shift;

	if ( defined($val) ) {
		$self->{data} = $val;
	}
	return $self->{data};
}

1;
