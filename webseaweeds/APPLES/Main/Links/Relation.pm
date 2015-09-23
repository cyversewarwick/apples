#!/usr/bin/perl

=head1 Link/asymmetric relation class to connect nodes

=cut

use MooseX::Declare;

class Links::Relation extends Serialization::Non_Serializable {

	use Scalar::Util qw(blessed);

	use Runtime;
	use Carp;

	require Links::Node;

## the source node
	has 'source' => (
					  is       => "rw",
					  isa      => "Links::Node",
					  required => 1,
#					  weak_ref => 1,
	);

## the target node
	has 'target' => (
					  is       => "rw",
					  isa      => "Links::Node",
					  required => 1,
#					  weak_ref => 1,
	);

## the data node
	has '_data' => (
					 is       => "rw",
					 isa      => "HashRef[Num|Str|Links::Node]",
					 default  => sub { return {} },
					 weak_ref => 1,
	);

## mark node as bidirectional
	has 'bidirectional' => (
							 is      => "ro",
							 isa     => "Int",
							 default => sub { return 0; },
	);

=head2 Update from another relation

This function updates this relation's attributes from another relation's
attributes.

Attributes already contained will not be changed.

 Parameters:
 $other : another relation 
 
 Returns:
 nothing

=cut

	method update(Links::Relation $other) {
		foreach my $attr ( keys %{ $other->{_data} } ) {
			if(!defined ($self->{_data}->{$attr})) {
				$self->{_data}->{$attr} = $other->{_data}->{$attr};
			}
		}
	}

=head2 Get a relation name/label

 Parameters:
 $self : a Links::Relation

 Returns:
 A name/label for the relation
 
=cut

	method label () {
		my $result = $self->{source}->label . " ->  " . $self->{target}->label;

		if ( defined( $self->{_data} ) ) {
			$result .= " (";

			while ( my ( $k, $v ) = each( %{ $self->{_data} } ) ) {
				$result .= "$k : $v ;";
			}
			$result .= ")";
		}

		return $result;
	}

=head2 Get the value associated with this node

 Parameters:
 $self : a Links::Relation
 $item : the attribute of the relation
 $value : (optional) the value to set for this attribute

 Returns :
 the data associated with this node or undef
 

=cut

	method data (Str $item, Any $value?) {
		if ( defined($value) ) {
			$self->{_data}->{$item} = $value;
			$self->meta->get_attribute('_data')
			  ->verify_against_type_constraint( $self->{_data} );
		}

		return $self->{_data}->{$item};
	}

=head2 Get a list of attribute names

 Returns: 
 an array of attributes

=cut

	method attributes () {
		if ( defined( $self->{_data} ) ) {
			return keys %{ $self->{_data} };
		} else {
			return ();
		}
	}
}
