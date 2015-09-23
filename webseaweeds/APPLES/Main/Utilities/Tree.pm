#!/usr/bin/perl 
use MooseX::Declare;

=head1 Tree class

 Contains functions for adding children and traversing top-down

=cut

class Utilities::Tree extends Serialization::Serializable {

	use Runtime;
	use Carp;
	use Data::Dumper;
	use Scalar::Util qw(blessed);

	# Parent node of this tree
	has 'parent' => (
		is      => 'rw',
		isa     => 'Undef|Utilities::Tree',
		default => undef,
	);

	# Array of children
	has 'children' => (
		is      => 'rw',
		isa     => 'ArrayRef[Utilities::Tree]',
		default => sub { [] },
	);

	# Each tree node can store a value
	has 'value' => (
		is      => 'rw',
		isa     => 'Any',
		default => undef,
	);

=head2 add a child node

 Parameters:
 $value : the value to add to the node

 Returns:
 The node that was created

=cut

	method add_child( Any $value ) {
		my $class = blessed($self);
		  my $n   = $class->new();
		  $n->value($value);
		  $n->parent($self);
		  push @{ $self->{children} }, $n;
		  return $n;
	  }

=head2 add a tree as a child node

 Parameters:
 $n : a subtree to append

 Returns: $self

=cut

	  method add_child_node( Utilities::Tree $n) {
		$n->parent($self);
		  push @{ $self->{children} }, $n;
		  return $self;
	  }

=head2 Traverse tree top to bottom

 Parameters:
 $c : a coderef that might return something for each node

 Returns:
 An array containing all values $c->($v) for all values $v in the tree

=cut

	  method traverse( CodeRef $c) {
		my $arr = [];
		{ # Moose Trivia: mess with $_ and Moose f**cks up $self. Just like that.
			local $_ = [ $self->value() ];
			push @$arr, $c->();
		}

		foreach ( @{ $self->{children} } ) {
			my $arr2 = $_->traverse($c);

			foreach (@$arr2) {
				push @$arr, $_;
			}
		}
		return $arr;
	  }

=head2 Traverse tree top to bottom, passing the value by reference

 Parameters:
 $c : a coderef that might return something for each node

 Returns:
 An array containing all values $c->($v) for all values $v in the tree

=cut

	  method traverse_update( CodeRef $c) {
		{ # Moose Trivia: mess with $_ and Moose f**cks up $self. Just like that.
			local $_ = [ $self->value() ];
			  $self->{value} = $c->();
		}

		foreach ( @{ $self->{children} } ) {
			$_->traverse_update($c);
		}
	  }

=head2 Convert tree to JSON-compatible tree.

 Returns:
 untyped, JSON-compatible tree (removing the parent links)
=cut

	  method to_hash() {
		my %h = %$self;
		  $h{parent}   = undef;
		  $h{children} = [];
		  foreach my $c ( @{ $self->{children} } ) {
			push @{ $h{children} },
			  Serialization::Serializable::to_hash( $c->TO_JSON() );
		}

		return \%h;
	  }

=head2 restore parent links and blessings after serializing to JSON

 Parameters:
 $hash : JSON-compatible hash (see TO_JSON)

 Returns:
 $self
=cut

	  method serialization_restore() {
		foreach my $c ( @{ $self->{children} } ){
			$c->{parent} = $self;
			  $c->serialization_restore;
		  } return $self;
	  };
};
