#!/usr/bin/perl

use MooseX::Declare;

=head1 Link node class

This class abstracts the concept of a linkable node for representing
relational links between data. Nodes can be connected using
the Relation class.

=cut

class Links::Node extends Datatypes::Datatype {

	use Scalar::Util qw(blessed);
	use Data::UUID;
	require Links::Relation;

	## every node has a unique id
	has 'unique_id' => (
		is      => "ro",
		isa     => "Str",
		default => sub {
			my $ug = Data::UUID->new();
			return $ug->create_str();
		},
	);

	## this is a list of links to other nodes
	has 'links' => (
		is      => "ro",
		isa     => "HashRef[Links::Relation]",
		default => sub { return {} },
	);

	## Node labels
	has 'label' => (
		is  => "ro",
		isa => "Str",
	);

	## a dot language style for pretty-printing
	has 'dot_style' => (
		is  => 'rw',
		isa => 'Str',
	);

	## nodes can be associated with a linkout URL
	has 'url' => (
		is => 'rw',
		isa => 'Str',
		builder => '_build_url',
	);

	## override this to change the color in derived classes
	method color () {
		return "#eee";
	}

=head2 Get a name/label

Override this method for custom labels.

 Parameters:
 $self : self object
 Returns:
 A name/label for the relation

=cut

	method label () {
		return $self->{label} || $self->unique_id;
	}

=head2 Add a link to another node

 Parameters:
 $self : self object
 $node : the node to link to

 Returns:
 a (new, if didn't exist) Links::Relation of the link to the other node

=cut

	method link_to (Links::Node $node) {
		my $rel = $self->{links}->{ $node->unique_id };

		if ( !defined($rel) ) {
			$rel = Links::Relation->new( source => $self, target => $node );
			$self->{links}->{ $node->unique_id } = $rel;

		}
		return $rel;
	}

=head2 Add a bidirectional link to another node

 Parameters:
 $self : self object
 $node : the node to link to and from

 Returns:
 a (new, if didn't exist) Links::Relation of the link to the other node

=cut

	method link_bidirectional (Links::Node $node) {
		my $rel = $self->{links}->{ "bidir_" . $node->unique_id }
		  || $node->{links}->{ "bidir_" . $self->unique_id };

		if ( !defined($rel) ) {
			$rel = Links::Relation->new(
				source        => $self,
				target        => $node,
				bidirectional => 1,
			);
			$self->{links}->{ "bidir_" . $node->unique_id } = $rel;
			$node->{links}->{ "bidir_" . $self->unique_id } = $rel;
		}
		return $rel;
	}

=head2 Remove a link

 Parameters:
 $self : self object
 $node : the node to remove the (directed) link to

 Returns:
 nothing

=cut

	method remove_link (Links::Node $node) {
		my $rel = $self->{links}->{ $node->unique_id };
		if ( defined($rel) ) {
			delete $self->{links}->{ $node->unique_id };
		}
	}

=head2 Remove a bidirectional link

 Parameters:
 $self : self object
 $node : the node to remove the bidirectional link to

 Returns:
 nothing

=cut

	method remove_bidirectional_link (Links::Node $node) {
		my $rel = $self->{links}->{ "bidir_" . $node->unique_id };
		if ( defined($rel) ) {
			delete $node->{links}->{ "bidir_" . $self->unique_id };
			delete $self->{links}->{ "bidir_" . $node->unique_id };
		}
	}

=head2 Check if we have a link to a given node

 Parameters:
 $self : self object
 $node : the node to check the link to

 Returns:
 The link to the node, or undef.

 If the node has a directed and a bidirectional link,
 the directed link is returned.

=cut

	method has_link (Links::Node $node) {
		return $self->{links}->{ $node->unique_id }
		  || $self->{links}->{ "bidir_" . $node->unique_id };
	}

=head2 Check if we have a bidirectional link to a given node

 Parameters:
 $self : self object
 $node : the node to check the link to

 Returns:
 The bidirectional link to the node, or undef.

=cut

	method has_bidirectional_link (Links::Node $node) {
		return $self->{links}->{ "bidir_" . $node->unique_id };
	}

=head2 Get the all neighbours

 Parameters:
 $self : self object
 $relationtype : regular expression to match
 				   against the relation type
 				   (or undef)
 $nodetype : regular expression to match
             the node type to return
             (or undef)

 Returns:
 an ARRAYREF with all neighbouring nodes (of given type)

=cut

	method get_neighbours (Str $relationtype = ".*", Str $nodetype = ".*") {
		my $n = [];
		while ( my ( $k, $v ) = each( %{ $self->{links} } ) ) {
			my $t = $v->target;
			
			if (!defined ($t) || $t == $self) {
				$t = $v->source;
			}
			
			if (!defined ($t)) {
				my $s = $v->source;
				my $t = $v->target;
				die "Invalid link : $v : $s -> $t";
			}
			
			push @$n, $t;
#			  if blessed($v) =~ m/$relationtype/ && blessed($t) =~ m/$nodetype/;
		}
		return $n;
	}

=head2 Get the out degree of this node

 Parameters:
 $self : self object

 Returns:
 The out degree

=cut
	method out_degree () {
		return scalar (keys %{$self->{links}});
	}

=head2 Get all the links to neighbours

 Parameters:
 $self : self object

 Returns:
 All links from this node in an ARRAYREF

=cut

	method get_links () {
		my @a = values %{ $self->{links} };
		return \@a;
	}

=head2 Reset all links (including bidirectional ones)

 Parameters:
 $self : self object

 Returns:
 Nothing.

=cut

	method reset_links () {
		while ( my ( $k, $v ) = each( %{ $self->{links} } ) ) {
			if ( $v->bidirectional ) {
				$self->remove_bidirectional_link( $v->target );
			}
		}
		$self->{links} = {};
	}

=head2 Restore from a hash

 We don't serialize links, this is done in Links_Network

 Parameters:
 $self : self object

 Returns:
 nothing

=cut

	method serialization_restore () {
		$self->{links} = {};
	}

=head2 URL element builder

Override this in subclasses to create the URL element.

=cut

	method _build_url () {
		return "";
	}
}
