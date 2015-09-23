#!/usr/bin/perl

=head1 Class Links::Links_Database

=cut

package Links::Links_Database;

use strict;

use Runtime;
use Carp;

our @ISA = qw(Serialization::Serializable);
require Serialization::Serializable;

=head2 Constructor

 Parameters: 
 $class : Links::Links_Database
 
 Returns 
 a new Links::Links_Database object 

=cut

sub new {
	my $class = shift;

	my $self = bless { 'nodelist' => {}, }, $class;

	return $self;
}

=head2 Add a node

 All nodes will be added to this database. 

 Parameters:
 $self : self object
 @nodelist : a list of Links::Node's 
 
 Returns:
 nothing

=cut

sub add_node {
	my ( $self, @nodelist ) = @_;

	my %nodes_checked = ();

	my %nodes_to_replace = ();

	while ( scalar @nodelist > 0 ) {
		my $n = shift @nodelist;
		$nodes_checked{ $n->unique_id } = $n;
		foreach my $othernode ( @{ $n->get_neighbours } ) {
			unless ( defined( $nodes_checked{ $othernode->unique_id } ) ) {
				push @nodelist, $othernode;
			}
		}

		my $old_n = $self->{nodelist}->{ $n->unique_id };

		## If the node is replacing a node with the same ID,
		## we need to copy all the links, too
		if ( defined($old_n) && $old_n != $n ) {
			confess(
"Cannot add node with the same ID which is already contained in the network."
			);
		}
		$self->{nodelist}->{ $n->unique_id } = $n;
	}

	while ( my ( $k, $node ) = ( each %nodes_to_replace ) ) {
	}
}

=head2 Get a node

 Returns:
 The node with the same unique_id in the network, 
 or undef

=cut

sub get_or_add_node {
	my $self = shift;
	my $node = shift;
	
	if(!defined($self->{nodelist}->{ $node->unique_id })) {
		$self->add_node ($node);
	}
	
	return $self->{nodelist}->{ $node->unique_id };
}

=head2 Get all nodes

 Parameters:
 $self : self object
 
 Returns:
 All nodes in the network

=cut

sub get_nodes {
	my $self = shift;
	if (!Serialization::Serializable::is_hash($self) ||
		!UNIVERSAL::isa( $self, 'Links::Links_Database' )
	  )
	{
		confess(" \$self must be a Links::Links_Database, is '$self'");
	}
	my @list = values( %{ $self->{nodelist} } );
	return \@list;
}

=head2 Get all links between nodes

 Parameters:
 $self : self object
 
 Returns:
 All relation links between nodes in the network

=cut

sub get_links {
	my $self = shift;

	my $list = $self->get_nodes;

	my @rels = ();

	my %bidir = ();
    
    foreach (@$list) {
        foreach my $l ( @{ $_->get_links } ) {
            if($l->bidirectional && !exists($bidir{ $l->source->unique_id . $l->target->unique_id }))
            {
                #print "\nBidir link from " . $l->source->accession . " to " . $l->target->accession . " (" . $l->bidirectional . ")";
                
				push @rels, $l;
				$bidir{ $l->source->unique_id . $l->target->unique_id } = 1;
            }
            
        }
    }

	foreach (@$list) {
		foreach my $l ( @{ $_->get_links } ) {
            #print "\nLink source " . $l->source->accession . " to " . $l->target->accession . " (" . $l->bidirectional . ")";

			if (
				!$l->bidirectional
				|| !exists(
					$bidir{ $l->source->unique_id . $l->target->unique_id }
				)
			  )
			{
                #print "\nAccepted";
				push @rels, $l;

				$bidir{ $l->source->unique_id . $l->target->unique_id } = 1;
			}
            else
            {
                #print "\nDenied";
            }
		}
	}

	return \@rels;
}

=head2 Merge two networks, merging nodes which are equal

For large networks, this may take O(n^2) time, where n is
the number of nodes. 

 Parameters:
 $net1 : the first network
 $net2 : the second network 

 Returns: 
 a network with all the nodes
 
=cut

sub merge_nodes {
	my $net1 = shift;
	my $net2 = shift;

	my %ser_nodes = ();

	foreach my $n ( @{ $net1->get_nodes } ) {
		$ser_nodes{ $n->unique_id } = Serialization::Serializable::to_hash($n);
	}

	foreach my $n ( @{ $net2->get_nodes } ) {
		my $nser = Serialization::Serializable::to_hash($n);
		if ( defined( $ser_nodes{ $n->unique_id } ) ) {
			## 'merge' the nodes -- update all values, keep the old ones
			while ( my ( $k, $v ) = each(%$nser) ) {
				$ser_nodes{ $n->unique_id }->{$k} = $v;
			}
		} else {
			$ser_nodes{ $n->unique_id } = $nser;
		}
	}

	my %ser_links = ();

	foreach my $l ( @{ $net1->get_links } ) {
		my $lser = $net1->_serialize_link($l);
		my $lsid =
		  $lser->{source} . "->" . $lser->{target} . " " . $lser->{bidir};
		$ser_links{$lsid} = $lser;
	}

	foreach my $l ( @{ $net2->get_links } ) {
		my $lser = $net2->_serialize_link($l);
		my $lsid =
		  $lser->{source} . "->" . $lser->{target} . " " . $lser->{bidir};

		if ( defined( $ser_links{$lsid} ) ) {

		} else {
			$ser_links{$lsid} = $lser;
		}
	}

	my @rnodes = values %ser_nodes;
	my @rlinks = values %ser_links;

	my $result = bless {
		'nodes' => \@rnodes,
		'links' => \@rlinks,
	  },
	  __PACKAGE__;

	$result->serialization_restore;

	return $result;
}

=head2 Make a serializable version of a link

 Parameters:
 $self : a self object
 $l    : a Links::Relation
 
 Returns: 
 a hash with the link data, relative to $self.
 
=cut

sub _serialize_link {
	my $self = shift;
	my $l    = shift;

	my $ln = {
		source => $l->source->unique_id,
		target => $l->target->unique_id,
		bidir  => $l->bidirectional,
	};

	my $data = $l->_data;
	if ( defined($data) ) {
		my %d  = %{$data};
		my $rd = \%d;
		bless $rd, "Serialization::Serializable";
		$ln->{node} = Serialization::Serializable::to_hash($rd);
	}
	return $ln;
}

=head2 Deserialize a link that has been serialized using _serialize_link
       into the current network

 Parameters:
 $self : a self object
 $l    : a link hash, relative to self;
 
 Returns: 
 The created Links::Relation object
=cut

sub _deserialize_link {
	my $self = shift;
	my $l    = shift;

	my $s = $self->{nodelist}->{ $l->{source} };
	my $t = $self->{nodelist}->{ $l->{target} };

	if ( !defined($s) || !defined($t) ) {
		confess("Inconsistent link with nonexistent nodes $s / $t");
	}

	my $rel;

	if ( defined( $l->{bidir} ) && $l->{bidir} ) {
		$rel = $s->link_bidirectional($t);
	} else {
		$rel = $s->link_to($t);
	}

	my $data = $l->{node};
	if ( defined($data) ) {
		$data = Serialization::Serializable::from_hash( $l->{node} );
		while ( my ( $k, $v ) = each %$data ) {
			next if $k eq "SERIAL_VERSIONID";
			$rel->data( $k, $v );
		}
	}
	return $rel;
}

=head2 To HASH for Serialization::Serializable

 Parameters:
 $self : self object
 
 Returns:
 a hash that can be serialized to JSON/XML

=cut

sub to_hash {
	my $self = shift;

	my $nodelist = $self->get_nodes;

	my $result = {
		nodes => [],
		links => [],
	};

	foreach my $n (@$nodelist) {
		push @{ $result->{nodes} }, Serialization::Serializable::to_hash($n);
	}

	my $linklist = $self->get_links;
	foreach my $l (@$linklist) {
		my $ln = $self->_serialize_link($l);
		push @{ $result->{links} }, $ln;
	}

	return $result;
}

=head2 Restore from a hash

 Parameters:
 $self : self object
 
 Returns:
 nothing

=cut

sub serialization_restore {
	my $self = shift;

	$self->{nodelist} = {};

	foreach my $n ( @{ $self->{nodes} } ) {
		my $nd = Serialization::Serializable::from_hash($n);
		$self->{nodelist}->{ $nd->unique_id } = $nd;
	}

	delete $self->{nodes};

	foreach my $l ( @{ $self->{links} } ) {
		$self->_deserialize_link($l);
	}
}

1;
