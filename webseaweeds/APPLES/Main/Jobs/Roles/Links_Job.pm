#!/usr/bin/perl

=head1 A links job is a job which queries a links network dynamically

This type of job will store a Links::Links_Database data item and use
versioning to merge it before it accesses it.

=cut

package Jobs::Roles::Links_Job;

use feature ':5.10';

use Runtime;

use Moose::Role;
use MooseX::Declare;

use Datatypes::Moose_Types;
use Digest::SHA;
use JSON;

require Links::Links_Database;

## this also requires the Dynamic_Data_Loader Role
## we need a method get_dataitem from job to
requires('get_dataitem');
requires('load_data');

has 'nc__links_network' => (
	is        => 'ro',
 	isa     => 'Str',
	default => sub { return 'Links::Links_Database' },
	documentation => "_ automatically filled",
);

has 'ns__links_network' => (
is    => 'rw',
  isa => 'Links::Links_Database', );

has 'ns__checksum__nc__links_network' => (
	is => 'rw',
	isa => 'Str',
);

has "node_selection" => (
	is            => "rw",
	isa           => "ArrayRef[Str]",
	default       => sub { return [] },
	documentation => "_ automatically filled: selection of nodes",
);


=head2 List all Links nodes of a certain type

 Parameters:
 $self : a self object
 $type : a classname

 Returns:
 all nodes of this type that are contained in this Job's Links network
 in an ArrayRef

=cut

sub list_nodes_by_type {
	  my $self = shift;
	  my $type = shift || 'Links::Node';

	  if ( !defined( $self->ns__links_network ) ) {
		  confess("You need to load data first before calling this sub.");
	  }

	  my $nodes = $self->ns__links_network->get_nodes;
	  my @ret   = ();
	  foreach my $n (@$nodes) {
		  if ( UNIVERSAL::isa( $n, $type ) ) {
			  push @ret, $n;
		  }
	  }
	  return \@ret;
}

=head2 List all Links nodes that were selected

 Parameters:
 $self : a self object

 Returns:
 the corresponding nodes

=cut

sub list_selected_nodes {
	  my $self = shift;

	  my $ids  = $self->node_selection;

	  if ( !defined( $self->ns__links_network ) ) {
		  confess("You need to load data first before calling this sub.");
	  }

	  my $nodes = $self->ns__links_network->get_nodes;
	  my @ret   = ();

	  foreach my $n (@$nodes) {
		  if ( scalar( grep ( { $n->unique_id eq $_ } @$ids ) ) == 1 ) {
			  push @ret, $n;
		  }
	  }
	  return \@ret;
}

=head2 Return the links network

 Parameters:
 $self : a self object

 Returns:
 a Links_Database

=cut

sub links_network {
	  my $self = shift;
	  if ( !defined( $self->ns__links_network ) ) {
		  $self->load_data('nc__links_network');
	  }
	  return $self->ns__links_network;
}

=head2 Merge previous versions, and store

 Parameters:
 $self : a self object
 $network : a network to merge with

 Returns:
 Nothing

=cut

sub merge_links_network {
	  my $self    = shift;
	  my $network = shift;

	  if ( !UNIVERSAL::isa( $network, "Links::Links_Database" ) ) {
		  confess(
"incorrect parameter type for network, must be Links::Links_Database"
		  );
	  }

	  $self->put_dataitem( $network, $self->nc__links_network );
	  $self->merge_dataitem_versions( $self->nc__links_network,
		  \&Links::Links_Database::merge_nodes );
}

=head2 Replace previous versions, and store

 Parameters:
 $self : a self object
 $network : a network to merge with

 Returns:
 Nothing

=cut

sub replace_links_network {
	  my $self    = shift;
	  my $network = shift;

	  if ( !UNIVERSAL::isa( $network, "Links::Links_Database" ) ) {
		  confess(
"incorrect parameter type for network, must be Links::Links_Database"
		  );
	  }

	  $self->put_dataitem( $network, $self->nc__links_network );
	  $self->purge_dataitem_versions( $self->nc__links_network );
}


1;
