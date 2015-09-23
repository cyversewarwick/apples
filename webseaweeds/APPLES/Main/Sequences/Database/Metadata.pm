#!/usr/bin/perl

=head1 Class Sequence_Database;

Generic sequence metadata database interface.



=cut

package Sequences::Database::Metadata;

use feature ':5.10';
use strict;

use Moose::Role;
use MooseX::Declare;

use Runtime;
use Carp;

use Serialization::Serializable;

requires('_get_meta_data');

use constant {
	TYPE_GENE => 1,
	TYPE_TRANSCRIPT => 2,
	TYPE_EXON => 3,
	TYPE_OTHER => 4,
};

=head2 Get all data available for an identifier

 Parameters:
 $self : a self object
 $id   : an id/accession

 Returns:
 {
 	url => an url (if exists),
 	common_name => common_name
 	type =>
 	...

 }

=cut

sub get_meta_data {
	my $self = shift;
	my $id   = shift;

	my $cid = bless {
		type => "SEQUENCE_METADATA",
		db => $self,
		id => $id,
	}, "Serialization::Serializable";

	my $val = cache->cache_get ($cid);

	if (!defined ($val)) {
		$val = $self->_get_meta_data ($id);

		cache->cache_put ($cid, $val);
	}

	return $val;
}


1;