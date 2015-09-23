#!/usr/bin/perl

use MooseX::Declare;

=head1 Compute reciprocal best BLAST hits

This job returns a network of BLAST scores, together with reciprocal
best hit information.

=cut

class Jobs::Subtasks::RBH extends Jobs::Job {
	use Runtime;
	use Orthology::Search::RBH;

	use Datatypes::Moose_Types;
	use Links::Links_Database;

	has 'startID' => (
		is  => 'rw',
		isa => 'Datatypes::Sequence::Location',
	);

	has "targetDB" => (
		is       => "rw",
		isa      => "Search_Database_ID",
		required => 1,
	);

	has "max_hits" => (
		is      => "rw",
		isa     => "Int",
		default => 1,
	);

=head2 Overridden run method

=cut

	method _run () {
		my $rbhs =
		  Orthology::Search::RBH->new( $self->startID->db, $self->targetDB, );

		my $cn = $rbhs->search( $self->startID->accession, $self->max_hits );

		return $cn;
	}
}
