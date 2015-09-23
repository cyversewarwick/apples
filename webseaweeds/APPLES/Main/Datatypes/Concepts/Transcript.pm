#!/usr/bin/perl

use MooseX::Declare;

=head1 Sequence Location Class

Stores sequence locations and databases.

=cut

class Datatypes::Concepts::Transcript extends Datatypes::Sequence::Location {

	use Runtime;

	require Sequences::Database::Relative_Location;

	## force metadata retrieval
	has 'metadata' => (
		is => "rw",
		isa => "Any",
		builder => "get_meta_data",
		lazy => 1,
	);

=head2 Make a new Transcript

=cut

	sub new_transcript {
		my $class = shift;
		my $db = shift;
		my $acc = shift;
		my $v = Datatypes::Concepts::Transcript->new(
			db => $db,
			location =>
			  Sequences::Database::Relative_Location->new( identifier => $acc )
		);
		$v->metadata();
		return $v;	
	}

	## override this to change the color in derived classes
	method color () {
		return "#fcc";
	}

=head2 Create a dot label

=cut

	method dot_style () {
		return
		    "label=\"Transcript: "
		  . $self->accession . " ("
		  . $self->db . ")"
		  . "\" shape=\"ellipse\"";
	}

}
