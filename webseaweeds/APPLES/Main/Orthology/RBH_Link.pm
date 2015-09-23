#!/usr/bin/perl

use MooseX::Declare;

=head1 Reciprocal best hit link 

=cut

class Orthology::RBH_Link extends Links::Node {

	## a RBH significance (normally the geometric mean of 
	## the individual scores)
	has 'significance' => (
		is => "rw",
		isa => 'Num',
		required => 1,
	);

	## the dot style
	method dot_style () {
		my $wstr = sprintf "%g", $self->weight;
		return "shape=diamond label=\"" . $self->weight . "\"";
	}
	
	## this returns the hit significance as the edge weight
	method weight () {
		return $self->significance;
	}
	
}