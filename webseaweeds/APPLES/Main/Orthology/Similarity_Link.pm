#!/usr/bin/perl

=head1 Similarity Link data 

=cut

use MooseX::Declare;

class Orthology::Similarity_Link extends Links::Node {

	use Runtime;
	
	## this is BLAST hit data
	has 'hit' => (
		is => "rw",
		isa => "HashRef",
		required => 1,
	);
	
	## the dot style
	method dot_style () {
		my $wstr = sprintf "%g",  $self->weight ;
		return "color=grey shape=hexagon label=\"$wstr\"";
	}
	
	## this returns the hit significance as the edge weight
	method weight () {
		return $self->{hit}->{significance} || 100;
	}

}
