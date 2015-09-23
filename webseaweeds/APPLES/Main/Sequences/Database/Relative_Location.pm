#!/usr/bin/perl

=head1 Class Sequences::Database::Relative_Location

Stores locations relative to a sequence identifier/accession.

This is used as an input for the retrieval functions in
Sequence databases.

=cut

use MooseX::Declare;

class Sequences::Database::Relative_Location extends Datatypes::Datatype {

	use Scalar::Util qw(blessed);
	use Datatypes::Moose_Types;

=head2 Members

 identifiers : sequence identifiers (ARRAYREF)
 anchor : the start anchor
 offset : [optional] offset to use
 length : [optional] limit the length
 stop_at_neighbours   : [optional] decrease upstream_length or downstream_length
                      to stop at neigbouring genes. default: 0

The length of the retrieved sequence is determined as follows.

All coordinates are relative to the interval defined by the length of the
object corresponding to the identifier, plus upstream and downstream lengths.

   [...upstream area...]<...identifier object...>[...downstream area...]
   |                    |                       |                      |
  [0]     ...      [ 5' end ]       ...     [ 3' end ]        ...    [-1]

* length specifies the length of sequence to retrieve, starting at anchor + offset
* offset and length must not refer to a sequence that reaches outside the
  interval shown above.
* Stopping at previous gene might not be supported by all database types

=cut

	has "identifier" => (
						  is            => "rw",
						  isa           => "Str",
						  documentation => "1. Identifier/Accession",
						  default       => sub { return "" },
	);

	has "anchor" => (
					  is            => "rw",
					  isa           => "SequenceAnchor",
					  default       => sub { return "5' end"; },
					  documentation => "2. Anchor type",
	);

	has "offset" => (
					  is            => "rw",
					  isa           => "Int",
					  default       => sub { return 0; },
					  documentation => "3. Offset",
	);

	has "length" => (
					  is            => "rw",
					  isa           => "Int",
					  default       => sub { return 1; },
					  documentation => "4. Length",
	);

	has "stop_at_neighbours" => (
								 is            => "rw",
								 isa           => "Bool",
								 default       => sub { return 0; },
								 documentation => "5. Stop at neigbouring gene",
	);

=head2 Helper to produce the same location, but on the opposite strand

 We switch the anchor, and multiply the offset and length with -1

=cut

	method mirror_strand () {
		$self->offset( $self->offset * -1 );
		$self->length( $self->length * -1 );
		if ( $self->anchor eq "5' end" ) {
			$self->anchor("3' end");
		} else {
			$self->anchor("5' end");
		}
	}

=head2 Return true if this is a 'simple' location

Simple locations correspond to what a sequence database will return
when only using the sequence identifier/accession.

 Parameters:
 $self : a Sequences::Database::Relative_Location

=cut

	method is_simple () {
		return
		     ( $self->anchor eq "5' end" )
		  && ( $self->offset == 0 )
		  && ( $self->length == -1 )
		  && ( $self->stop_at_neighbours == 0 );
	}

=head2 Unique ID override

=cut

	method unique_id () {
		return
		    blessed($self)
		  . ": $self->{identifier},$self->{anchor}"
		  . "$self->{offset},$self->{length},$self->{stop_at_neighbours}";
	}

}
