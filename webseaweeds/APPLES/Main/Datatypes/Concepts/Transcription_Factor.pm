#!/usr/bin/perl

use MooseX::Declare;

=head1 Transcription Factor class

Represents a transcription factor and stores related information.

=cut

class Datatypes::Concepts::Transcription_Factor extends Links::Node {

    has 'id' => (
      is => 'rw',
      isa => 'Str',
      default => sub {return 'no id'}
    );
    
    has 'name' => (
      is => 'rw',
      isa => 'Str',
      default => sub {return 'no name'}
    );

    has 'species' => (
      is => 'rw',
      isa => 'Str',
      default => sub {return 'no species'}
    );

=head2 Create a dot label

=cut

	method dot_style () {
		return
		    "label=\"TF: "
		  . $self->name . " ("
		  . $self->species . ")"
		  . "\" shape=\"triangle\"";
	}

=head2 Unique ID override

=cut

	method unique_id () {
		return blessed($self) . " : $self->{id}";
	}


}