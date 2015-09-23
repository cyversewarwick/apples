#!/usr/bin/perl

use MooseX::Declare;

=head1 Everything that comes out from comparing two sequences using an alignment plot

=cut

class Datatypes::Results::Alignment_Plot extends Datatypes::Sequence::Set {

	## the job that created this plot
	has 'job' => (
		is => 'rw',
		isa => 'Jobs::Subtasks::Seaweed_Job',
		required =>  1,
	);

	## a histogram of scores
	has 'histogram' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
		required =>  1,
	);

	## conservation profile for sequence 1
	has 'profile_1' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
		required =>  1,
	);

	## conservation profile for sequence 2
	has 'profile_2' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
		required =>  1,
	);

	## all the scores above threshold in a plot
	has 'plot' => (
		is => 'rw',
		isa => 'HashRef[Num]',
		required =>  1,
	);

	## first sequence
	has 'sequence_1' => (
		is => 'rw',
		isa => 'Sequences::Genomic_Sequence',
		required =>  1,
	);

	## second sequence
	has 'sequence_2' => (
		is => 'rw',
		isa => 'Sequences::Genomic_Sequence',
		required =>  1,
	);

=head2 Override and do nothing.

=cut
	method compact_sequences() {
		
	}

=head2 Overloaded sequences method

Datatypes::Sequence::Set has an element "sequences", we only return
the two sequences here.

=cut
	method sequences () {
		return [
			$self->sequence_1(),
			$self->sequence_2(),
		];
	}

};