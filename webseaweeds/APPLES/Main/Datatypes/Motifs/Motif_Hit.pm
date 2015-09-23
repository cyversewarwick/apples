#!/usr/bin/perl

use MooseX::Declare;

=head1 PSSM Motif hit

Stores a PSSM and related information.

=cut

class Datatypes::Motifs::Motif_Hit extends Links::Node {

	## Five prime position
	has 'five_prime_pos' => (
		is => "rw",
		isa => "Int",
        required => 1,
	);

	## Three prime position
	has 'three_prime_pos' => (
		is => "rw",
		isa => "Int",
        required => 1,
	);

	## p value
	has 'pvalue' => (
		is => "rw",
		isa => "Num",
        required => 1,
	);

	## p value
	has 'qvalue' => (
		is => "rw",
		isa => "Num|Undef",
        required => 0,
	);

	## Score
	has 'score' => (
		is => "rw",
		isa => "Num",
        required => 1,
	);

	## PSSM name
	has 'name' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

	## The matched sequence
	has 'sequence' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

	## The strand it was matched on.
	has 'strand' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

}
