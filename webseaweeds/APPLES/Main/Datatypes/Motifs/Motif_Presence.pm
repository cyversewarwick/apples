#!/usr/bin/perl

use MooseX::Declare;

=head1 PSSM Motif Presence (detected by PSSMTool)

Stores a PSSM and related information.

=cut

class Datatypes::Motifs::Motif_Presence extends Links::Node {
	use Datatypes::Moose_Types;
	
	## P value from statistical test
	has 'pvalue' => (
		is => "rw",
		isa => "Num",
        required => 1,
	);

	## information about the test that was performed to create this data item
	has 'test' => (
		is => "rw",
		isa => "Overrep_Test",
        required => 1,
	);
	
	## Number of scores that were collected for this test
	has 'collected_scores' => (
		is => "rw",
		isa => "Int",
        required => 1,
	);

	## Number of scores that were deemed significant
	has 'count' => (
		is => "rw",
		isa => "Int",
        required => 1,
	);

	## name of the pssm that was matched
	has 'pssm_name' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

	## bifa/transfac accession of the matched pssm
	has 'pssm_accession' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

	## information about the sequence used 
	has 'sequence_name' => (
		is => "rw",
		isa => "Str",
        required => 1,
	);

}
