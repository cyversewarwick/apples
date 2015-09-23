### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Molecule Class ###

use MooseX::Declare;

class Molecule {
    use APPLES_Datatypes qw (Probability);

    has 'likelihood' => (is => 'rw', isa => Probability);

    method unique_ID() {
	# return type is Str

	die 'method not implemented yet.'
    } # unique_ID #

} # Molecule #

class Protein extends Molecule {
    has 'uniprot_id' => (is => 'ro', isa => 'Str');

    override unique_ID() {
	my $result;
	if (defined $self->uniprot_id) {
	    $result = $self->uniprot_id;
	} else {
	    die 'no uniprot-ID was provided to this object - cannot find a unique ID.';
	}
	return $result;
    } # unique_ID #

} # Protein #
