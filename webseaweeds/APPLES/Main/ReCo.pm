### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReCo: Regulatory Complex class ###
# set of molecules (such as proteins) which can have likelihoods
# can be used to summarise potential regulators present on a gene's ReRe

use MooseX::Declare;

class ReCo {
    use Parameters;
    use Molecule;
    
    has 'molecules' => (is => 'rw', isa => 'ArrayRef[Molecule]');

    method add_molecule(Molecule $molecule) {
	push(@{$self->molecules},$molecule);
    } # add_molecule #

    method lower_bound_on_likelihood(Str $unique_ID, Num $sampling_likelihood) { 
	# for molecule(s) matching the ID, sets likelihood to given value unless existing likelihood is already greater

	my $number_of_molecules = @{$self->molecules};
	for (my $i=0;$i<$number_of_molecules;$i++) {
	    my $molecule = ${$self->molecules}[$i];
	    my $molecule_ID = $molecule->unique_ID();
	    if ($molecule_ID eq $unique_ID) {
		my $current_likelihood = -1;
		if (defined $molecule->likelihood) {
		    $current_likelihood = $molecule->likelihood;
		}
		if ($sampling_likelihood > $current_likelihood) {
		    ${$self->molecules}[$i]->likelihood($sampling_likelihood);
		}
	    }
	}
    } # lower_bound_on_likelihood #

} # ReCo #
