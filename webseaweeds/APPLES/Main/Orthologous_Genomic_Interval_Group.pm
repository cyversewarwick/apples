### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Orthologous_Genomic_Interval_Group class ###

# this class is not implemented yet, just illustrating an idea
# rename to Homologous_Genomic_Interval_Group?

use MooseX::Declare;

class Orthologous_Genomic_Interval_Group {
	use Parameters;
	use APPLES_Datatypes qw (HomologyType);
	
	has 'genomic_interval' => (is => 'rw', isa => 'Genomic_Interval');
	has 'orthologous_genomic_interval_list' => (is => 'rw', isa => 'ArrayRef[Genomic_Interval]');
	has 'homology_type' => (is => 'rw', isa => HomologyType);
											
	method find_orthologous_genomic_interval() {
	
		die 'method not implemented yet!';
	
	} # find_orthologous_genomic_interval #

} # Orthologous_Genomic_Interval_Group #
