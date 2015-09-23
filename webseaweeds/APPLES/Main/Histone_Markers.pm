### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Histone Markers class ###
use MooseX::Declare;

# this class is not implemented yet, just illustrating an idea

class Histone_Markers {
	use APPLES_Datatypes qw(HistoneModification);
	use Parameters;
	
	has 'tissue'		=> 	(is => 'rw', isa => 'Str');
	has 'chip'			=> (is => 'rw', isa => 'Obj');
	has 'aa_of_histone_modification'	=>	(is => 'rw', isa => 'Str');
	has 'histone_modification_type'	=>	(is => 'rw', isa => HistoneModification);

} # Histone_Markers #
