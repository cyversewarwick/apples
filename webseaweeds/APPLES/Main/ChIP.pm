### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ChIP Class: parent class for ChIP (Chromatin Immunoprecipitation) data ###

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class ChIP {
	use Parameters;
	
	has 'protein_id'	=> (is => 'rw', isa => 'Str');
	has 'antibody'		=> (is => 'rw', isa => 'Str'); 
	has 'bisi_list'		=> (is => 'rw', isa => 'HashRef[BiSi]'); 

} # ChIP #
