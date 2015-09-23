### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ChIP_Chip Class: describes ChIP-chip data

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class ChIP_Chip {
	use Parameters;
	
	has 'protein_id'  	=> (is => 'rw', isa => 'Str');
	has 'bisi_list'		=> (is => 'rw', isa => 'HashRef[BiSi]');			
	has 'evidence' 		=> (is => 'rw', isa => 'Str');		 # e.g. what microarray experiment

} # ChIP_Chip #

