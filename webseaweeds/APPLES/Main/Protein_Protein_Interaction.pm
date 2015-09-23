### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Protein_Protein_Interaction Class: describes interactions between 2 or more proteins ###

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class Protein_Protein_Interaction {
	use Parameters;
	
	has 'PPI_members' => (is => 'rw', isa => 'ArrayRef[Protein_ID]'); # list proteins in the interaction

} # Protein_Protein_Interaction #

# Does this need to include ordering of binding (e.g. A+B, then C)?
#							sites of binding?
#							evidence of binding?
