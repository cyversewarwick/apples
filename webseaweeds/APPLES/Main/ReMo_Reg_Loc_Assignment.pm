### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReMo_Reg_Loc_Assignment: a pair of a ReMo and a Reg_Loc ###

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class ReMo_Reg_Loc_Assignment {
	use Parameters;
	use Reg_Loc;
	use ReMo;
	
	use Moose;
	has 'remo' => (is => 'rw', isa => 'ReMo', required => 1);
	has 'reg_loc' => (is => 'rw', isa => 'Reg_Loc', required => 1);
	has 'assignment_score' => (is => 'rw', isa => 'Int', required => 1);

} # ReMo_Reg_Loc_Assignment #
