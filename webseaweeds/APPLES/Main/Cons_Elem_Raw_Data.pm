### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Cons_Elem_Raw_Data class ###
# technical class playing a role in Bundling only
use MooseX::Declare;

class Cons_Elem_Raw_Data {
	use APPLES_Datatypes qw (Boolean);
	use constant {FALSE => 0,
		      TRUE	=> 1};
	has 'conservation' => (is => 'rw', isa => 'Num');
	has 'repeatratio' => (is => 'rw', isa => 'Num');
	has 'startbase' => (is => 'rw', isa => 'Int');
	has 'endbase' => (is => 'rw', isa => 'Int');
	has 'targetstartbase' => (is => 'rw', isa => 'Int', required => 0);# --NOT DEFINED FOR REFERENCE SEQUENCE
	has 'targetendbase' => (is => 'rw', isa => 'Int', required => 0); # --NOT DEFINED FOR REFERENCE SEQUENCE
	has 'beliefscoreisdefined' => (is => 'rw', isa => Boolean);
	has 'beliefscore' => (is => 'rw', isa => 'Num');
	has 'integrate' => (is => 'rw', isa => 'Num');
} # Cons_Elem_Raw_Data #

