### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReMo_Data class ###
# technical class only used in Bundling

use MooseX::Declare;

class ReMo_Data_Object {
	use Cons_Elem_Raw_Data;

	has 'firstsequence' => (is => 'rw', isa => 'Sequence');
	has 'secondsequence' => (is => 'rw', isa => 'Sequence');
	has 'windowlength' => (is => 'rw', isa => 'Int');
	has 'stepwidth1' => (is => 'rw', isa => 'Int');
	has 'stepwidth2' => (is => 'rw', isa => 'Int');
	has 'windowpairs' => (is => 'rw', isa => 'ArrayRef[WindowPair]'); 
  # --HOLDS A RECORD

} # ReMo_Data_Object #

class Sequence {
	use APPLES_Datatypes qw(APPLESSpeciesName);

	has 'maskedsequence' => (is => 'rw', isa => 'Str');
	has 'unmaskedsequence' => (is => 'rw', isa => 'Str');
	has 'species' => (is => 'rw', isa => APPLESSpeciesName);
	has 'id' => (is => 'rw', isa => 'Str'); # REDUNDANT - only for print_bundles when debugging
  # --HOLDS SEQUENCE INFORMATION

} # Sequence #

class WindowPair {
	has 'offset1' => (is => 'rw', isa => 'Int');
	has 'offset2' => (is => 'rw', isa => 'Int');
	has 'score' => (is => 'rw', isa => 'Num');
  # --HOLDS WINDOW LOCATIONS AND SCORE OF WINDOW PAIR

} # WindowPair #

