### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Evidence Class ###
# Evidence may be experimental or computational, and should have an associated score

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class Evidence {
	use Parameters;
	
	has 'evidence_type' => (is => 'rw', isa => 'Str');
	has 'evidence_description' => (is => 'rw', isa => 'Str');
	has 'evidence_score' => (is => 'rw', isa => 'Int');

} # Evidence #
