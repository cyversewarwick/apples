### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Version_Number Class ###
use MooseX::Declare;

class Version_Number {

	has 'major' => (is => 'rw', isa => 'Int');
	has 'minor' => (is => 'rw', isa => 'Int');
	has 'patch' => (is => 'rw', isa => 'Int');
	#has 'svn' => (is => 'rw', isa => 'Int');
	
	method compare_version_number(Version_Number $version_to_compare) {
		
		die 'method not implemented yet!';
		# (method may not be required)
		# would take another Version_Number as input, compare to self, and work out which is greater
	} # compare_version_number #

} # Version_Number #
