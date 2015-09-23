## (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

=head1 Class XML_Reader

Class for reading serializable objects from XML output.

=cut

use MooseX::Declare;

class XML_Reader {
	use DataInterface::XML::Serializable;
	use XML::Simple qw (XMLin);
	
	method read_serializable_objects(Str $filename) {
		
	}
	
}
