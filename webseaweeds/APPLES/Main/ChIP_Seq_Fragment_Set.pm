### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ChIP_Seq_Fragment_Set Class ###

# this class is not implemented yet, just illustrating an idea

use MooseX::Declare;

class ChIP_Seq_Fragment_Set {
	use Parameters;
	
	has 'genomic_interval'		=>	(is => 'rw', 
						 isa => 'Genomic_Interval',
						 required => 1); 				# must have a Genomic_Interval
	has 'frequency_plot'		=>	(is => 'rw', isa => 'HashRef');
	has 'protein_id'			=>	(is => 'rw', isa => 'Str');

} # ChIP_Seq_Fragment_Set #

