### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genomic_Data_Miner Class ###
use MooseX::Declare;

class Genomic_Data_Miner {
	use Parameters;
	use APPLES_Datatypes;
	use Statistics;
	use Data::Dumper;
	use Generic_Pattern_Matching_Model;	
	use BiSi;
	use General_Utilities;
	use WM_Scores;

	my $GU = General_Utilities->new();

	### all methods that were here before have found better homes now
	### keep this class anyway as a place for higher level scripts around the
	### topic of "genomic data mining"

}; # Genomic_Data_Miner #

	
