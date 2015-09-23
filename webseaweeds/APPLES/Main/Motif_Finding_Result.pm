### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Motif_Finding_Result Class ###

use MooseX::Declare;

class Motif_Finding_Result {
	use Parameters;
	use Generic_Sequence_Pattern;
} # Motif_Finding_Result #

class MEME_Motif_Finding_Result extends Motif_Finding_Result {
    use APPLES_Datatypes qw (Probability);

    has 'num_sites' => (is => 'rw', isa => 'Int', required => 1);
    has 'motif_width' => (is => 'rw', isa => 'Int', required => 1);
    has 'occurrence_ratio' => (is => 'rw', isa => 'Value', required => 1);
    has 'e_value' => (is => 'rw', isa => 'Str', required => 1);
    has 'positional_bias_pvalue' => (is => 'rw', isa => Probability, required => 0);
    has 'strand_bias_pvalue' => (is => 'rw', isa => Probability, required => 0);
    has 'WM' => (is => 'rw', isa => 'Generic_Weight_Matrix', required => 1);
	has 'which_genes' => (is => 'rw', isa => 'ArrayRef', required => 0, default => sub {[]});
    #has 'closest_matching_WM' => (is => 'rw', isa => 'HashRef[]', required => 0); Unsure as to whether this attribute is needed
		
} # MEME_Motif_Finding_Result extends Motif_Finding_Result #
