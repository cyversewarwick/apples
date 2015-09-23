### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Reg_Loc class ###
# Describes a regulated locus, for example: the site (may be exactly defined TSS, or a more general region)
# where polymerase begins transcription
# A regulated locus could be an Ensembl gene, an unannotated gene, miRNA etc.

use lib "/home/grannysmith/webseaweeds/APPLES/Main/";
use MooseX::Declare;
use Serialization::Serializable;

class Reg_Loc extends Serialization::Serializable {
	use Parameters;
	use Genome_DB_Utilities;
	use constant {FALSE => 0,
		      TRUE => 1};	

        use Moose;
	use APPLES_Datatypes qw(StrandType RegLocLocusType);
	has 'genome_db'	=>	(is => 'ro', isa => 'Genome_Sequence_Database_Parameters', required => 1);
	has 'coord_sys_name' => (is => 'ro', isa => 'Str');
	has 'region'	=>	(is => 'ro', isa => 'Str', required => 1); # must have a region identifier
	has 'position'	=>	(is => 'ro', isa => 'Int', required => 1); # must have a position
	has 'strand'	=>  (is => 'rw', isa => StrandType, required => 1);
	has 'locus_type'			=> (is => 'rw', isa => RegLocLocusType);
	has 'gene_ID'               => (is => 'rw', isa => 'Str');
	has 'transcript_ID'			=> (is => 'rw', isa => 'Str');

	method get_length_of_regulation_target() {
	    # regulation target is transcript in case of genes, miRNA in case of miRNAs, may
	    # not be sensible notion for all types of Reg_Locs
	    
	    my $result;
	    if (!defined $self->locus_type) {
		die 'locus type must be defined for this operation to work.';
	    } elsif ($self->locus_type eq 'gene') {		
		if (defined $self->gene_ID) {
		    my $gdbu = Genome_DB_Utilities->new();
		    $result = $gdbu->get_gene_length($self->genome_db,$self->gene_ID);
		} else {
		    die 'gene-ID must be defined for this operation to work.';
		}
	    } else {
		die 'method not implemented for this locus type yet: '.$self->locus_type;
	    }	 
	    return $result;
	} # get_length_of_regulation_target #

} # Reg_Loc #
