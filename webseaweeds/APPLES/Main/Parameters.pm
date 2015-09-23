### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Parameter classes: Defines parameter structures and default parameters for APPLES ###

use MooseX::Declare;

class Parameters {
	
} # Parameters #

class Genome_Sequence_Database_Parameters {
	use APPLES_Datatypes qw (APPLESSpeciesName);

	has 'orthology_methods' => (is => 'rw', isa => 'ArrayRef[Orthology_Method_And_Species_Restriction]'); # that can be used when this database is the TARGET
	has 'dbname' => (is => 'ro', isa => 'Str', required => 1); # essential for ensembl for connecting AND proper caching in all cases. If FASTA, needs to be a unique filename

} # Genome_Sequence_Database_Parameters #

class Genome_CDS_GOterm_Data {
	has 'go_process' => (is => 'rw', isa => 'ArrayRef[Str]');
	has 'go_component' => (is => 'rw', isa => 'ArrayRef[Str]');
	has 'go_function' => (is => 'rw', isa => 'ArrayRef[Str]');
} # Genome_CDS_GOterm_Data #


class FASTA_Sequence_Database_Parameters extends Genome_Sequence_Database_Parameters{
	use APPLES_Datatypes qw (APPLESSpeciesName);
	use General_Utilities;

	has 'filename' => (is => 'ro', required => 1, isa => 'Str'); # must be the full pathname
	                                                             # Sascha: had problems with Bioperl when filename ended in '.txt' instead of '.fasta', so better make it '.fasta'
	has 'natural_species_name' => (is => 'ro', isa => APPLESSpeciesName); # optional, but will be necessary if dealing with phylogenetic trees
	
	method get_md5sum() {
	  my $GU = General_Utilities->new();
	  my $filepath = $self->filename;
	  my $result = $GU->md5sum($filepath);
	  return $result;
	} # get_md5sum #
 
} # FASTA_Sequence_Database_Parameters extends Genome_Sequence_Database_Parameters #

class Ensembl_Database_Parameters extends Genome_Sequence_Database_Parameters {
	use APPLES_Datatypes qw (APPLESSpeciesName Boolean EnsemblLocation);
	use constant {FALSE => 0,
		      TRUE  => 1};	

	has 'alias' => (is => 'ro', required => 1, isa => APPLESSpeciesName);
	has 'location' => (is => 'ro', isa => EnsemblLocation, required => 1);

} # Ensembl_Database_Parameters extends Genome_Sequence_Database_Parameters #

class Genbank_Sequence_Database_Parameters extends Genome_Sequence_Database_Parameters {
    use APPLES_Datatypes qw (Boolean APPLESSpeciesName);
    use constant {FALSE => 0,
		  TRUE	=> 1};	

    has 'filenames' => (is => 'rw', required => 1, isa => 'ArrayRef', trigger => \&sort_filenames); # in addition to, but can be same as, 'dbname' attribute. 
    has 'block_sort_filename_trigger' => (is => 'rw', isa => Boolean, required => 0, default => FALSE);
    has 'natural_species_name' => (is => 'ro', isa => APPLESSpeciesName); # optional, but will be necessary if dealing with phylogenetic trees
    
    # By a wonderful peice of Perlian logic, a Moose trigger has to refer to a sub, not a method, despite the fact that
    # Moose has evolved to use methods rather than subs within classes. 
    # (Yes Nigel, this is because 'trigger' is old-skool Moose that only used subs, but we generally use an extension, MooseX::Declare, which is what enables us to use 'method'. 
    # It makes perfect sense!)

    sub sort_filenames() {
	my ($self, $arg) = @_; 
	
	if (!$self -> block_sort_filename_trigger()) {
	     $self -> block_sort_filename_trigger(TRUE);
	     my $unsorted_filenames = $self->filenames();
	     my @sorted_filenames = sort @{$unsorted_filenames};
	     $self->filenames(\@sorted_filenames);
	     $self -> block_sort_filename_trigger(FALSE);
	}
    } # sort_filenames #
    
} # Genbank_Sequence_Database_Parameters # 

class Compara_Database_Parameters {
  has 'host'    => (is => 'ro', isa => 'Str', required => 1);
  has 'user'    => (is => 'ro', isa => 'Str', required => 1);
  has 'port'    => (is => 'ro', isa => 'Int', required => 1);
  has 'species' => (is => 'ro', isa => 'Str', required => 1);
  has 'dbname'  => (is => 'ro', isa => 'Str', required => 1);

} # Compara_Database_Parameters #

class Genomic_Interval_Parameters {
	use APPLES_Datatypes qw (GenomicRegionType StrandType);

	has 'geneid' 		  		=> (is => 'rw', required => 0, isa => 'Str');
	has 'genomic_region_type' 	=> (is => 'ro', required => 1, isa => GenomicRegionType);
	has 'database'				=> (is => 'ro', required => 1, isa => 'Genome_Sequence_Database_Parameters');
	has 'region'                => (is => 'ro', required => 1, isa => 'Str'); 
	has 'five_prime_pos'		=> (is => 'ro', required => 1, isa => 'Int');
	has 'three_prime_pos'		=> (is => 'ro', required => 1, isa => 'Int');
	has 'strand'			=> (is => 'ro', required => 1, isa => StrandType);

} # Genomic_Interval_Parameters #

class Window_Pair_Algorithm_Parameters {
	use APPLES_Datatypes qw (AlignmentAlgorithm);

} # Window_Pair_Algorithm_Parameters #

# not sure about these 2 classes
class Ott_Algorithm_Parameters extends Window_Pair_Algorithm_Parameters {
	has 'stepwidth1' => (is => 'ro', isa => 'Int', required => 1, default => 1);
	has 'stepwidth2' => (is => 'ro', isa => 'Int', required => 1, default => 1);
	has 'windowlength' => (is => 'ro', isa => 'Int', required => 1, default => 100); 
	has 'cutoff_for_uninteresting_alignments' => (is => 'ro', isa => 'Int', default => 55);
	has 'nw_matchscore' => (is => 'ro', isa => 'Num', default => '1.0');
	has 'nw_mismatchscore' => (is => 'ro', isa => 'Num', default => '0.0');
	# the next three values must fulfill: (effectivegapscore/effectivematchscore) = gapscore
	has 'gapscore' => (is => 'ro', isa => 'Num', default => '-0.5'); # must be negative
	has 'effectivematchscore' => (is => 'ro', isa => 'Num', default => '2');
	has 'effectivegapscore' => (is => 'ro', isa => 'Num', default => '-1');
	has 'maxconsregionsoffset' => (is => 'ro', isa => 'Int', default => 3000);
	has 'maxconsregionsperbase' => (is => 'ro', isa => 'Num', default => '0.4');

} # Ott_Algorithm_Parameters #

class Seaweed_Algorithm_Parameters extends Window_Pair_Algorithm_Parameters {
	has 'stepwidth' => (is => 'rw', isa => 'Int', required => 1, default => 1);
	has 'windowlength' => (is => 'ro', isa => 'Int', required => 1, default => 100); 
	has 'cutoff_for_uninteresting_alignments' => (is => 'ro', isa => 'Int', default => 57);
} # Seaweed_Algorithm_Parameters #

class Sequence_Parameters {
	use APPLES_Datatypes qw (Boolean SequenceRegionType PositiveInt);
	use constant {FALSE => 0,
		      TRUE  => 1};	

	has 'region' => (is => 'ro', isa => SequenceRegionType, required => 1);
	has 'min_length_to_return' => (is => 'ro', isa => PositiveInt, required => 1, default => 100);
	has 'max_length_to_search' => (is => 'ro', isa => PositiveInt, required => 1, default => 10000);
	has 'stop_at_neighbouring_gene' => (is => 'ro', isa => Boolean, required => 1, default => TRUE);
	has 'ignore_pseudogenes' => (is => 'ro', isa => Boolean, required => 1, default => FALSE);
	has 'ignore_RNAgenes' => (is => 'ro', isa => Boolean, required => 1, default => FALSE);
} # Sequence_Parameters #

class ReMo_Constructor_Parameters {
	use APPLES_Datatypes qw (Boolean);
	use constant {FALSE => 0, TRUE => 1};
    has 'ignore_pseudogenes'    => (is => 'ro', required => 1, isa => Boolean, default => TRUE);
    has 'ignore_RNAgenes'           => (is => 'ro', required => 1, isa => Boolean, default => FALSE); # [HM] introduced in May 2010
} # ReMo_Constructor_Parameters #

class ReMo_Core_Promoter_Constructor_Parameters	extends ReMo_Constructor_Parameters {
	use APPLES_Datatypes qw (Boolean PositiveInt);
    use constant {FALSE => 0, TRUE => 1};
	has 'length'			=> (is => 'ro', required => 1, isa => PositiveInt, default => 500);
	has 'stop_at_neighbouring_gene' => (is => 'ro', required => 1, isa => Boolean, default => TRUE);
	has 'minimum_length'            => (is => 'ro', required => 1, isa => PositiveInt, default => 100);
} # ReMo_Core_Promoter_Constructor_Parameters #

class ReMo_Set_Constructor_Parameters {
	has 'sequence_databases_to_use_for_homologs' => (is => 'ro', required => 1, isa => 'ArrayRef[Genome_Sequence_Database_Parameters]');
	#has 'orthology_finding_parameters' => (is => 'ro', required => 1, isa => 'Generic_Orthology_Finding_Parameters'); # redundant

} # ReMo_Set_Constructor_Parameters #

class ReMo_Set_Core_Promoter_Constructor_Parameters extends ReMo_Set_Constructor_Parameters {
	has 'remo_core_promoter_constructor_parameters' => (is => 'ro', required => 1, isa => 'ReMo_Core_Promoter_Constructor_Parameters');

} # ReMo_Set_Core_Promoter_Constructor_Parameters #

class ReMo_Set_Phylogenetic_Constructor_Parameters extends ReMo_Set_Constructor_Parameters {
	use APPLES_Datatypes qw (Boolean);
	use constant {FALSE => 0,
		      TRUE	=> 1};	

	has 'window_pair_algorithm_parameters' => (is => 'ro', isa => 'Window_Pair_Algorithm_Parameters', required => 1);
	has 'star_bundler_parameters' => (is => 'ro', isa => 'Star_Bundler_Parameters', required => 1);
	has 'sequence_parameters' => (is => 'ro', isa => 'Sequence_Parameters', required => 1);

} # ReMo_Set_Phylogenetic_Constructor_Parameters #

class ReRe_Constructor_Parameters {
	has 'remo_set_constructor_parameters'	=> (is => 'ro', required => 1, isa => 'ArrayRef[ReMo_Set_Constructor_Parameters]');

} # ReRe_Constructor_Parameters #

class ReRe_Set_Constructor_Parameters {
	has 'rere_constructor_parameters' => (is => 'ro', required => 1, isa => 'ReRe_Constructor_Parameters');

} # ReRe_Set_Constructor_Parameters #

class Orthologous_Genomic_Interval_Group_Region_to_Take {
	use APPLES_Datatypes qw (OrthologousGenomicIntervalGroupRegionToTake);

	has 'region' => (is => 'ro', isa => OrthologousGenomicIntervalGroupRegionToTake);

} # Orthologous_Genomic_Interval_Group_Region_to_Take #

#class BiFa_Parameters_Core {
#	has 'threshold' => (is => 'rw', isa => 'Int', required => 1);
#	has 'wm_sets' => (is => 'rw', isa => 'ArrayRef[Str]', required => 1);
#	has 'scoring_method' => (is => 'rw', isa => 'Str', required => 1); 
#	has 'wm_query' => (is => 'rw', isa => 'Str', required => 1); # default is all, but user can specify subset
#}

class BiFa_Parameters_Core {
	use APPLES_Datatypes qw (BiFaAlgorithm);

	#has 'algorithm' => (is => 'ro', isa => BiFaAlgorithm, required => 1);
	#has 'significance_threshold' => (is => 'ro', isa => 'Num', required => 1, default => 0.01);

} # BiFa_Parameters_Core #

class BiFa_Parameters_Multiple {
	has 'bifa_parameters_core' => (is => 'ro', isa => 'BiFa_Parameters_Core', required => 1);
	has 'phylogenetic_threshold' => (is => 'rw', isa => 'Int', required => 1);

} # BiFa_Parameters_Multiple #

class MEME_Parameters {
    use APPLES_Datatypes qw (Boolean);
    use constant {FALSE => 0,
		  TRUE	=> 1};	

    has '-oc' => (is => 'rw', isa => 'Str'); # name of directory for output files will replace existing directory
    has '-dna' => (is => 'rw', isa => Boolean, default => TRUE); #sequences use DNA alphabet
    has '-mod' => (is => 'rw', isa => 'Str', default => 'zoops');		# [-mod oops|zoops|anr]		distribution of motifs
    has '-nmot ifs' => (is => 'rw', isa => 'Int', default => '3');	# [-nmotifs <nmotifs>]		maximum number of motifs to find
    has '-evt' => (is => 'rw', isa => 'Int');		# [-evt <ev>]				stop if motif e-value greater than <evt>
    has '-nsites' => (is => 'rw', isa => 'Int');	# [-nsites <sites>]			number of sites for each motif
    has '-minsites' => (is => 'rw', isa => 'Int');	# [-minsites <minsites>]	minimum number of sites for each motif
    has '-maxsites' => (is => 'rw', isa => 'Int');	# [-maxsites <maxsites>]	maximum number of sites for each motif
    has '-w' => (is => 'rw', isa => 'Int');			# [-w <w>]					motif width
    has '-minw' => (is => 'rw', isa => 'Int');		# [-minw <minw>]			minumum motif width
    has '-maxw' => (is => 'rw', isa => 'Int');		# [-maxw <maxw>]			maximum motif width
    has '-revcomp' => (is => 'rw', isa => Boolean, default => TRUE); 	# [-revcomp]				allow sites on + or - dna strands
    has '-pal' => (is => 'rw', isa => Boolean, default => FALSE);		# [-pal]					force palindromes (requires -dna)
    has '-maxiter' => (is => 'rw', isa => 'Int', default => 5000);	# [-maxiter <maxiter>]		maximum em iterations to run
    has '-distance' => (is => 'rw', isa => 'Int');	# [-distance <distance>]	em convergence criterion
    has '-cons' => (is => 'rw', isa => 'Str'); 		# [-cons <cons>]			consensus sequence to start EM from
    has '-maxsize' => (is => 'rw', isa => 'Int', default => 500000); # must be above size of sequence file for MEME to accept input
	has '-bfile' => (is => 'rw', isa => 'Str');		# [-bfile <file> ] bckground model file
} # MEME_Parameters #

class Empirical_Scoring_Method {
	
}

class Biobase_Additive extends Empirical_Scoring_Method {
	has 'entropy_param' => (is =>'ro', isa => 'Num', required => 1, default => 4);
}

class Multiplicative extends Empirical_Scoring_Method {
	
}

class Generic_WM_Scoring_Parameters {

} # Generic_WM_Scoring_Parameters #

class BiFa_WM_Scoring_Parameters extends Generic_WM_Scoring_Parameters {
	use APPLES_Datatypes qw (BiFaAlgorithm);
	
	has 'algorithm' => (is => 'ro', isa => BiFaAlgorithm, default => 'ALG_OTT');
		
} # BiFa_WM_Scoring_Parameters #

class Empirical_WM_Scoring_Parameters extends Generic_WM_Scoring_Parameters {
	
	has 'scoring_method' => (is =>'ro', isa => 'Empirical_Scoring_Method', required => 1);
	
} # Empirical_WM_Scoring_Parameters #

class Genome_Based_Empirical_WMS_Parameters extends Empirical_WM_Scoring_Parameters{
	has 'genome_db' => (is => 'ro', isa => 'Genome_Sequence_Database_Parameters', required => 1);

} # Genome_Based_Empirical_WMS_Parameters #

class Model_Based_Empirical_WMS_Parameters extends Empirical_WM_Scoring_Parameters{	
	has 'sequence_model_source' => (is => 'ro', isa => 'Sequence_Model_Source', required => 1);
	has 'training_sequence_length' => (is => 'ro', isa => 'Int', required => 1, default => '100000000');# Default seq len is 100 million
	has 'burn_in' => (is => 'ro', isa => 'Int', default => 100);
	has 'seed' => (is => 'ro', isa => 'Int', required => 1, default => 13);

} # Model_Based_Empirical_WMS_Parameters #

class Sequence_Model_Source {
	
} # Sequence_Model_Source #

class Partial_Sequence_Model_Learning_Params extends Sequence_Model_Source {
	use APPLES_Datatypes qw(GenomicSequenceClass);
	
	has 'model' => (is => 'ro', isa => 'Generic_Statistical_Sequence_Model_Type', required => 1);
	has 'learning_restrictions' => (is => 'rw', isa => GenomicSequenceClass, required => 1);

} # Partial_Sequence_Model_Learning_Params #

class Full_Sequence_Model_Learning_Params extends Partial_Sequence_Model_Learning_Params {	
	has 'genome_db' => (is => 'ro', isa => 'Genome_Sequence_Database_Parameters', required => 1);

} # Full_Sequence_Model_Learning_Params #

class Actual_Sequence_Model extends Sequence_Model_Source{	
	has 'sequence_model' => (is => 'ro', isa => 'Sequence_Model', required => 1);

} # Actual_Sequence_Model #

class Generic_Statistical_Sequence_Model_Type {
		
} # Generic_Statistical_Sequence_Model_Type #

class Markov_Model extends Generic_Statistical_Sequence_Model_Type {	
	has 'order' => (is => 'ro', isa => 'Int', default => 3);

} # Markov_Model #

class WM_Pseudocount_Parameters{
	use APPLES_Datatypes qw(PseudocountType NonNegativeNum);
	
	has 'pseudocount_type' => (is => 'ro', isa => PseudocountType, required => 1);
	has 'p' => (is => 'ro', isa => 'NonNegativeNum', required => 1);
} # WM_Pseudocount_Parameters #

class Job_Mode_Parameters {
	use APPLES_Datatypes qw (RunningMode);

	has 'mode' => (is => 'rw', isa => RunningMode, required => 1);
} # Job_Mode_Parameters #

class Star_Bundler_Parameters {
    use APPLES_Datatypes qw (Probability);

    has 'overlap_tolerance' => (is => 'ro', isa => 'Int', required => 1, default => 20);
    has 'belief_value' => (is => 'ro', isa => Probability, required => 1, default => 0.1);
    has 'partial_threshold_matrix' => (is => 'ro', isa => 'Partial_Threshold_Matrix', required => 1);
} # Star_Bundler_Parameters #

class Generic_Orthology_Finding_Parameters {

} # Generic_Orthology_Finding_Parameters #

class Compara_Orthology_Finding_Parameters extends Generic_Orthology_Finding_Parameters {
    has 'compara_db' => (is => 'ro', isa => 'Str', required => 1); # name of compara database

} # Compara_Orthology_Finding_Parameters #

class Synteny_Orthology_Finding_Parameters extends Generic_Orthology_Finding_Parameters {

} # Synteny_Orthology_Finding_Parameters #

class RBH_Orthology_Finding_Parameters extends Generic_Orthology_Finding_Parameters {

} # RBH_Orthology_Finding_Parameters #

class Random_Assignment_Orthology_Finding_Parameters extends Generic_Orthology_Finding_Parameters {
  has 'seed' => (is => 'rw', isa => 'Int', required => 1, default => 13); # random seed for srand
} # Random_Assignment_Orthology_Finding_Parameters #

class Orthology_Method_And_Species_Restriction {
  has 'generic_orthology_finding_parameters' => (is => 'ro', isa => 'Generic_Orthology_Finding_Parameters', required => 1);
  has 'allowed_source_dbnames' => (is => 'ro', isa => 'ArrayRef', required => 1); # filter for acceptable species. Was species, is now dbnames, to be generic
} # Orthology_Method_And_Species_Restriction #

class Alignment_Parameters {
} # Alignment_Parameters #

class NCBI_BLAST_Parameters extends Alignment_Parameters {
    use APPLES_Datatypes qw (Boolean StrandChoice);
    use constant {FALSE => 0,
		  TRUE	=> 1};
    use Bio::Tools::Run::StandAloneBlast;

    has 'word_length'    => (is => 'rw', isa => 'Int', required => 1);
    has 'alignments_to_do' => (is => 'rw', isa => 'Int', required => 1, default => 200);
    has 'alignments_to_show' => (is => 'rw', isa => 'Int', required => 1, default => 200);
    has 'filter' => (is => 'rw', isa => Boolean, required => 1, default => FALSE);
    has 'program' => (is => 'rw', isa => 'Str', required => 1);
    has 'strands' => (is => 'rw', isa => StrandChoice, required => 1);

    method create_blast_factory(Int $mtype, Str $outputfilename, Str $BLASTdatabase) {
	my $GU = General_Utilities->new();
	$GU->user_info(3,"creating blast factory\n");
	my @factoryparameters = (
	    program     => $self->program,
	    _READMETHOD => "Blast",
	    m => $mtype,
	    o => $outputfilename,
	    d => $BLASTdatabase,
	    W => $self->word_length,
	    v => $self->alignments_to_do,
	    b => $self->alignments_to_show,
	);
	my $factory = Bio::Tools::Run::StandAloneBlast->new(@factoryparameters);
	if ($self->filter) {
	    $factory->F('T');
	}
	else {
	    $factory->F('F');
	}
	if ($self->strands eq 'both') {
	    $factory->S(3);
	}
	elsif ($self->strands eq 'positive') {
	    $factory->S(1);
	}
	elsif ($self->strands eq 'negative') {
	    $factory->S(2);
	}
	else {
	    die 'behaviour for this strand choice is not defined.';
	}
	return $factory;
    } # create_blast_factory #

} # NCBI_BLAST_Parameters #

class Job_Parameters {
    use APPLES_Datatypes qw (Boolean);
    use constant {FALSE => 0,TRUE => 1};

    has 'memory_requirement_high' => (is => 'ro', isa => Boolean, required => 1); 
    has 'wall_time_estimate' => (is => 'ro', isa => 'Num', required => 1); # seconds
    has 'recache' => (is => 'ro', isa => Boolean, default => FALSE);
    has 'cache_in_memory' => (is => 'ro', isa => Boolean, default => FALSE); # keep result in memory when retrieving from cache database (to avoid repetitive retrieval)

} # Job_Parameters #

class P_Value_Sampling_Likelihood_Conversion_Parameters {
    
    # to be developed

} # P_Value_Sampling_Likelihood_Conversion_Parameters #

class ReCo_Parameters {
    use APPLES_Datatypes qw (ProteinProteinNetworks Probability PositiveInt);

    has 'p_p_network' => (is => 'ro', isa => ProteinProteinNetworks, required => 1);
    has 'minimum_node_probability' => (is => 'ro', isa => Probability, required => 1, default => 0.01);
    has 'number_of_sampling_iterations' => (is => 'rw', isa => PositiveInt, required =>1, default => 10000);
    has 'max_number_of_sites' => (is => 'ro', isa => PositiveInt, required => 1);  # for binomial overrepresentation test
    has 'conversion_parameters' => (is => 'ro', isa => 'P_Value_Sampling_Likelihood_Conversion_Parameters', required => 1);
    has 'scale_number_of_experiments' => (is => 'ro', isa => 'Num', required => 1); # for logistic function (edge weights)
    has 'scale_number_of_experimental_methods' => (is => 'ro', isa => 'Num', required => 1); # for logistic function (edge weights)
    has 'weight_for_two_logistic_functions' => (is => 'ro', isa => Probability, required => 1); # weighting when integrating two logistic functions (edge weights)

} # ReCo_Parameters #
