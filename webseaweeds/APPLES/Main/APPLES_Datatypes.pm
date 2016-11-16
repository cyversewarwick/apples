### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### APPLES_Datatypes Class ###
package APPLES_Datatypes; 
use MooseX::Types
	-declare => [qw(
		Boolean 
		GenomicRegionType
		OrthologousGenomicIntervalGroupRegionToTake
		AminoAcid
		HistoneModification
		LocusType
		StrandType
		ScientificSpeciesName
		APPLESSpeciesName
		RegLocLocusType
		HomologyType
		GenomicDirection
		BiFaAlgorithm   
		SequenceChemicalType
		RunningMode
		AlignmentAlgorithm
		JobHandlerIdentifier
		Status
		WebService
		EvolutionaryTreeKingdom
		SequenceRegionType
		PrintLevel
		WeightMatrixSource
		ScoringModel
		JobMode
		ResultType
                ListenerType
		MethodType
		CommandLine
		EnsemblLocation
		StrandChoice
                SequenceFunctionalType
		TranscriptChoice
		GenomicSequenceClass
		PerseveranceCategory
		PositiveNum
		NonNegativeNum
		PseudocountType
		PositiveInt
		NonNegativeInt
                Probability
                ProteinProteinNetworks
                WorkingSequenceType
                IntPercentage        
	)]; # list all APPLES data types;
use MooseX::Types::Moose 'Int';
use MooseX::Types::Moose 'Str';
use MooseX::Types::Moose 'Num';

use constant {FALSE => 0,
	      TRUE	=> 1};

subtype Boolean,
	as Int,
	where { ( $_ == 0 ) || ( $_ == 1 ) },
	message { "$_ invalid boolean value" };

subtype GenomicRegionType,
	as Str,
	where { ( $_ eq 'upstream_region' ) || ( $_ eq 'downstream_region' ) || ( $_ eq 'mRNA' ) 
	|| ( $_ eq 'gene_space' ) || ( $_ eq '5primeUTR' ) || ( $_ eq '3primeUTR' ) },
	message { "$_ invalid GenomicRegionType" };
	
subtype OrthologousGenomicIntervalGroupRegionToTake,
	as Str,
	where { ( $_ eq 'wide' ) || ( $_ eq 'narrow' ) }, # rule to deal with case where syntenic boundary is undefined
	message { "$_ invalid OrthologousGenomicIntervalGroupRegionToTake" };
	
subtype AminoAcid,
	as Str,
	where { ( $_ eq 'alanine' ) || ( $_ eq 'arginine' ) || ( $_ eq 'asparagine' ) || ( $_ eq 'aspartic acid' ) || ( $_ eq 'cysteine' ) 
	|| ( $_ eq 'glutamic_acid' ) || ( $_ eq 'glutamine' ) || ( $_ eq 'glycine' ) || ( $_ eq 'histidine' ) 
	|| ( $_ eq 'isoleucine' ) || ( $_ eq 'leucine' ) || ( $_ eq 'lysine' ) || ( $_ eq 'methionine' ) || ( $_ eq 'phenylalanine' ) 
	|| ( $_ eq 'proline' ) || ( $_ eq 'serine' ) || ( $_ eq 'threonine' ) || ( $_ eq 'tryptophan' ) || ( $_ eq 'tyrosine' ) || ( $_ eq 'valine' ) },
	message { "$_ invalid AminoAcid" };
	
subtype HistoneModification,
	as Str,
	where { ( $_ eq 'monomethylation' ) || ( $_ eq 'dimethylation' ) || ( $_ eq 'trimethylation' ) || ( $_ eq 'methylation' ) || ( $_ eq 'acetylaction' ) 
	|| ( $_ eq 'citrullination' ) || ( $_ eq 'phosphorylation' ) || ( $_ eq 'Sumoylation' ) 
	|| ( $_ eq 'ubiquitination' ) || ( $_ eq 'ADP-ribosylation' ) },
	message { "$_ invalid HistoneModification" };

subtype LocusType,
	as Str,
	where { ( $_ eq 'gene' ) ||  ( $_ eq 'microRNA' ) || ( $_ eq 'unannotated_gene' ) 
	|| ( $_ eq 'alternatively_spliced_exon' ) },
	message { "$_ invalid LocusType" };

subtype StrandType,
	as Str,
	where { ( $_ eq 'positive' ) || ( $_ eq 'negative' ) },
	message { "$_ invalid StrandType" };

my @Aensembl_scientific_names = ('Homo sapiens', 'Pan troglodytes', 'Gorilla gorilla', 'Pongo pygmaeus', 'Macaca mulatta', 'Papio hamadryas', 'Colobus guereza', 'Tarsius syrichta', 'Aotus trivirgatus', 'Callithrix jacchus', 'Callicebus moloch', 'Microcebus murinus', 'Otolemur garnettii', 'Tupaia belangeri', 'Oryctolagus cuniculus', 'Ochotona princeps', 'Spermophilus tridecemlineatus', 'Mus musculus', 'Rattus norvegicus', 'Dipodomys ordii', 'Cavia porcellus', 'Bos taurus', 'Tursiops truncatus', 'Vicugna pacos', 'Sus scrofa', 'Canis familiaris', 'Felis catus', 'Equus caballus', 'Myotis lucifugus', 'Pteropus vampyrus', 'Erinaceus europaeus', 'Sorex araneus', 'Procavia capensis', 'Loxodonta africana', 'Echinops telfairi', 'Choloepus hoffmanni', 'Dasypus novemcinctus', 'Macropus eugenii', 'Monodelphis domestica', 'Ornithorhynchus anatinus', 'Gallus gallus', 'Taeniopygia guttata', 'Anolis carolinensis', 'Xenopus tropicalis', 'Tetraodon nigroviridis', 'Takifugu rubripes', 'Oryzias latipes', 'Gasterosteus aculeatus', 'Danio rerio', 'Ciona savignyi', 'Ciona intestinalis', 'Odontophorus guttatus', 'Ambystoma mexicanum', 'Petromyzon marinus');

my %Hensembl_scientific_names;

@Hensembl_scientific_names{@Aensembl_scientific_names}=();

my $ensembl_scientific_names;

subtype ScientificSpeciesName,
	as Str,
	where { ( $_ eq 'ArabidopsisTAIR10' ) || ( $_ eq 'Niben101' ) || ( $_ eq 'Arabidopsis thaliana' ) || ( $_ eq 'Arabidopsis lyrata' ) || ( $_ eq 'Oryza sativa japonica' ) || ( $_ eq 'Oryza sativa Indica Group' ) || ( $_ eq 'Oryza glaberrima' ) || ( $_ eq 'Populus trichocarpa' ) || ( $_ eq 'Vitis vinifera' ) || ( $_ eq 'Zea mays' ) || ( $_ eq 'Sorghum bicolor' ) || ( $_ eq 'Homo sapiens' ) || ( $_ eq 'Medicago truncatula' ) || ( $_ eq 'Glycine max' ) || ( $_ eq 'Ricinus communis' ) || ( $_ eq 'Carica papaya' ) || ( $_ eq 'Lycopersicon esculentum' ) || ( $_ eq 'Gossypium hirsutum' ) || (exists $Hensembl_scientific_names{$_}) },
	message { "$_ invalid ScientificSpeciesName" };

my @Aensembl_species = ('human', 'chimpanzee', 'gorilla', 'orangutan', 'macaque', 'baboon', 'colobus', 'tarier', 'night monkey', 'marmoset', 'dusky titi monkey', 'mouse lemur', 'bushbaby', 'tree shrew', 'rabbit', 'pika', 'squirrel', 'mouse', 'rat', 'kangaroo rat', 'guinea pig', 'cow', 'dolphin', 'alpaca', 'pig', 'dog', 'cat', 'horse', 'microbat', 'megabat', 'hedgehog', 'shrew', 'hyrax', 'elephant', 'lesser hedgehog', 'sloth', 'armadillo', 'wallaby', 'opossum', 'platypus', 'chicken', 'zebrafinch', 'lizard', 'xenopus', 'tetraodon', 'fugu', 'medaka', 'stickleback', 'zebrafish', 'savignyi', 'intestinalis', 'quail', 'salamander', 'lamprey', 'agaricus', 'magnaporthe', 'sclerotinia', 'botrytis');

my %Hensembl_species;

@Hensembl_species{@Aensembl_species}=();

my $ensembl_species;

subtype APPLESSpeciesName,
	as Str,
	where { ( $_ eq 'ArabidopsisTAIR10' ) || ( $_ eq 'Niben101' ) || ( $_ eq 'arabidopsis' ) || ( $_ eq 'lyrata' ) || ( $_ eq 'poplar' ) || ( $_ eq 'grape' ) || ( $_ eq 'rice' ) || ( $_ eq 'maize' ) || ( $_ eq 'sorghum' ) || ( $_ eq 'rice' ) || ( $_ eq 'indica' ) ||  ( $_ eq 'glaberrima' ) || ( $_ eq 'brassica' ) || ( $_ eq 'basidio' ) || ( $_ eq 'medicago' )  || ( $_ eq 'soybean' ) || ( $_ eq 'papaya' ) || ( $_ eq 'tomato' ) || ( $_ eq 'cotton' ) || ( $_ eq 'castor' ) || (exists $Hensembl_species{$_}) || ( $_ eq 'vitis_vinifera' ) || ( $_ eq 'swamp she-oak' ) || $_ eq "nasonia vitripennis" || $_ eq "apis mellifera" || $_ eq "drosophila melanogaster" || $_ eq "atta cephalotes" || $_ eq "bombyx mori" || $_ eq "homo sapiens" || $_ eq "oryza sativa" || $_ eq "vitis vinifera" || $_ eq "musa acuminata" || $_ eq "arabidopsis thaliana" || $_ eq "sorghum bicolor"|| $_ eq "populus trichocarpa" || $_ eq "medicago truncatula" || $_ eq "oryza sativa" || $_ eq "vitis vinifera" || $_ eq "phoenix dactylifera" ||  $_ eq ""},
	message { "$_ invalid APPLESSpeciesName" };
	
subtype RegLocLocusType,
	as Str,
	where { ( $_ eq 'gene' ) || ( $_ eq 'mirna') || ($_ eq 'chip_fragment') },
	message { "$_ invalid RegLocLocusType" };

subtype HomologyType,
	as Str,
	where { ( $_ eq 'ortholog' ) || ( $_ eq 'paralog' ) || ( $_ eq 'homeolog' )
	|| ( $_ eq 'homolog' ) },
	message { "$_ invalid HomologyType" };

subtype GenomicDirection,
	as Str,
	where { ( $_ eq 'towards_five_prime' ) || ( $_ eq 'towards_three_prime' ) },
	message { "$_ invalid GenomicDirection" };

subtype BiFaAlgorithm,
	as Str,
	where { ( $_ eq 'ALG_OTT' ) || ( $_ eq 'ALG_BAYESIAN' )},
	message { "$_ invalid BiFaAlgorithm" };

subtype SequenceChemicalType,
	as Str,
	where { ( $_ eq 'dna' ) || ( $_ eq 'protein' ) || ( $_ eq 'rna' ) },
	message { "$_ invalid SequenceChemicalType" };
	
subtype RunningMode,
	as Str,
	where { ($_ eq 'statistics') || ( $_ eq 'statistics_and_retrieval' ) || ($_ eq 'preparation') || ($_ eq 'normal') },
	message { "$_ invalid RunningMode" };

subtype JobHandlerIdentifier,
	as Str,
	where { ( $_ eq 'seaweed' ) || ( $_ eq 'APPLES_function' ) },
	message { "$_ invalid JobHandlerIdentifier" };

subtype AlignmentAlgorithm,
	as Str,
	where { ( $_ eq 'seaweed' ) || ( $_ eq 'ott' ) },
	message { "$_ invalid AlignmentAlgorithm" };

subtype Status, 
	as Str,
	where { ( $_ eq 'queued' ) || ( $_ eq 'running' ) || ( $_ eq 'complete' ) || ( $_ eq 'failed' ) || ( $_ eq 'never' ) },
	message { "$_ invalid Status" };

subtype WebService,
	as Str,
	where { ( $_ eq 'meme' ) },
	message { "$_ invalid WebService" };

subtype EvolutionaryTreeKingdom,
        as Str,
        where { ( $_ eq 'plants' ) || ( $_ eq 'fungi' ) || ( $_ eq 'vertebrates' ) || ( $_ eq 'all' ) || ( $_ eq 'invertebrates' ) },
        message { "$_ invalid EvolutionaryTreeTypes" };

subtype WeightMatrixSource,
	as Str,
	where { ($_ eq 'BiFa') || ($_ eq 'custom') },
	message { "$_ invalid WeightMatrixSource" };

subtype SequenceRegionType,
    as Str,
    where { ( $_ eq 'upstream' ) || ( $_ eq 'downstream' ) || ( $_ eq 'generegion' ) || ( $_ eq 'upstream_and_generegion' ) || ( $_ eq 'downstream_and_generegion' ) || ( $_ eq 'all' ) },
    message { "$_ invalid SequenceRegionType" };

subtype PrintLevel,
  as Int,
  where { ( $_ == 0 ) || ( $_ == 1 ) || ( $_ == 2 ) || ( $_ == 3 ) },
  message { "$_ invalid PrintLevel" };

subtype ScoringModel,
	as Str,
	where { ( $_ eq 'BiFa' ) || ( $_ eq 'empirical' ) },
	message { "$_ invalid ScoringModel" };

subtype JobMode,
  as Str,
  where { ( $_ eq 'session' ) || ( $_ eq 'queue' ) },
  message { "$_ invalid JobMode" };

subtype ResultType,
  as Str,
  where { ($_ eq 'applesthing' ) }, # all Class names, or file?
  message { "$_ invalid ResultType" };

subtype ListenerType,
  as Str,
  where { ($_ eq 'request') || ($_ eq 'running') || ($_ eq 'complete') || ($_ eq 'error') || ($_ eq 'trash') },
  message { "$_ invalid ListenerType" };

subtype MethodType,
  as Str,
  where { ( $_ eq 'APPLES_Ott_Job' ) || ( $_ eq 'APPLES_Seaweed_Job' ) || ( $_ eq 'pick_random_gene_ids') },
  message { "$_ invalid MethodType" };

subtype CommandLine,
  as Str,
  where { ( $_ eq 'perl perlscript.pl' ) || ( $_ eq 'APPLES_Ott_Job' ) || ( $_ eq 'APPLES_Seaweed_Job' ) || ( $_ eq 'APPLES_MEME_Job' ) },
  message { "$_ invalid CommandLine" };

subtype EnsemblLocation, 
  as Str, 
  where { ( $_ eq 'ensembl' ) || ( $_ eq 'ensemblgenomes' ) || ( $_ eq 'local' ) },
  message { "$_ invalid EnsemblLocation" };

subtype StrandChoice,
  as Str,
  where { ( $_ eq 'positive' ) || ( $_ eq 'negative' ) || ($_ eq 'both') },
  message { "$_ invalid StrandChoice" };

subtype SequenceFunctionalType,
  as Str,
  where { ( $_ eq 'cdna' ) || ( $_ eq 'genomic_dna' ) || ( $_ eq 'micro_rna' ) || ( $_ eq 'protein' ) },
  message { "$_ invalid SequenceFunctionalType" };

subtype TranscriptChoice,
  as Str,
  where { ( $_ eq 'one' ) || ( $_ eq 'all'  ) },
  message { "$_ invalid TranscriptChoice" };

subtype GenomicSequenceClass,
	as Str,
	where { ( $_ eq 'intergenic' ) || ( $_ eq '500bp_promoter' ) || ( $_ eq '1000bp_promoter' ) || ( $_ eq 'coding' ) || ( $_ eq 'whole_genome' )},
	message { "$_ invalid ScoringModel" };

subtype PerseveranceCategory,
  as Str,
  where { ( $_ eq 'Ensembl' ) || ( $_ eq 'BiFa' ) || ( $_ eq 'Unknown') || ( $_ eq 'TempDirs') || ( $_ eq 'CacheDB') },
  message { "$_ invalid PerseveranceCategory" };

subtype PositiveNum,
as Num,
where { ($_ > 0) },
message { "The number you provided, $_, was not a positive number" };

subtype 'NonNegativeNum', # used quotes here in order to be able to specify an attribute in the class Exception 
                          # as 'HashRef[NonNegativeNum]', otherwise this would not work neither with or without
                          # quotes around this attribute type
as Num,
where { ($_ >= 0) },
message { "The number you provided, $_, was not a non-negative number" };

subtype PseudocountType,
	as Str,
	where { ( $_ eq 'p_over_n' ) || ( $_ eq 'p_over_sqrt_n' ) || ( $_ eq 'p' ) },
	message { "$_ empirical scoring metric" };

subtype PositiveInt,
  as Int,
  where { ($_ > 0) },
  message { "The number you provided, $_, was not a positive integer" };

subtype 'NonNegativeInt', # reason for quotes see above
                          # it may be best to always use quotes
    as Int,
    where { ($_ >= 0) },
    message { "The number you provided, $_, was not a non-negative integer" };

subtype Probability,
    as Num,
    where { ($_ >= 0.0) && ($_ <= 1.0) },
    message { "The number you provided, $_, was not a probability"};

subtype ProteinProteinNetworks,
    as Str,
    where { ($_ eq 'human_intact') || ($_ eq 'no_network') }, # choose 'no_network' if you do not have a protein-protein network
    message { "Chosen protein-protein network is not supported."};

subtype WorkingSequenceType,
    as Str,
    where { ($_ eq 'ref_sequence') || ($_ eq 'ref_sequence_repeat_masked') || ($_ eq 'modified_sequence') },
    message { "The sequence type you have specified is an invalid type"};

subtype IntPercentage,
    as Int,
    where { ($_ >= 0) && ($_ <= 100) },
    message { "The number you provided is not an integer percentage." };

1;

