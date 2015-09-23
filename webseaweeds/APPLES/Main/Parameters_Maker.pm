### (c) copyright University of Warwick 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Parameters_Maker Class ###
use MooseX::Declare;

class Parameters_Maker {
    use Parameters;
    use Genome_DB_Utilities;
    use General_Utilities;
    use APPLES_Datatypes qw (Boolean PositiveInt APPLESSpeciesName SequenceRegionType Probability);
    use Data::Dumper;
    
    use constant {FALSE => 0,
		  TRUE	=> 1};
    my $GDBU = Genome_DB_Utilities->new();
    my $GU = General_Utilities->new();

    method default_NCBI_BLAST_parameters_for_DNA () {
	my $result = NCBI_BLAST_Parameters->new();
	$result->word_length(7);
	$result->program('blastn');
	$result->strands('both');
	return $result;
    } # default_NCBI_BLAST_parameters_for_DNA #

    method default_NCBI_BLAST_parameters_for_protein () {
	my $result = NCBI_BLAST_Parameters->new();
	$result->word_length(2);
	$result->program('blastp');
	$result->strands('positive');
	return $result;
    } # default_NCBI_BLAST_parameters_for_protein #


	# Method for constructing default ReRe set parameters excluding core promoters, around source database as a 
	# core species and target database(s). The phylogenetic tree will also be pruned to the specified target species.
	# Input:    $source_db_parameters - Genome_Sequence_Database_Parameters for the source database.
	#			@target_db_parameters - Genome_Sequence_Database_Parameters for the target database(s).
	#			$max_length_for_phylogenetic_search - PositiveInt for length of the phylogenetic search.
	#			$stop_at_neighbouring_gene_for_phylogenetic_search - Boolean whether to stop at neighbouring genes.
	#			$belief_score - Num a belief score to use.
	#			$ortholog_method - Str that specifies what orthology method to use.
	#			$kingdom - Kingdom of interest to use in tree pruning and Compara orthology if specified.
	#			$region - Str that specifies the direction of the search: 'upstream' or 'downstream'.
	#			$alignment_algorithm - Window_Pair_Algorithm_Parameters for the alignment algorithm to use.
	method default_rere_set_constructor_parameters(Genome_Sequence_Database_Parameters $source_db_parameters,
									ArrayRef[Genome_Sequence_Database_Parameters] $target_db_parameters,
									PositiveInt $max_length_for_phylogenetic_search, 
									Boolean $stop_at_neighbouring_gene_for_phylogenetic_search, 
									Probability $belief_score,
									Str $ortholog_method,
									Str $kingdom,
									SequenceRegionType $direction,
									Window_Pair_Algorithm_Parameters $alignment_algorithm) {

        # Establish orthology method(s)
        my $orthology_method;
        
        # Terrible way of doing this, but switch doesn't work with MooseX (5.10.0) and elseif don't work too.
        if ($ortholog_method eq "RBH") {
		print "RBH\n";
        	$orthology_method = RBH_Orthology_Finding_Parameters->new();
        } else {
        	if ($ortholog_method eq "Random") {
			print "random\n";
        		$orthology_method = Random_Assignment_Orthology_Finding_Parameters->new();
        	} else {
			print "compara\n";
        		if ($ortholog_method eq "Compara") {
        			$orthology_method = Compara_Orthology_Finding_Parameters->new(compara_db => $kingdom);
        		} else {
				print "OOPS none\n";
			}
        	}
        }


 		# Create an array of database parameters
		push(my @genome_sequence_database_parameters, $source_db_parameters);
		my @selected_species;
		push(@selected_species, $source_db_parameters->alias);
		
		my @allowed_source_dbnames;
		push(@allowed_source_dbnames, $source_db_parameters->dbname);
		
		for my $target_db(@$target_db_parameters) {
			push(@genome_sequence_database_parameters, $target_db);
			push(@allowed_source_dbnames, $target_db->dbname);
			push(@selected_species, $target_db->alias);
		}


                my $orthology_method_and_species_restriction = Orthology_Method_And_Species_Restriction->new(generic_orthology_finding_parameters => $orthology_method, allowed_source_dbnames => \@allowed_source_dbnames);

		my @orthology_methods;
		push (@orthology_methods, $orthology_method_and_species_restriction);
		
		for my $target_db(@$target_db_parameters) {
			$target_db->{orthology_methods} = \@orthology_methods;
		}
								
		#my $orthology_method_and_species_restriction = Orthology_Method_And_Species_Restriction->new(generic_orthology_finding_parameters => $orthology_method, allowed_source_dbnames => \@allowed_source_dbnames);
		
		#print "ORTHOLOGY\n";
		#print Dumper($orthology_method_and_species_restriction);
		
		#my @orthology_methods;
		#push (@orthology_methods, $orthology_method_and_species_restriction);

		my $tree = Evolutionary_Tree_Maker->new()->make_evolutionary_tree($kingdom);
		my $pruned_tree = $tree->prune_tree(\@selected_species);
                
        # need to make a PTM
		my $PTM = Partial_Threshold_Matrix_Maker->new()->make_partial_threshold_matrix($pruned_tree, $kingdom);	

        # needs bundler parameters
		my $star_bundler_parameters = Star_Bundler_Parameters->new(partial_threshold_matrix => $PTM,
					belief_value => $belief_score); # default settings
		
		my $sequence_parameters = Sequence_Parameters->new( region => $direction,
					stop_at_neighbouring_gene => $stop_at_neighbouring_gene_for_phylogenetic_search,
					max_length_to_search => $max_length_for_phylogenetic_search); # + default settings
	
			
		my $remo_set_phylogenetic_constructor_parameters =
	    		ReMo_Set_Phylogenetic_Constructor_Parameters->new(sequence_databases_to_use_for_homologs => \@genome_sequence_database_parameters,
					window_pair_algorithm_parameters => $alignment_algorithm,
					star_bundler_parameters => $star_bundler_parameters, 
					sequence_parameters => $sequence_parameters);
					
		# set up for Ensembl API 55
		my @remo_set_constructor_parameters;
					
	   	push (@remo_set_constructor_parameters, $remo_set_phylogenetic_constructor_parameters); 


		my $rere_constructor_parameters = ReRe_Constructor_Parameters->new(remo_set_constructor_parameters=>
					\@remo_set_constructor_parameters);
		my $rere_set_constructor_parameters = ReRe_Set_Constructor_Parameters->new(rere_constructor_parameters=>$rere_constructor_parameters);
		
		return $rere_set_constructor_parameters;
    } # default_no_core_rere_set_constructor_parameters_arabidopsis #
    

    method default_rere_set_constructor_parameters_arabidopsis (PositiveInt $window_length, Boolean $stop_at_neighbouring_gene_for_core_promoter, PositiveInt $core_promoter_max_length, Boolean $use_phylogenetic_ReMos, PositiveInt $max_length_for_phylogenetic_search, Boolean $stop_at_neighbouring_gene_for_phylogenetic_search, Num $belief_score, PositiveInt $core_promoter_min_length) {
        # set up for Ensembl API 55

	my @remo_set_constructor_parameters;

        ### constructor for core promoters ###
	my $remo_core_promoter_cons_params = ReMo_Core_Promoter_Constructor_Parameters->new(length => $core_promoter_max_length,
											    stop_at_neighbouring_gene => $stop_at_neighbouring_gene_for_core_promoter,
											    minimum_length => $core_promoter_min_length);
	my @homolog_databases = ();
	my $remo_set_core_promoter_cons_params = ReMo_Set_Core_Promoter_Constructor_Parameters->new(remo_core_promoter_constructor_parameters => $remo_core_promoter_cons_params,
												    sequence_databases_to_use_for_homologs => \@homolog_databases);
	push(@remo_set_constructor_parameters,$remo_set_core_promoter_cons_params);

        ### constructor for phylogenetically conserved regions ###
	# Establish database parameters for the core species (the source)
	my $genome_sequence_database_parameters = Ensembl_Database_Parameters->new( alias => 'arabidopsis', location => 'ensemblgenomes', dbname => 'arabidopsis_thaliana_core_3_55_9'); 

        # Create an array of database parameters
	push(my @genome_sequence_database_parameters, $genome_sequence_database_parameters);

        # Establish orthology method(s)
	my $orthology_method = Compara_Orthology_Finding_Parameters->new(compara_db => 'plants');

	#my @allowed_source_species = ('arabidopsis', 'poplar', 'grape');
	my @allowed_source_dbnames = ('arabidopsis_thaliana_core_3_55_9', 
								'vitis_vinifera_core_3_55_1', 
								'populus_trichocarpa_core_3_55_11');
	my $orthology_method_and_species_restriction = Orthology_Method_And_Species_Restriction->new(generic_orthology_finding_parameters => $orthology_method, allowed_source_dbnames => \@allowed_source_dbnames);
	my @orthology_methods;
	push (@orthology_methods, $orthology_method_and_species_restriction);

        # Establish database parameters for the homolog species (target species)
	my $poplarEDP = Ensembl_Database_Parameters->new( alias => 'poplar', location => 'ensemblgenomes', dbname => 'populus_trichocarpa_core_3_55_11', orthology_methods => \@orthology_methods); 
	push(@genome_sequence_database_parameters, $poplarEDP);

	my $grapeEDP = Ensembl_Database_Parameters->new( alias => 'grape', location => 'ensemblgenomes', dbname => 'vitis_vinifera_core_3_55_1', orthology_methods => \@orthology_methods);
	push(@genome_sequence_database_parameters, $grapeEDP);

        # Establish ReMo_Set_Phylogenetic_Constructor_Parameters
        # needs alignment algorithm parameters
	my $cutoff = int(0.57*$window_length);
	my $window_pair_alignment_parameters = Ott_Algorithm_Parameters->new(windowlength => $window_length,
									     cutoff_for_uninteresting_alignments => $cutoff);
	my $kingdom = 'plants';
	my $tree = Evolutionary_Tree_Maker->new()->make_evolutionary_tree($kingdom);

        # need to make a PTM
	my $PTM = Partial_Threshold_Matrix_Maker->new()->make_partial_threshold_matrix($tree, $kingdom);

        # needs bundler parameters
	my $star_bundler_parameters = Star_Bundler_Parameters->new(partial_threshold_matrix => $PTM,
								   belief_value => $belief_score); # default settings
	my $sequence_parameters = Sequence_Parameters->new( region => 'upstream',
							    stop_at_neighbouring_gene => $stop_at_neighbouring_gene_for_phylogenetic_search,
							    max_length_to_search => $max_length_for_phylogenetic_search); # + default settings

	my $remo_set_phylogenetic_constructor_parameters =
	    ReMo_Set_Phylogenetic_Constructor_Parameters->new(sequence_databases_to_use_for_homologs => \@genome_sequence_database_parameters,
							      window_pair_algorithm_parameters => $window_pair_alignment_parameters,
							      star_bundler_parameters => $star_bundler_parameters, 
							      sequence_parameters => $sequence_parameters);

	if ($use_phylogenetic_ReMos) {
	    push (@remo_set_constructor_parameters, $remo_set_phylogenetic_constructor_parameters); 
	}

	my $rere_constructor_parameters = ReRe_Constructor_Parameters->new(remo_set_constructor_parameters=>
									   \@remo_set_constructor_parameters);
	my $rere_set_constructor_parameters = ReRe_Set_Constructor_Parameters->new(rere_constructor_parameters=>$rere_constructor_parameters);
	return $rere_set_constructor_parameters;
    } # default_rere_set_constructor_parameters_arabidopsis #

    method default_rere_set_constructor_parameters_human (ArrayRef $all_parameters) {
	# put all parameters into one array as a quick work-around

	my @parameters_array = @{$all_parameters};
	my $length = shift @parameters_array;
	my $minimum_length = shift @parameters_array;
	my $stop_at_neighbouring_gene = shift @parameters_array;
	my $ignore_pseudogenes = shift @parameters_array;
	my $ignore_RNAgenes = shift @parameters_array;
	my $dbnames_array_ref = shift @parameters_array;
	my $window_length = shift @parameters_array;
	my $belief_threshold = shift @parameters_array;
	my $max_sequence_length = shift @parameters_array;
	my $min_length_phylogenetic_search = shift @parameters_array; 
	my $sequence_region = shift @parameters_array;
	my $use_seaweed_algorithm = shift @parameters_array;
	
	# set up for Ensembl API 55

	my @dbnames_to_use_for_phylogenetic_conservation = @{$dbnames_array_ref};

	# ---------------------
	# orthology parameters
	# ---------------------
	my $compara = Compara_Orthology_Finding_Parameters->new(compara_db => 'ensembl_compara_55');
	my @allowed_source_dbnames = ('homo_sapiens_core_55_37'); 
	push(@allowed_source_dbnames,@dbnames_to_use_for_phylogenetic_conservation);
	my $one_orthology_method = Orthology_Method_And_Species_Restriction->new(generic_orthology_finding_parameters => $compara,
										 allowed_source_dbnames => \@allowed_source_dbnames);
	my @orthology_methods = ($one_orthology_method);

	# -----------------------------------------
        # pick Ensembl databases for a few species
	# -----------------------------------------
	my $mouse_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'mouse', location => 'ensembl', dbname => 'mus_musculus_core_55_37h', orthology_methods => \@orthology_methods);
	my $opossum_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'opossum', location => 'ensembl', dbname => 'monodelphis_domestica_core_55_5i', orthology_methods => \@orthology_methods);
	my $xenopus_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'xenopus', location => 'ensembl', dbname => 'xenopus_tropicalis_core_55_41n', orthology_methods => \@orthology_methods);
	my $tetraodon_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'tetraodon', location => 'ensembl', dbname => 'tetraodon_nigroviridis_core_55_8b', orthology_methods => \@orthology_methods);
	my $fugu_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'fugu', location => 'ensembl', dbname => 'takifugu_rubripes_core_55_4k', orthology_methods => \@orthology_methods);
	my $zebrafish_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'zebrafish', location => 'ensembl', dbname => 'danio_rerio_core_55_8a', orthology_methods => \@orthology_methods);
	my %ensembl_db_hash;
	$ensembl_db_hash{'mus_musculus_core_55_37h'} = $mouse_ensembl_db_params;
	$ensembl_db_hash{'monodelphis_domestica_core_55_5i'} = $opossum_ensembl_db_params;
	$ensembl_db_hash{'xenopus_tropicalis_core_55_41n'} = $xenopus_ensembl_db_params;
	$ensembl_db_hash{'tetraodon_nigroviridis_core_55_8b'} = $tetraodon_ensembl_db_params;
	$ensembl_db_hash{'takifugu_rubripes_core_55_4k'} = $fugu_ensembl_db_params;
	$ensembl_db_hash{'danio_rerio_core_55_8a'} = $zebrafish_ensembl_db_params;

	# ----------------------------
	# check user's species choice
	# ----------------------------
	my @homolog_databases_for_phylogenetic_remos;
	foreach my $chosen_species (@dbnames_to_use_for_phylogenetic_conservation) {
	    my $hash_entry = $ensembl_db_hash{$chosen_species};
	    if ($hash_entry eq '') {
		die 'add entry for species '.$chosen_species.' to define a choice of Ensembl-database.';
	    }
	    push(@homolog_databases_for_phylogenetic_remos,$hash_entry);
	}
	
        # ---------------------------------------------------------
        # put together one or more ReMo_Set_Constructor_Parameters
        # ---------------------------------------------------------

        ### constructor for core promoters ###
	my $remo_core_promoter_cons_params = ReMo_Core_Promoter_Constructor_Parameters->new(length => $length,
											    stop_at_neighbouring_gene => $stop_at_neighbouring_gene,
											    ignore_pseudogenes => $ignore_pseudogenes,
											    ignore_RNAgenes => $ignore_RNAgenes,
											    minimum_length => $minimum_length);
	my @homolog_databases = ();
	my $remo_set_core_promoter_cons_params = ReMo_Set_Core_Promoter_Constructor_Parameters->new(remo_core_promoter_constructor_parameters => $remo_core_promoter_cons_params,
												    sequence_databases_to_use_for_homologs => \@homolog_databases);

        ### constructor for phylogenetically conserved regions ###

        # alignment algorithm
	my $cutoff = int(0.57*$window_length);
	my $seaweed_parameters = Seaweed_Algorithm_Parameters->new(windowlength => $window_length,
								   cutoff_for_uninteresting_alignments => $cutoff);
	my $ott_algo_parameters = Ott_Algorithm_Parameters->new(windowlength => $window_length,
								cutoff_for_uninteresting_alignments => $cutoff);
	my $window_pair_parameters;
	if ($use_seaweed_algorithm) {
	    $window_pair_parameters = $seaweed_parameters;
	} else {
	    $window_pair_parameters = $ott_algo_parameters;
	} 

        # bundling
	my $kingdom = 'vertebrates';
	my $evolutionary_tree_maker = Evolutionary_Tree_Maker->new();
	my $tree = $evolutionary_tree_maker->make_evolutionary_tree($kingdom);
	my @allowed_source_species = $GDBU->tell_species_names_from_database_hash(\%ensembl_db_hash);

	my $pruned_tree = $tree->prune_tree(\@allowed_source_species);
	my $partial_threshold_matrix_maker = Partial_Threshold_Matrix_Maker->new();
	my $partial_threshold_matrix = $partial_threshold_matrix_maker->make_partial_threshold_matrix($pruned_tree, $kingdom);
	my $star_bundler_params = Star_Bundler_Parameters->new(belief_value => $belief_threshold,
							       partial_threshold_matrix => $partial_threshold_matrix);

        # sequence choice
	my $sequence_params = Sequence_Parameters->new(max_length_to_search => $max_sequence_length,
						       region => $sequence_region,
						       min_length_to_return => $min_length_phylogenetic_search,
						       ignore_RNAgenes => $ignore_RNAgenes,
						       ignore_pseudogenes => $ignore_pseudogenes);

        # put all these together
	my $remo_set_phylogenetic_cons_params = ReMo_Set_Phylogenetic_Constructor_Parameters->new(sequence_databases_to_use_for_homologs => \@homolog_databases_for_phylogenetic_remos,
												  window_pair_algorithm_parameters => $window_pair_parameters,
												  star_bundler_parameters => $star_bundler_params,
												  sequence_parameters => $sequence_params);

        # ------------------------------------------------------
        # add all ReMo_Set_Constructor_Parameters into an array
        # ------------------------------------------------------
	my @remo_set_cons_params = ($remo_set_core_promoter_cons_params,$remo_set_phylogenetic_cons_params);

        # -----------------------------------
        # get ReRe_Set and genomic intervals
        # -----------------------------------
	my $rere_cons_params = ReRe_Constructor_Parameters->new(remo_set_constructor_parameters => \@remo_set_cons_params);
	my $rere_set_const_params = ReRe_Set_Constructor_Parameters->new(rere_constructor_parameters => $rere_cons_params);
	
	return $rere_set_const_params;
	
    } # default_rere_set_constructor_parameters_human #

    method full_default_rere_set_constructor_parameters_human (Boolean $use_phylogenetic_remos) {
	# core promoter
	my $length = 1000;
	my $minimum_length = 100;
	my $stop_at_neighbouring_gene = TRUE;
	my $ignore_pseudogenes = TRUE;
	my $ignore_RNAgenes = TRUE;
	# phylogenetic ReMos
	my @species_to_use_for_phylogenetic_conservation;
	if ($use_phylogenetic_remos) {
	    @species_to_use_for_phylogenetic_conservation = qw (mus_musculus_core_55_37h monodelphis_domestica_core_55_5i xenopus_tropicalis_core_55_41n tetraodon_nigroviridis_core_55_8b  takifugu_rubripes_core_55_4k danio_rerio_core_55_8a);
	} else {
	    @species_to_use_for_phylogenetic_conservation = ();
	}
	my $window_length = 100;
	my $belief_threshold = 0.5;
	my $max_sequence_length = 40000;
	my $min_length_phylogenetic_search = 500;
	my $sequence_region = 'upstream';
	my $use_seaweed_algorithm = FALSE;
	my @parameters_to_pass_on = ($length,$minimum_length,$stop_at_neighbouring_gene,$ignore_pseudogenes,$ignore_RNAgenes,\@species_to_use_for_phylogenetic_conservation,$window_length,$belief_threshold,$max_sequence_length,$min_length_phylogenetic_search,$sequence_region,$use_seaweed_algorithm); 
	my $result = $self->default_rere_set_constructor_parameters_human(\@parameters_to_pass_on);
	return $result;
    } # full_default_rere_set_constructor_parameters_human #

    method default_rere_set_constructor_parameters_ostreococcus (
	PositiveInt $core_promoter_max_length,
	PositiveInt $core_promoter_min_length,
	Boolean $stop_at_neighbouring_gene_for_core_promoter
	) {
	# core promoters only

	my $remo_core_promoter_cons_params = ReMo_Core_Promoter_Constructor_Parameters->new(length => $core_promoter_max_length,
											    stop_at_neighbouring_gene => $stop_at_neighbouring_gene_for_core_promoter,
											    minimum_length => $core_promoter_min_length);
	my @homolog_databases = ();
	my $remo_set_core_promoter_cons_params = ReMo_Set_Core_Promoter_Constructor_Parameters->new(remo_core_promoter_constructor_parameters => $remo_core_promoter_cons_params,
												    sequence_databases_to_use_for_homologs => \@homolog_databases);
	push(my @remo_set_constructor_parameters,$remo_set_core_promoter_cons_params);
	my $rere_constructor_parameters = ReRe_Constructor_Parameters->new(remo_set_constructor_parameters=>
									   \@remo_set_constructor_parameters);
	my $rere_set_constructor_parameters = ReRe_Set_Constructor_Parameters->new(rere_constructor_parameters=>$rere_constructor_parameters);
	return $rere_set_constructor_parameters;
    } # default_rere_set_constructor_parameters_ostreococcus #
    
    method default_model_based_empirical_wms_parameters (Int $seed) {
	my $biobase = Biobase_Additive->new();

	my $model_type = Markov_Model->new();
	my $sequence_model_source = Partial_Sequence_Model_Learning_Params->new(model => $model_type,  learning_restrictions => "whole_genome");
	my $empirical_scoring_params = Model_Based_Empirical_WMS_Parameters->new(
	    sequence_model_source => $sequence_model_source,
	    seed => $seed,
	    scoring_method => $biobase);
	return $empirical_scoring_params;
    } # default_model_based_empirical_wms_parameters #


    method default_rere_set_constructor_parameters_zebrafish (ArrayRef $all_parameters) {
	# put all parameters into one array as a quick work-around

	my @parameters_array = @{$all_parameters};
	my $length = shift @parameters_array;
	my $minimum_length = shift @parameters_array;
	my $stop_at_neighbouring_gene = shift @parameters_array;
	my $ignore_pseudogenes = shift @parameters_array;
	my $ignore_RNAgenes = shift @parameters_array;
	my $dbnames_array_ref = shift @parameters_array;
	my $window_length = shift @parameters_array;
	my $belief_threshold = shift @parameters_array;
	my $max_sequence_length = shift @parameters_array;
	my $min_length_phylogenetic_search = shift @parameters_array; 
	my $sequence_region = shift @parameters_array;
	my $use_seaweed_algorithm = shift @parameters_array;
	
	# set up for Ensembl API 55

	my @dbnames_to_use_for_phylogenetic_conservation = @{$dbnames_array_ref};

	# ---------------------
	# orthology parameters
	# ---------------------
	my $compara = Compara_Orthology_Finding_Parameters->new(compara_db => 'ensembl_compara_55');
	my @allowed_source_dbnames = ('danio_rerio_core_55_8a'); 
	push(@allowed_source_dbnames,@dbnames_to_use_for_phylogenetic_conservation);
	my $one_orthology_method = Orthology_Method_And_Species_Restriction->new(generic_orthology_finding_parameters => $compara,
										 allowed_source_dbnames => \@allowed_source_dbnames);
	my @orthology_methods = ($one_orthology_method);

	# -----------------------------------------
        # pick Ensembl databases for a few species
	# -----------------------------------------
	my $mouse_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'mouse', location => 'ensembl', dbname => 'mus_musculus_core_55_37h', orthology_methods => \@orthology_methods);
	my $opossum_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'opossum', location => 'ensembl', dbname => 'monodelphis_domestica_core_55_5i', orthology_methods => \@orthology_methods);
	my $xenopus_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'xenopus', location => 'ensembl', dbname => 'xenopus_tropicalis_core_55_41n', orthology_methods => \@orthology_methods);
	my $tetraodon_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'tetraodon', location => 'ensembl', dbname => 'tetraodon_nigroviridis_core_55_8b', orthology_methods => \@orthology_methods);
	my $fugu_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'fugu', location => 'ensembl', dbname => 'takifugu_rubripes_core_55_4k', orthology_methods => \@orthology_methods);
	my $human_ensembl_db_params = Ensembl_Database_Parameters->new( alias => 'zebrafish', location => 'ensembl', dbname => 'homo_sapiens_core_55_37', orthology_methods => \@orthology_methods);
	my %ensembl_db_hash;
	$ensembl_db_hash{'mus_musculus_core_55_37h'} = $mouse_ensembl_db_params;
	$ensembl_db_hash{'monodelphis_domestica_core_55_5i'} = $opossum_ensembl_db_params;
	$ensembl_db_hash{'xenopus_tropicalis_core_55_41n'} = $xenopus_ensembl_db_params;
	$ensembl_db_hash{'tetraodon_nigroviridis_core_55_8b'} = $tetraodon_ensembl_db_params;
	$ensembl_db_hash{'takifugu_rubripes_core_55_4k'} = $fugu_ensembl_db_params;
	$ensembl_db_hash{'homo_sapiens_core_55_37'} = $human_ensembl_db_params;

	# ----------------------------
	# check user's species choice
	# ----------------------------
	my @homolog_databases_for_phylogenetic_remos;
	foreach my $chosen_species (@dbnames_to_use_for_phylogenetic_conservation) {
	    my $hash_entry = $ensembl_db_hash{$chosen_species};
	    if ($hash_entry eq '') {
		die 'add entry for species '.$chosen_species.' to define a choice of Ensembl-database.';
	    }
	    push(@homolog_databases_for_phylogenetic_remos,$hash_entry);
	}
	
        # ---------------------------------------------------------
        # put together one or more ReMo_Set_Constructor_Parameters
        # ---------------------------------------------------------

        ### constructor for core promoters ###
	my $remo_core_promoter_cons_params = ReMo_Core_Promoter_Constructor_Parameters->new(length => $length,
											    stop_at_neighbouring_gene => $stop_at_neighbouring_gene,
											    ignore_pseudogenes => $ignore_pseudogenes,
											    ignore_RNAgenes => $ignore_RNAgenes,
											    minimum_length => $minimum_length);
	my @homolog_databases = ();
	my $remo_set_core_promoter_cons_params = ReMo_Set_Core_Promoter_Constructor_Parameters->new(remo_core_promoter_constructor_parameters => $remo_core_promoter_cons_params,
												    sequence_databases_to_use_for_homologs => \@homolog_databases);

        ### constructor for phylogenetically conserved regions ###

        # alignment algorithm
	my $cutoff = int(0.57*$window_length);
	my $seaweed_parameters = Seaweed_Algorithm_Parameters->new(windowlength => $window_length,
								   cutoff_for_uninteresting_alignments => $cutoff);
	my $ott_algo_parameters = Ott_Algorithm_Parameters->new(windowlength => $window_length,
								cutoff_for_uninteresting_alignments => $cutoff);
	my $window_pair_parameters;
	if ($use_seaweed_algorithm) {
	    $window_pair_parameters = $seaweed_parameters;
	} else {
	    $window_pair_parameters = $ott_algo_parameters;
	} 

        # bundling
	my $kingdom = 'vertebrates';
	my $evolutionary_tree_maker = Evolutionary_Tree_Maker->new();
	my $tree = $evolutionary_tree_maker->make_evolutionary_tree($kingdom);
	my @allowed_source_species = $GDBU->tell_species_names_from_database_hash(\%ensembl_db_hash);

	my $pruned_tree = $tree->prune_tree(\@allowed_source_species);
	my $partial_threshold_matrix_maker = Partial_Threshold_Matrix_Maker->new();
	my $partial_threshold_matrix = $partial_threshold_matrix_maker->make_partial_threshold_matrix($pruned_tree, $kingdom);
	my $star_bundler_params = Star_Bundler_Parameters->new(belief_value => $belief_threshold,
							       partial_threshold_matrix => $partial_threshold_matrix);

        # sequence choice
	my $sequence_params = Sequence_Parameters->new(max_length_to_search => $max_sequence_length,
						       region => $sequence_region,
						       min_length_to_return => $min_length_phylogenetic_search,
						       ignore_RNAgenes => $ignore_RNAgenes,
						       ignore_pseudogenes => $ignore_pseudogenes);

        # put all these together
	my $remo_set_phylogenetic_cons_params = ReMo_Set_Phylogenetic_Constructor_Parameters->new(sequence_databases_to_use_for_homologs => \@homolog_databases_for_phylogenetic_remos,
												  window_pair_algorithm_parameters => $window_pair_parameters,
												  star_bundler_parameters => $star_bundler_params,
												  sequence_parameters => $sequence_params);

        # ------------------------------------------------------
        # add all ReMo_Set_Constructor_Parameters into an array
        # ------------------------------------------------------
	my @remo_set_cons_params = ($remo_set_core_promoter_cons_params,$remo_set_phylogenetic_cons_params);

        # -----------------------------------
        # get ReRe_Set and genomic intervals
        # -----------------------------------
	my $rere_cons_params = ReRe_Constructor_Parameters->new(remo_set_constructor_parameters => \@remo_set_cons_params);
	my $rere_set_const_params = ReRe_Set_Constructor_Parameters->new(rere_constructor_parameters => $rere_cons_params);
	
	return $rere_set_const_params;
	
    } # default_rere_set_constructor_parameters_zebrafish #
    
    method default_reco_parameters () {
	my $p_value_likelihood_conv_params = P_Value_Sampling_Likelihood_Conversion_Parameters->new();
	my $reco_params = ReCo_Parameters->new(
	    p_p_network => 'human_intact',
	    minimum_node_probability => 0.001,
	    number_of_sampling_iterations => 10000,
	    max_number_of_sites => 7,
	    conversion_parameters => $p_value_likelihood_conv_params,
	    scale_number_of_experiments => 1,
	    scale_number_of_experimental_methods => 1,
	    weight_for_two_logistic_functions => 0.5);
	return $reco_params;
    } # default_reco_parameters #

} # Parameters_Maker #
