### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReMo_Set_Maker class ###
# class for constructing groups of ReMo_Sets, for example, by accessing a ReMo-DB
use MooseX::Declare;

class ReMo_Set_Maker {
        use Parameters;
	use ReMo_Set;
	use ReMo_Maker;
	use Star_Bundler;
	use Data::Dumper;
	use Job_Handler;
	use General_Utilities;
	use Genomic_Interval_Set_Maker;
	use Genome_DB_Utilities;
	use constant {FALSE => 0,
		      TRUE	=> 1};
			  	  
	my $GU = General_Utilities->new();
	
	method make_remo_set_from_core_promoter(ArrayRef[Reg_Loc] $all_reg_locs, ReMo_Core_Promoter_Constructor_Parameters $remo_core_promoter_constructor_parameters){
		
		my $remo_maker = ReMo_Maker->new();
		#my $remo_core_promoter_constructor_parameters = ReMo_Core_Promoter_Constructor_Parameters->new(stop_at_neighbouring_gene => FALSE,
		#length => $upstream_length);
		my @remos;
		
		my $count = 0;
		foreach my $reg_loc (@{$all_reg_locs}) {
			my $remo = $remo_maker->make_remo_from_core_promoter($remo_core_promoter_constructor_parameters,$reg_loc,TRUE);
			push(@remos, $remo);
			#$GU->user_info(3,"Got remo\n");
			#$count++;
			#if( $count == 1 ){
					
			#	last;
			#}
		}
		
		my $remo_set = ReMo_Set->new(remo_set => \@remos);
		
		return $remo_set;
	} # make_remo_set_from_core_promoter #
	
	method make_remo_set_from_core_promoter_through_job_handler(ArrayRef[Reg_Loc] $all_reg_locs, ReMo_Core_Promoter_Constructor_Parameters $remo_core_promoter_constructor_parameters) {
	    my $object = ReMo_Set_Maker->new();
	    my $function = 'make_remo_set_from_core_promoter';
	    my @parameters = ($all_reg_locs,$remo_core_promoter_constructor_parameters);
	    my $high_memory = TRUE;
	    my @result = $GU->standard_call_to_job_handler($object,$function,\@parameters,$high_memory,FALSE);
	    return $result[0];
	} # make_remo_set_from_core_promoter_through_job_handler #

	method make_remo_sets_from_core_promoter(ReMo_Set_Core_Promoter_Constructor_Parameters $parameters, Reg_Loc $reg_loc) {
	  my $source_dbname = $reg_loc->genome_db->{dbname};
	  my @remos;
		my $remo_maker = ReMo_Maker->new();
		if ($reg_loc->locus_type ne 'gene') {
		  die 'treatment of non-Ensembl genes not implemented yet'
		}
		my $own_promoter = $remo_maker->make_remo_from_core_promoter($parameters->remo_core_promoter_constructor_parameters, $reg_loc,TRUE);
		push(@remos,$own_promoter);
		
		$GU->user_info( 3, $own_promoter->five_prime_pos.":".$own_promoter->three_prime_pos."\n" );
		$GU->user_info( 3, "own promoter:". $own_promoter->meta->name."\n" );
		
		foreach my $homolog_database (@{$parameters->sequence_databases_to_use_for_homologs}) {
		  if ($homolog_database->{alias} eq $reg_loc->genome_db->{alias}) {
		    $GU->user_info( 3, "skipping self!\n" );
		  }
		  else {
		    # loop through orthology methods here, and push into one list
		    foreach my $orthology_method_and_species_restriction (@{$homolog_database->orthology_methods}) { # deref ArrayRef of orthology methods
			my @remos_one_species = $remo_maker->make_remos_from_orthologous_core_promoters($parameters->remo_core_promoter_constructor_parameters, $homolog_database, \$reg_loc, $orthology_method_and_species_restriction, $source_dbname);
			push(@remos, @remos_one_species);
		    }
		  }
		}		
	  
	  my $remo_set = ReMo_Set->new(remo_set => \@remos);
	  return $remo_set;
	} #  make_remo_sets_from_core_promoter #
	
	method make_phylogenetic_remo_sets (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters,
										Reg_Loc $reg_loc) {
	  if ($reg_loc->locus_type ne 'gene') {
	    die 'treatment of locus types other than genes is not implemented yet';
	  }
	  
	  my $gdbu = Genome_DB_Utilities->new();
	  # get genomic intervals (method of genome_db_utilities) and pass to Star_Bundler
	  my $centre_gi = $gdbu->get_genomic_sequence($reg_loc, $parameters->sequence_parameters);
	  my $gi_set_maker = Genomic_Interval_Set_Maker->new();
	  my $gi_set = $gi_set_maker->make_gi_set_for_phylogenetic_remos($parameters, $reg_loc); 
	  
	  my @remo_sets = Star_Bundler->new()->make_remo_sets($centre_gi, $gi_set, $parameters);
	  return @remo_sets;
	} # make_phylogenetic_remo_sets #

} # ReMo_Set_Maker #
