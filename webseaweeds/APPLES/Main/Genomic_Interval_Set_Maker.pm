### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genomic_Interval_Set_Maker class ###
use MooseX::Declare;

class Genomic_Interval_Set_Maker {
	use Parameters;
	use Genome_DB_Utilities;
	use Genomic_Interval_Set;
	use Data::Dumper;
	use General_Utilities;
	use ReRe;
	use constant {FALSE => 0,
		      TRUE	=> 1};	

	my $GU = General_Utilities->new();
	
	method make_gi_set_from_fasta_file (FASTA_Sequence_Database_Parameters $parameters) {
		
		my $utility = Genome_DB_Utilities->new();
		my @sequences = $utility->retrieve_all_sequences_from_DB($parameters);
		
		my $gi;
		my @gi_array;
		my $length;
		foreach my $sequence(@sequences) {
			
		    # turn into a gi and push into array @gi_set
		    $gi = Genomic_Interval->new(genome_db => $parameters,
						coord_sys_name => $parameters->filename,
						region => $sequence->{NAME},
						five_prime_pos => 1,
						three_prime_pos => length($sequence->{SEQUENCE}),
						strand => 'positive',
			);
		    $gi->gi_sequence($sequence->{SEQUENCE});
		    push (@gi_array, $gi);
		    #$GU->user_info ( 3, Dumper (@gi_array) );
		}
		my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
		return $gi_set; # Genomic_Interval_Set
	} # make_gi_set_from_fasta_file #

method make_gi_set_from_fasta_file_using_specified_ids_and_length (FASTA_Sequence_Database_Parameters $parameters, ArrayRef $ids, Int $length) {
  my @ids = @{$ids};
  my %ids;
  foreach my $item (@ids) { $ids{$item} = 13 };
  my $utility = Genome_DB_Utilities->new();
  my @sequences = $utility->retrieve_all_sequences_from_DB($parameters);
  my $gi;
  my @gi_array;
  # take specified length of sequence, unless this is longer than the database sequence
  foreach my $sequence(@sequences) {
     my $length_to_take = $length;

     my $db_seq_length = length($sequence->{SEQUENCE});
#	print "db seq length = ". $db_seq_length."\n";	
     if ($db_seq_length < $length) {
       $length_to_take = $db_seq_length;
     }

    my $five_prime_pos = $db_seq_length - $length_to_take +1;
 #   print "5 prime pos: ".$five_prime_pos."\n";
    # if id is specified in the list (arrayref->hash), include in gi_set   
    if (exists $ids{$sequence->{NAME}}) {
     	 # turn into a gi and push into array @gi_set
      	$gi = Genomic_Interval->new(genome_db => $parameters,
				  coord_sys_name => $parameters->filename,
				  region => $parameters->dbname,
				  label => $sequence->{NAME},
				  five_prime_pos => $five_prime_pos,
				  three_prime_pos => length($sequence->{SEQUENCE}),
				  strand => 'positive',
				 );
		my $subseq = substr $sequence->{SEQUENCE}, -$length_to_take;
		print "Retrieved gene id ".$sequence->{NAME}." from the FASTA database\n";
      	print $subseq."\n";
      	$gi->gi_sequence($subseq);
      	push (@gi_array, $gi);
      	#$GU->user_info ( 3, Dumper (@gi_array) );
    }
  }
  my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
  return $gi_set; # Genomic_Interval_Set
} # make_gi_set_from_fasta_file_using_specified_ids_and_length #



	method make_gi_set_for_phylogenetic_remos (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters,
						   Reg_Loc $reg_loc) {
	  
	  # utilises GDBU functions to extract sequences for input into ReMo search
	  my $gdbu = Genome_DB_Utilities->new();
	  my @gi_array;
	  my $source_dbname = $reg_loc->genome_db->{dbname};
	  my $aggregate_exception = Aggregate_Exception->new();
	  my $stats_required = FALSE;
	  
	  # also loop through orthology method array here to make separate lists? (LB) then push into one list
	  
	  foreach my $homolog_database_parameters (@{$parameters->sequence_databases_to_use_for_homologs}) {	      
	    my @new_gis = ();
	    eval {
	      if ($homolog_database_parameters->{dbname} ne $reg_loc->genome_db->{dbname}) {
		foreach my $orthology_method_and_species_restriction (@{$homolog_database_parameters->orthology_methods}) { # deref ArrayRef[]
		  # check $source_dbname is in the list of allowed_source_dbnames
		  my @allowed_source_dbnames = @{$orthology_method_and_species_restriction->allowed_source_dbnames};
		  my $exists = FALSE;
		  foreach my $species (@allowed_source_dbnames) {
		    if ($source_dbname eq $species) {
		      $exists = TRUE;
		    }
		  }
		  if (!$exists) {
		    $GU->user_info(3,"Species you allowed:\n".Dumper(\@allowed_source_dbnames));
		    die 'this orthology finding method is not valid for this source database: '.$source_dbname;
		  }
		  
		  # get orthologs, foreach ortholog, get relevant gi
		  my $orthology_method = $orthology_method_and_species_restriction->generic_orthology_finding_parameters;
		  my @orthologous_gene_ids = $gdbu->get_orthologous_gene_IDs( $reg_loc->gene_ID, $reg_loc->genome_db, $homolog_database_parameters, $orthology_method );
		  my $reg_loc_maker = Reg_Loc_Maker->new();
		  my @orthologous_reg_locs = $reg_loc_maker->make_reg_locs_from_list( $homolog_database_parameters, \@orthologous_gene_ids );
		  foreach my $orthologous_reg_loc (@orthologous_reg_locs) {
		    # get relevant gi
		    my $gi = $gdbu->get_genomic_sequence($orthologous_reg_loc, $parameters->sequence_parameters);
		    push @new_gis, $gi;
		  }
		}
	      }
	      else {
		$GU->user_info (3, "skipping self!\n");
	      }
	    };
	    if ($@) {
	      my $exception_content = $@;
	      my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
	      if ($outcome[0]) {
		$stats_required = TRUE;
	      }
	      $aggregate_exception = $outcome[1];
	    }
	    else {
	      push(@gi_array,@new_gis);
	    }
	  }
	  if ($stats_required) {
	    $aggregate_exception->print_statistics_and_die;
	  }
	  my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
	  return $gi_set; # Genomic_Interval_Set
	} # make_gi_set_for_phylogenetic_remos #
	
	method make_genomic_interval_set_from_gene_id_list (ArrayRef $genelistref, Sequence_Parameters $sequence_parameters, Genome_Sequence_Database_Parameters $gsdb_parameters) {
	  #die 'method not implemented yet!';
	  if ($gsdb_parameters->isa('Ensembl_Database_Parameters')) {
	    my @geneIDlist = @{$genelistref};
	    my $gdbu = Genome_DB_Utilities->new();
	    my @gi_array;
	    foreach my $geneid (@geneIDlist) {
	      my $reg_loc = Reg_Loc_Maker->new()->make_reg_loc_from_gene_id ($gsdb_parameters, $geneid);
	      my $gi = $gdbu->get_genomic_sequence($reg_loc, $sequence_parameters);
	      $gi->get_sequence(); # fill out sequences
	      push (@gi_array, $gi);
	      $GU->user_info(3,"adding gi\n");
	    }
	    my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
	    return $gi_set;
	  }
	  elsif ($gsdb_parameters->isa('FASTA_Sequence_Database_Parameters')) {
	    die 'cannot make genomic interval set using FASTA file gene ids';
	  }
	  else {
	    die 'cannot deal with this kind of Genome_Sequence_Database_Parameter';
	  }
	} # make_genomic_interval_set_from_gene_id_list #
	
	method make_gi_set_from_core_promoters () {
	  die 'method not implemented yet!';
	} # make_gi_set_from_core_promoters #

	method make_gi_set_from_gi_array (ArrayRef $gi_array) {
	  my $GI_set = Genomic_Interval_Set->new(genomic_interval_set => $gi_array);
	  return $GI_set;
	} # make_gi_set_from_gi_array #

} # Genomic_Interval_Set_Maker #
