### (c) copyright University of Warwick 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Ortholog_Mapping_Maker Class ###
# Maker class for Ortholog_Mapping objects
use MooseX::Declare;

class Ortholog_Mapping_Maker {
  use Ortholog_Mapping;
  use Alignment_Utilities;
  use Parameters;
  use Genome_DB_Utilities;
  use Data::Dumper;
  use APPLES_Datatypes qw (Boolean);
  use constant {FALSE => 0,
		  TRUE	=> 1};
  my $gdbu = Genome_DB_Utilities->new();
  my $GU = General_Utilities->new();

  method create_ortholog_mappings (Genome_Sequence_Database_Parameters $genome_a, Genome_Sequence_Database_Parameters $genome_b) {
    my $ncbi_blast_parameters = NCBI_BLAST_Parameters->new(
							   word_length => 0,
							   alignments_to_do => 1,
							   alignments_to_show => 1,
							   program => 'blastp',
							   strands => 'positive'
							  );

    my $max_matches = 1;
    my $transcript_choice = 'one';
    my $AU = Alignment_Utilities->new();

    my $first_alignment_set = $AU->align_protein_sequences_in_genome_DBs($genome_a, $genome_b, $ncbi_blast_parameters, $max_matches, $transcript_choice);
    $GU->user_info(2, "\nfirst set derived\n");
    my $second_alignment_set = $AU->align_protein_sequences_in_genome_DBs($genome_b, $genome_a, $ncbi_blast_parameters, $max_matches, $transcript_choice);

    # Via AU/GDBU: <<takes 2 genomes, of any type or source,
    # derives a set of proteins from each,
    # calculates reciprocal best hits using BLAST (very basic approach to ortholog assignment) >>
    
    # get gene id lists    
    my @all_ids_A = $gdbu->list_stable_ids_for_a_given_genome_database($genome_a);
    my @all_ids_B = $gdbu->list_stable_ids_for_a_given_genome_database($genome_b);

    # derive reciprocal best hit hashes from $alignment_set_1 and $alignment_set_2
    my %A_reciprocal_B = $self->private_get_reciprocal_best_hits($first_alignment_set, $second_alignment_set, \@all_ids_A);
    #$GU->user_info(3,"first hashed derived\n");
    my %B_reciprocal_A = $self->private_get_reciprocal_best_hits($second_alignment_set, $first_alignment_set, \@all_ids_B);
    

    # Then, uses alignment set pair to produce an Ortholog_Mapping object
    my $ortholog_mapping = Ortholog_Mapping->new(genome_a => $genome_a,
						 genome_b => $genome_b,
						 a_to_b => \%A_reciprocal_B,
						 b_to_a => \%B_reciprocal_A
						);
    return $ortholog_mapping;
  } # create_ortholog_mappings #
  
  method create_ortholog_mappings_through_job_handler(Genome_Sequence_Database_Parameters $genome_a, Genome_Sequence_Database_Parameters $genome_b) {
    my $ortholog_mapping_maker = Ortholog_Mapping_Maker->new();
    my $cache = TRUE;
    my $cache_duration = 360;
    my $job_handler = Job_Handler->new();
    $job_handler->get_config_settings();
    my $function = 'create_ortholog_mappings';
    
    # Note: We only want to generate keyfiles containing the minimal information, so no superfluous attributes should be passed.
    # Make 'minimal info' copies of the Genome_Sequence_Database_Parameters objects
    
    my $minimal_info_genome_a = $self->private_make_minimal_info_GSDP_object($genome_a);
    my $minimal_info_genome_b = $self->private_make_minimal_info_GSDP_object($genome_b);

    push (my @parameters, $minimal_info_genome_a, $minimal_info_genome_b);
    my $job_parameters = Job_Parameters->new(memory_requirement_high => TRUE,
					     wall_time_estimate => 172800
					    );
    my $aggregate_exception = Aggregate_Exception->new();
    my @cache_result = eval {
      $job_handler->handle_APPLES_function($function, $ortholog_mapping_maker, \@parameters, $cache, $cache_duration, $job_parameters);
    };
    if ($@) {
      my $exception_content = $@;
      my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
      # no parallelism here, though
      $aggregate_exception = $outcome[1];
      if ($outcome[0]) {
	$aggregate_exception->print_statistics_and_die;
      } else {
	die 'did not expect this could happen at design time of this code.';
      }
    }
    my $ortholog_mapping_object = $cache_result[1];
    return $ortholog_mapping_object;
  } # create_ortholog_mappings_through_job_handler #

  method private_get_reciprocal_best_hits(Alignment_Set $alignment_set_1, Alignment_Set $alignment_set_2, ArrayRef $ids_ref) {
    
    my @ids = @{$ids_ref};
    #$GU->user_info(3,Dumper (@ids));
    my %HA_reciprocal_B = (); # initialise hash
    foreach my $id(@ids) {
      $HA_reciprocal_B{$id} = 'not found';
    } # fill keys and set values to 'not found'
    
    foreach my $id(@ids) {
      # look up key $id in alignment_set
      my @hits_in_B = $alignment_set_1->give_hits($id);
      my $top_hit_in_B = $hits_in_B[0]; # take first element of array (there's only going to be one anyway)
      $GU->user_info (3, "top hit of ". $id. " is ". $top_hit_in_B. "\n");

      # look up that result in opposite alignment set
      my @hits_in_A =  $alignment_set_2->give_hits($top_hit_in_B);
      my $top_hit_in_A = $hits_in_A[0];
      $GU->user_info(3,  "top hit of ". $top_hit_in_B. " is ". $top_hit_in_A. "\n");
      #set value in %HA_reciprocal_B to value if they are the same
      if ($id eq $top_hit_in_A) {
	$GU->user_info (3, "$id and $top_hit_in_B are reciprocal best hits\n");
	$HA_reciprocal_B{$id} = $top_hit_in_B;
      }
      else {
	$GU->user_info (3, "no reciprocal hit found for ". $id."\n");
      }
    }
    
    #$GU->user_info(3,"HASH:\n");
    while ( my ($key_query, $value_target) = each(%HA_reciprocal_B) ) {
      #$GU->user_info(3,"$key_query => $value_target\n");
    }
    return %HA_reciprocal_B;
  } # private_get_reciprocal_best_hits #

  method private_make_minimal_info_GSDP_object (Genome_Sequence_Database_Parameters $genome) {
    # make new 'minimal info' object, containing only those attributes that will affect the result of the computation

    my $minimal_info_genome;
    
    my $genome_dbname = $genome->dbname; # required by all

    # WE NEED TO TELL WHAT KIND OF DATABASE_PARAMETERS CLASS THE OBJECTS ARE, AND MAKE ONE OF THOSE
    if ($genome->isa('Ensembl_Database_Parameters')) {
      my $genome_alias = $genome->alias; # required by ensembl
      my $genome_location = $genome->location; # required by ensembl
      $minimal_info_genome = Ensembl_Database_Parameters->new( alias=>$genome_alias, 
							       dbname=>$genome_dbname, 
							       location=>$genome_location );  
    }
    elsif ($genome->isa('Genbank_Sequence_Database_Parameters')) {
      my $block_sort_filename_trigger = $genome->block_sort_filename_trigger; # required by genbank
      my $filenames = $genome->filenames; # required by genbank
      $minimal_info_genome = Genbank_Sequence_Database_Parameters->new( dbname=>$genome_dbname, 
									filename=>$filenames,
									block_sort_filename_trigger=>$block_sort_filename_trigger
								      );  
    }
    elsif ($genome->isa('FASTA_Sequence_Database_Parameters')) {
      my $genome_filename = $genome->filename; # required by fasta
      $minimal_info_genome = FASTA_Sequence_Database_Parameters->new( dbname=>$genome_dbname, 
								      filename=>$genome_filename );  
    }
    else {
      die 'function private_make_minimal_info_GSDP_object cannot handle this type of Genome_Sequence_Database_Parameters object yet';
    }
    return $minimal_info_genome;
  }


} # Ortholog_Mapping_Maker #
