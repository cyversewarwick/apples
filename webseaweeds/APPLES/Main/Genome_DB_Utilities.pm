### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genome_DB_Utilities Class ###
# class abstracting genome databases, most functions will only support Ensembl-database handling, though
use MooseX::Declare;
use lib "/home/grannysmith/ensembl-api/ensembl/modules/";
use lib "/home/grannysmith/ensembl-api/ensembl-compara/modules/";
class Genome_DB_Utilities {
  use APPLES_Datatypes qw (Boolean GenomicDirection APPLESSpeciesName StrandType ScientificSpeciesName SequenceChemicalType TranscriptChoice IntPercentage);
  use constant {FALSE => 0,
		TRUE	=> 1};	
  use Parameters;
  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  use Storable;
  use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
  use Bio::EnsEMBL::Registry;
  use diagnostics;
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Data::Dumper;
  use Bio::Seq;
  use Bio::DB::GenBank;
  use General_Utilities;
  use Genomic_Interval;
  use Genomic_Interval_Set;
  use Reg_Loc;
  use Bio::SeqIO;
  use Exception;
  use File::Copy;
  
  use constant DBVERSION => 3;
  use constant FORCE_REFRESH => FALSE; 

  use Moose;
  has 'registry' => (is => 'rw', isa => 'Any', clearer => 'private_clear_registry');
  has 'current_registry_location' => (is => 'rw', isa => 'Str', clearer => 'private_clear_location'); # is private clearer method necessary?

  my $GU = General_Utilities->new();
  my $gb = new Bio::DB::GenBank;

##########################
##### public methods #####
##########################

  method get_distance_to_end_of_sequence_region(Reg_Loc $reg_loc, GenomicDirection $direction) {
      # sequence region can be a chromosome, scaffold, contig, ...
      # $direction interpreted relative to $reg_loc

      my $result;
      if ($reg_loc->genome_db->isa('Ensembl_Database_Parameters')) {
	  $result = $self->private_get_distance_to_end_of_sequence_region_ensembl($reg_loc,$direction);
      } else {
	  my $type = $reg_loc->genome_db->meta->name;
	  die 'method not yet implemented for this type of genome database: '.$type;
      }
      return $result;
  } # get_distance_to_end_of_sequence_region #

  method check_existence_of_neighbouring_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
    # strand of $reg_loc is ignored!
    if ($reg_loc->genome_db->meta->name eq 'Ensembl_Database_Parameters') {
      my $result = $self->private_check_existence_of_neighbouring_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
      return $result;
    } 
    elsif ($reg_loc->genome_db->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'cannot resolve gene neighbourhood from FASTA file!';
    }
    elsif ($reg_loc->genome_db->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      my $result = $self->private_check_existence_of_neighbouring_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
    }
    else {
      $GU->user_info( 1, $reg_loc->genome_db->meta->name."\n" ); 
      die 'this method is not implemented for this database type yet!';
    }
  } # check_existence_of_neighbouring_gene #

  method get_distance_to_neighbouring_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
    # call check_existence_of_neighbouring_gene first as distance cannot be evaluated if no neighbour exists
    # (this method must throw exception No_Neighbouring_Gene_Exception in this case)
    # strand of $reg_loc is ignored!

      # $GU->user_info(3,"************ Entered get_distance_to_neighbouring_gene: ".$reg_loc->gene_ID.", ".$ignore_pseudogenes.", ".$ignore_RNAgenes.", ".$direction."\n");

    if ($reg_loc->genome_db->isa('Ensembl_Database_Parameters')) {
      my $distance = $self->private_get_distance_to_neighbouring_ensembl_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
      return $distance;
    }
    elsif ($reg_loc->genome_db->isa('FASTA_Sequence_Database_Parameters')) {
      die 'cannot resolve gene distances from FASTA file!';
    }
    elsif ($reg_loc->genome_db->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      my $distance = $self->private_get_distance_to_neighbouring_genbank_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
    }
    else {
      die 'method not implemented yet!';
    }
  } # get_distance_to_neighbouring_gene #

  method get_gene_id_of_neighbouring_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
    # call check_existence_of_neighbouring_gene first as distance cannot be evaluated if no neighbour exists	   
    # strand of $reg_loc is ignored!
    if ($reg_loc->genome_db->isa('Ensembl_Database_Parameters')) {
      my $result = $self->private_get_gene_id_of_neighbouring_ensembl_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
      return $result;
    }
    elsif ($reg_loc->genome_db->isa('FASTA_Sequence_Database_Parameters')) {
      die 'cannot resolve neighbour gene ID from FASTA file!';
    }
    elsif ($reg_loc->genome_db->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      my $result = $self->private_get_gene_id_of_neighbouring_genbank_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
    }
    else {
      die 'method not implemented yet!';
    }
  } # get_gene_id_of_neighbouring_gene #

  method get_id_and_distance_of_neighbouring_gene (Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes) {
  # call check_existence_of_neighbouring_gene first as distance cannot be evaluated if no neighbour exists	
  # strand of $reg_loc is ignored!
    if ($reg_loc->genome_db->isa('Ensembl_Database_Parameters')) {
      my @result = $self->private_get_id_and_distance_of_neighbouring_ensembl_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes);
      return @result;
    } elsif ($reg_loc->genome_db->isa('FASTA_Sequence_Database_Parameters')) {
      die 'cannot resolve neighbour gene ID and distance from FASTA file!';
    } elsif ($reg_loc->genome_db->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      my $result = $self->private_get_id_and_distance_of_neighbouring_genbank_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes);
    } else {
      die 'method not implemented yet!';
    }
  } # get_id_and_distance_of_neighbouring_gene #

  method list_stable_ids_for_a_given_genome_database (Genome_Sequence_Database_Parameters $parameters) {
    if ( $parameters->meta->name eq 'Ensembl_Database_Parameters' ) {
	my @geneIDlist = $self->private_list_stable_ids_for_ensembl_database ($parameters);
	return @geneIDlist;
    } 
    elsif ( $parameters->meta->name eq 'Genbank_Sequence_Database_Parameters' ) {
	my @geneIDlist = $self->private_list_stable_ids_for_genbank_database ($parameters);
	return @geneIDlist;
    }
    elsif ( $parameters->meta->name eq 'FASTA_Sequence_Database_Parameters' ) {
	my @geneIDlist = $self->private_list_stable_ids_for_fasta_database ($parameters);
	return @geneIDlist;
    }
    else {
	die 'Method cannot list stable IDs for this kind of database';
    }
  } # list_stable_ids_for_a_given_genome_database # 

  method get_genomic_sequence(Reg_Loc $reg_loc, Sequence_Parameters $parameters) {
      # to provide flexible retrieval of any genomic sequence, according to parameters
      
      my $five_prime_pos; 
      my $three_prime_pos;
      my $region = $parameters->region;

      # work out boundaries of $reg_loc
      my $reg_loc_start = $reg_loc->position;
      my $reg_loc_length = $reg_loc->get_length_of_regulation_target();
      my $reg_loc_end;      
      if ($reg_loc->strand eq 'positive') {
	  $reg_loc_end = $reg_loc_start+($reg_loc_length-1);
      }
      else {
	  $reg_loc_end = $reg_loc_start-($reg_loc_length-1);
      }
      # work out five prime and three prime positions
      if ($region eq 'generegion') {
	  $five_prime_pos = $reg_loc_start;
	  $three_prime_pos = $reg_loc_end;
      } 
      else {
	  # work out boundaries for upstream and downstream first	 
	  my $upstream_max = $self->get_distance_to_end_of_sequence_region($reg_loc,'towards_five_prime');
	  my $downstream_max = $self->get_distance_to_end_of_sequence_region($reg_loc,'towards_three_prime');
	  my @upstream_distances = ($parameters->max_length_to_search,$upstream_max);
	  my @downstream_distances = ($parameters->max_length_to_search,$downstream_max);
	  if ($parameters->stop_at_neighbouring_gene) {
	      # nota bene: functions for neighbouring gene existence/distance do not interpret genomic direction
	      # relative to strand of $reg_loc (as opposed to get_distance_to_end_of_sequence_region)
	      my $genomic_direction_upstream;
	      my $genomic_direction_downstream;
	      if ($reg_loc->strand eq 'positive') {
		  $genomic_direction_upstream = 'towards_five_prime';
		  $genomic_direction_downstream = 'towards_three_prime';
	      } else {
		  $genomic_direction_upstream = 'towards_three_prime';
		  $genomic_direction_downstream = 'towards_five_prime';
	      }
	      my $upstream_neighbour_gene_exists = $self->check_existence_of_neighbouring_gene($reg_loc,$parameters->ignore_pseudogenes, $parameters->ignore_RNAgenes,$genomic_direction_upstream);
	      if ($upstream_neighbour_gene_exists) {
		  my $distance_to_upstream_neighbour_gene = $self->get_distance_to_neighbouring_gene($reg_loc, $parameters->ignore_pseudogenes, $parameters->ignore_RNAgenes, $genomic_direction_upstream);
		  push(@upstream_distances,$distance_to_upstream_neighbour_gene);		  
	      }
	      my $downstream_neighbour_gene_exists = $self->check_existence_of_neighbouring_gene($reg_loc,$parameters->ignore_pseudogenes, $parameters->ignore_RNAgenes,$genomic_direction_downstream);
	      if ($downstream_neighbour_gene_exists) {
		  my $distance_to_downstream_neighbour_gene = $self->get_distance_to_neighbouring_gene($reg_loc, $parameters->ignore_pseudogenes, $parameters->ignore_RNAgenes, $genomic_direction_downstream);
		  push(@downstream_distances,$distance_to_downstream_neighbour_gene);
	      }	       
	  }
	  my $upstream_distance = $GU->minimum(\@upstream_distances);
	  my $downstream_distance = $GU->minimum(\@downstream_distances);
	  if ($upstream_distance<$parameters->min_length_to_return) {
	      $upstream_distance = $parameters->min_length_to_return;
	  }
	  if ($downstream_distance<$parameters->min_length_to_return) {
	      $downstream_distance = $parameters->min_length_to_return;
	  }
	  my $upstream_boundary_position;
	  my $downstream_boundary_position;
	  if ($reg_loc->strand eq 'positive') {
	      $upstream_boundary_position = $reg_loc_start-($upstream_distance-1);
	      $downstream_boundary_position = $reg_loc_end+($downstream_distance-1);
	  } else {
	      $upstream_boundary_position = $reg_loc_start+($upstream_distance-1);
	      $downstream_boundary_position = $reg_loc_end-($downstream_distance-1);
	  }
	  # sort out cases of other sequence region types
	  if ($region eq 'upstream') {
	      $five_prime_pos = $upstream_boundary_position;
	      $three_prime_pos = $reg_loc_start;
	  }
	  elsif ($region eq 'downstream') {
	      $five_prime_pos = $reg_loc_end;
	      $three_prime_pos = $downstream_boundary_position;
	  } 
	  elsif ($region eq 'upstream_and_generegion') {
	      $five_prime_pos = $upstream_boundary_position;
	      $three_prime_pos = $reg_loc_end;
	  } 
	  elsif ($region eq 'downstream_and_generegion') {
	      $five_prime_pos = $reg_loc_start;
	      $three_prime_pos = $downstream_boundary_position;
	  } 
	  elsif ($region eq 'all') {
	      $five_prime_pos = $upstream_boundary_position;
	      $three_prime_pos = $downstream_boundary_position;
	  } else {
	      die 'method not yet implemented for this sequence region type: '.$region;
	  }
      } 
      # create genomic interval
      $GU->user_info(3, "five prime: ".$five_prime_pos .", three prime: ". $three_prime_pos."\n");
      my $result = Genomic_Interval->new('genome_db' => $reg_loc->genome_db,
					 'region' => $reg_loc->region,
					 'five_prime_pos' => $five_prime_pos,
					 'three_prime_pos' => $three_prime_pos,
					 'strand' => $reg_loc->strand,
					 'label' => $reg_loc->gene_ID);
      if (defined $reg_loc->coord_sys_name) {
	  $result->coord_sys_name($reg_loc->coord_sys_name);
      }
      return $result; # Genomic_Interval 
  } # get_genomic_sequence #

# method get_genomic_unmasked_sequence_of_gi returns unmasked sequence only
  
  method get_genomic_unmasked_sequence_of_gi (Genome_Sequence_Database_Parameters $parameters, Str $coord_sys_name, Str $region, Int $fiveprime, Int $threeprime, StrandType $strand) {
      my $result;
      if ($parameters->isa('Ensembl_Database_Parameters')) {
      	my $reconnect = FALSE;
          my $try_to_connect = 1;
      	while ($try_to_connect) {
      		eval {
      			if ($reconnect eq TRUE){
      				$self->private_reconnect_to_registry($parameters->location,$parameters->dbname);
				$reconnect = FALSE;
      			}
	  			$result = $self->private_get_genomic_sequence_of_ensembl_gi($parameters, $coord_sys_name, $region, $fiveprime, $threeprime, $strand, 'unmasked');
      		};
      		my $error = $@;
            $try_to_connect=0;
      		if ($error) {
                die ($error);
      		}
      	}
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  $result = $self->private_get_genomic_unmasked_sequence_of_genbank_gi($parameters, $coord_sys_name, $region, $fiveprime, $threeprime, $strand);
      }
      elsif ($parameters->isa('FASTA_Sequence_Database_Parameters')) {
	  $result = $self->private_get_genomic_unmasked_sequence_of_fasta_gi($parameters, $coord_sys_name, $region, $fiveprime, $threeprime, $strand);
      }
      else {
	  die 'cannot get genomic sequence for this type of Genome_Sequence_Database_Parameter yet';
      }
      
      return $result;
  } # get_genomic_unmasked_sequence_of_gi #
  
# method get_genomic_repeatmasked_sequence_of_gi returns repeatmasked sequence only  
  
  method get_genomic_repeatmasked_sequence_of_gi (Genome_Sequence_Database_Parameters $parameters, Str $coord_sys_name, Str $region, Int $fiveprime, Int $threeprime, StrandType $strand) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      my $result = $self->private_get_genomic_sequence_of_ensembl_gi($parameters, $coord_sys_name, $region, $fiveprime, $threeprime, $strand, 'repeatmasked');
      return $result;
    }
    elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
#      my $result = $self->private_get_genomic_repeatmasked_sequence_of_genbank_gi ($parameters, $coord_sys_name, $region, $fiveprime, $threeprime, $strand);
#      return $result;
        die "\n   get_genomic_repeatmasked_sequence_of_gi: Repeat masking is not implemented in genbank database!\n   Terminating!";
    }
    else {
      die 'cannot get genomic sequence for this type of Genome_Sequence_Database_Parameter yet';
    } 

  } # get_genomic_repeatmasked_sequence_of_gi #

  method get_orthologous_gene_IDs (Str $query_gene_id, Genome_Sequence_Database_Parameters $source_db_parameters, Genome_Sequence_Database_Parameters $target_db_parameters, Generic_Orthology_Finding_Parameters $ortholog_parameters) { 
    # changes: 1 - needs to take both target and source gsdp as input (e.g. in case we need to produce an RBH mapping): DONE
    # changes: 2 - first check what ortholog parameters are being requested (at the moment we're assuming we want to use Ensembl+Compara).: DONE
    # changes: 3 - call RBH method via job handler in this function : DONE

    # $GU->user_info(3, "getting orthologous gene ids\n");
    # test for orthology method type:
    if ($ortholog_parameters->isa('Compara_Orthology_Finding_Parameters')) {
      # deal with Compara method (first check that target and source are both Ensembl type databases)
      if ($target_db_parameters->isa('Ensembl_Database_Parameters')) { 
	if ($source_db_parameters->isa('Ensembl_Database_Parameters')) {
	  my @result = $self->private_get_orthologous_ensembl_gene_IDs ($query_gene_id, $target_db_parameters, $ortholog_parameters);
	  return @result; # return an array of orthologous ensembl gene IDs
	}
	else {
	  die '$source_db_parameters->dbname is not an Ensembl database, therefore cannot use Compara to find orthologs as requested';
	}
      }
      else {
	die '$target_db_parameters->dbname is not an Ensembl database, therefore cannot use Compara to find orthologs as requested';
      }
    }

    elsif ($ortholog_parameters->isa('RBH_Orthology_Finding_Parameters')) {
      # deal with RBH method (call method here, via Job_Handler)
      my @result = $self->private_get_orthologous_gene_ID_using_RBH_method ($query_gene_id, $source_db_parameters, $target_db_parameters, $ortholog_parameters);
      return @result; # return an orthologous gene ID using RBH method
    }
    elsif ($ortholog_parameters->isa('Random_Assignment_Orthology_Finding_Parameters')) {
      my @result = $self->private_get_orthologous_gene_ID_using_Random_Assignment_Method ($query_gene_id, $source_db_parameters, $target_db_parameters, $ortholog_parameters);
      return @result;
    }    
    elsif ($ortholog_parameters->isa('Synteny_Orthology_Finding_Parameters')) {
      die 'method not implemented yet!';
    }
    else {
      die 'unknown orthology finding parameter type';
    }
  } # get_orthologous_gene_IDs #

  method retrieve_all_sequences_from_DB (Genome_Sequence_Database_Parameters $parameters) {
    # returns all sequences as an array of records, each record has entries NAME and SEQUENCE
    # this method could be further parameterised
    if ( $parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
	my @sequences = $self->private_retrieve_all_sequences_from_fasta_DB($parameters);
	return @sequences;
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
	die 'method not implemented yet'; # 
    }
    elsif ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
	die 'method not implemented yet'; # 
    }
    else {
	die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # retrieve_all_sequences_from_DB #

  method get_all_chromosomes (Genome_Sequence_Database_Parameters $parameters) {
      # returns a genomic interval set, for some genome databases mitochondrial and chloroplast
      # sequences may be excluded

      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  return $self->private_get_all_chromosomes_ensembl($parameters);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  return $self->private_get_all_chromosomes_genbank($parameters);
      }
      else {
	  die 'method not yet implemented for this type of genome database.'
      }
  } # get_all_chromosomes #

  method make_random_gi_set (Genomic_Interval_Set $gi_set, IntPercentage $max_n_ratio, IntPercentage $max_repeat_ratio, Int $seed) {
      # takes a gi_set as input, with n sequences of m lengths, and creates another gi_set of n sequences with m lengths, randomly from the same genome.
	  if (!$gi_set->all_intervals_are_on_same_genome_database()) {
	  die 'did not implement randomisation of mixed GI-sets yet.';
      }
      my @lengths = $gi_set->return_all_gi_lengths();
      my @random_gi_array = ();
      if ($#lengths >= 0) {
	  my $first_gi = ${$gi_set->genomic_interval_set}[0];
	  @random_gi_array = $self->get_set_of_random_genomic_intervals($first_gi->genome_db,\@lengths,$max_n_ratio,$max_repeat_ratio,$seed);
      }
      my $result = Genomic_Interval_Set->new(genomic_interval_set => \@random_gi_array);
      return $result;

  } # make_random_gi_set #

  method random_genomic_interval_picker (Genome_Sequence_Database_Parameters $parameters, Int $number, IntPercentage $max_n_ratio, IntPercentage $max_repeat_ratio, Int $length, Int $seed ) {
      # number of sequences to pick
      # max % Ns in assembly to allow
      # max % Ns in repeatmasked sequence to allow
      # length of sequence to pick
      # the random seed
      
      my @lengths = $GU->repeat_array_element($length,$number);
      return $self->get_set_of_random_genomic_intervals($parameters,\@lengths,$max_n_ratio,$max_repeat_ratio,$seed);
  } # random_genomic_interval_picker #
  
  method all_transcripts_to_fasta_file (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self->private_ensembl_all_transcripts_to_fasta_file ($parameters);
    }
    elsif ($parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    else {
      die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # all_transcripts_to_fasta_file #

  method one_transcript_per_gene_to_fasta_file (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self->private_ensembl_all_transcripts_to_fasta_file ($parameters);
    }
    elsif ($parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      $self->private_genbank_all_transcripts_to_fasta_file ($parameters);
    }
    else {
      die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # one_transcript_per_gene_to_fasta_file #

 method one_CDS_per_gene_to_fasta_file (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self->private_ensembl_one_CDS_per_gene_to_fasta_file ($parameters);
    }
    elsif ($parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    else {
      die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # one_CDS_per_gene_to_fasta_file #

  method genome_to_fasta_file (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self -> private_ensembl_genome_to_fasta_file ($parameters);
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      $self -> private_genbank_genome_to_fasta_file($parameters) 
    }
    elsif ($parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    else {
      die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # genome_to_fasta_file #

  method masked_genome_to_fasta_file (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self->private_ensembl_masked_genome_to_fasta_file ($parameters);
    }
    elsif ($parameters->meta->name eq 'FASTA_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    elsif ($parameters->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      die 'method not implemented yet'; # 
    }
    else {
      die 'cannot retrieve sequences from this kind of genome sequence database'; 
    }
  } # masked_genome_to_fasta_file #

  method fasta_2_raw (Str $filename) {
	
    my $in = Bio::SeqIO->new (-file => $filename , '-format' => 'Fasta') || die 'not FASTA';
    my $out = Bio::SeqIO->new (-file => '>'. $filename.'.raw' , '-format' => 'raw');
    
    while ( my $seq = $in->next_seq() ) {
      $out->write_seq($seq);
    }
    
  } # fasta_2_raw #
  
  method database_stats (Genome_Sequence_Database_Parameters $parameters) {
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      $self->private_ensembl_database_stats ($parameters);
    }
    else {
      die 'cannot produce stats for this kind of Genome_Sequence_Database_Parameters object yet';
    }
  } # database_stats #

  method get_natural_name (Str $geneid, Genome_Sequence_Database_Parameters $parameters) {	  
    if ($parameters->isa('Ensembl_Database_Parameters'))  {
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);

      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
      my $gene = $gene_adaptor->fetch_by_stable_id($geneid);
      my $natural_name = $gene->external_name();
      if (!$natural_name || $natural_name eq '') {
	$natural_name = 'no natural name';
      }
      return $natural_name;
    }
    else {
      die 'cannot get natural name for gene ids other than ensembl format';
    }
  } # get_natural_name #
  
  method pick_random_gene_ids (Int $seed, Int $number_to_pick, Genome_Sequence_Database_Parameters $parameters) {
    
    if ($parameters->meta->name eq 'Ensembl_Database_Parameters') {
      #my $alias = $parameters->alias; # redundant
      
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      
      # get array of all gene ids in database
      my @all_gene_ids = $self->list_stable_ids_for_a_given_genome_database($parameters);
      my @random_ids;
      my $available_gene_ids = scalar(@all_gene_ids);
      #$GU->user_info(3,"$available_gene_ids available\n");
      
      if ($number_to_pick > $available_gene_ids) {
	die "you asked for $number_to_pick gene ids, but there are only $available_gene_ids gene ids in the specified database";
      }
      
      # using user-specified random seed, generate a random number based on the available_gene_ids number
      srand($seed);
      my $i;
      for ($i = 0; $i < $number_to_pick; $i++) {
	#$GU->user_info(3,"count = $i\n");
	$available_gene_ids = scalar(@all_gene_ids); # how many genes are left
	#$GU->user_info(3,"$available_gene_ids ids are left\n");
	my $index = int(rand($available_gene_ids)); # select a random number to provide the index
	#$GU->user_info(3,"random index = $index\n");
	my $selected_id = $all_gene_ids[$index]; # select id of relevant indexprint "SELECTED: $selected_id\n";
	#$GU->user_info(3,"selected id = $selected_id\n");
	push (@random_ids, $selected_id); # add to random_id array
	splice (@all_gene_ids, $index, 1);# remove (splice) single id from all gene id array
      }
      return @random_ids; # return type: array
    }
    else {
      die 'cannot pick random gene ids for non-Ensembl databases yet';
    }
  } # pick_random_gene_ids #

  method proteins_to_fasta_file (Genome_Sequence_Database_Parameters $parameters, Str $filename, TranscriptChoice $transcript_choice) {
    # calls individual private methods for different GSDP types. 
    # returns 1; when FASTA file has been written
    # this method is 'direct'; see protein_intervals_to_fasta_file for a similar method that first creates GI_Set.
    if ($parameters->isa('Ensembl_Database_Parameters')) {
      $self->private_ensembl_proteins_to_fasta_file ($parameters, $filename, $transcript_choice);
    }
    elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
      $self->private_genbank_proteins_to_fasta_file ($parameters, $filename,$transcript_choice);
    }
    elsif ($parameters->isa('FASTA_Sequence_Database_Parameters')) {
      copy($parameters->filename, $filename);
    } else {
    	die 'Cannot make a FASTA file for proteins from this kind of database yet - method not implemented';
    }
    return 1;
  } # proteins_to_fasta_file #

  method check_dna_fasta_file (Str $filename) {
    # checks that the file is all IUPAC
    my $result = TRUE;
    my $seqio = Bio::SeqIO->new(-file => $filename);
    while (my $seq = $seqio->next_seq() ) {
      if ( $seq->seq =~ /[^ATGCRYSWKMDHBVN]/i ) {			
	$result = FALSE;
      }
    }
    return $result;
  } # check_dna_fasta_file #
  
  method count_number_of_characters_fasta_file (Str $filename) {
    my $seqio = Bio::SeqIO->new(-file => $filename);
    my $result = 0;
    while (my $seq = $seqio->next_seq() ) {
      $result = $result+length ($seq->seq);
    }
    return $result;
  } # count_number_of_characters_fasta_file # 
  
  method roughly_check_protein_fasta_file (Str $filename) {
    # tests rough criterion to decide if sequence is protein or not
    my $result = TRUE;
    my $seqio = Bio::SeqIO->new(-file => $filename);
    while (my $seq = $seqio->next_seq() ) {
	my $ATGCNcount = ( $seq->seq =~ tr/nNaAtTgGcC//);
	my $length = length ($seq->seq);
	my $percentage = ($ATGCNcount/$length)*100;
	if (($percentage > 90)&&($length>100)) {		
	  $result = FALSE;
	}
      }
    return $result;
  } # roughly_check_protein_fasta_file #

  method get_coord_system_name_from_gene_id (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  	return $self->private_get_coord_system_name_from_ensembl_gene_id($parameters,$geneid);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  	  my $geneinf = $self->private_get_genbank_geneinfo($geneid);
    	  return $geneinf->{COORD_SYSTEM};
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_coord_system_name_from_gene_id #

  method get_seq_region_name_from_gene_id (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
		  return $self->private_get_seq_region_name_from_ensembl_gene_id($parameters,$geneid);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  	  my $geneinf = $self->private_get_genbank_geneinfo($geneid);
    	  return $geneinf->{REGION};
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_seq_region_name_from_gene_id #

  method get_gene_start_position (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  return $self->private_get_ensembl_gene_start_position($parameters,$geneid);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
    	  return $self->private_get_genbank_gene_start_position($geneid);
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_gene_start_position #

  method get_gene_end_position (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  return $self->private_get_ensembl_gene_end_position($parameters,$geneid);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  die 'method not yet implemented for GenBank databases';
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_gene_end_position #

  method get_gene_length (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      my $gene_start =  $self->get_gene_start_position($parameters,$geneid);
      my $gene_end = $self->get_gene_end_position($parameters,$geneid);
      my $result = $gene_end-$gene_start;
      if ($result<0) {
	  $result = -$result;
      }
      $result++;
      return $result;
  } # get_gene_length #

  method get_gene_parameters (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  die 'method not yet implemented for this type of genome database';
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
    	  return $self->private_get_genbank_gene_parameters($parameters,$geneid);
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_gene_parameters #

  method get_gene_go_terms (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
    # Gets the go terms for each of the genes in the geneome database
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  die 'method not yet implemented for this type of genome database';
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
    	  return $self->private_get_genbank_go_terms($parameters,$geneid);
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_gene_go_terms #

  method get_gene_proteins (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
    # Gets the go terms for each of the genes in the geneome database
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  die 'method not yet implemented for this type of genome database';
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
    	  return $self->private_get_genbank_proteins($parameters,$geneid);
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_gene_go_terms #

  method get_strand_of_gene (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  return $self->private_get_strand_of_ensembl_gene($parameters,$geneid);
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  my $geneinf = $self->private_get_genbank_geneinfo($geneid);
    	  return $geneinf->{STRAND};
      }
      else {
	  die 'method not yet implemented for this type of genome database';
      }
  } # get_strand_of_gene #

  method protein_intervals_to_fasta_file (Genome_Sequence_Database_Parameters $parameters, Str $filename, TranscriptChoice $transcript_choice) {
    # calls individual private methods for different GSDP types. 
    # returns 1; when FASTA file has been written 

    # converts ensembl proteins to a FASTA file, of the name provided as input
    # need to re-implement as an alternative (GI_set-based) method, as follows:
    # makes an array of Genomic_Intervals
    # makes a GI_Set
    # calls $GI_Set->render_as_FASTA (or create_FASTA_from_genomic_interval_set)
    # returns filename to FASTA file (full path)

    die 'method not implemented yet!';

    if ($parameters->isa('Ensembl_Database_Parameters')) {
      $self->private_ensembl_protein_intervals_to_fasta_file ($parameters, $filename, $transcript_choice);
    }
    elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
      $self->private_genbank_protein_intervals_to_fasta_file ($parameters, $filename);
    }
    else {
      die 'Cannot make a FASTA file for proteins from this kind of database yet - method not implemented';
    }
    return 1;

  } # protein_intervals_to_fasta_file #

  method get_display_name_for_database(Genome_Sequence_Database_Parameters $parameters) {
    # "nice" label for use in plots

    my $result;
    if ($parameters->isa('Ensembl_Database_Parameters')) {
      if (defined $parameters->alias) {
	$result = $parameters->alias();
      }
      else {
	$result = $parameters->dbname;
      }
    }
    elsif ($parameters->isa('FASTA_Sequence_Database_Parameters')) {
      if (defined $parameters->natural_species_name) {
	$result = $parameters->natural_species_name;
      } else {
	$result = $parameters->dbname;
      }
    }
    elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
      $result = $parameters->dbname;
    }
    else {
      die 'not implemented display name for this type of genome database yet!';
    }
    return $result;
  } # get_display_name_for_database #

  method get_ID_for_database(Genome_Sequence_Database_Parameters $parameters) {
    my $result;
    if ($parameters->isa('Ensembl_Database_Parameters')) {
      $result = $parameters->dbname();
    }
    elsif ($parameters->isa('FASTA_Sequence_Database_Parameters')) {
      $result = $parameters->filename();
    }
    elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
      $result = $parameters->dbname();
    }
    else {
      die 'not implemented display name for this type of genome database yet!';
    }
    return $result;
  } # get_ID_for_database #

  method tell_species_name_from_database ($database_array_ref) {
    #die 'method not implemented yet';
    my @APPLESSpeciesNames;
    foreach my $database (@{$database_array_ref}) {
      if ($database->isa('Ensembl_Database_Parameters')) {
	if (defined $database->alias) {
	  push (@APPLESSpeciesNames, $database->alias);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no alias has been provided for this Ensembl database ($database->dbname)';# is dying too harsh?
	}
      }
      elsif ($database->isa('FASTA_Sequence_Database_Parameters')) {
	if (defined $database->natural_species_name) {
	  push (@APPLESSpeciesNames, $database->natural_species_name);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no natural_species_name has been provided for this FASTA database ($database->dbname)';# is dying too harsh?
	}
      }
      elsif ($database->isa('Genbank_Sequence_Database_Parameters')) {
	if (defined $database->natural_species_name) {
	  push (@APPLESSpeciesNames, $database->natural_species_name);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no natural_species_name has been provided for this Genbank database ($database->dbname)';# is dying too harsh?
	}
      }
      else {
	die 'not implemented tell_species_name_from_database method for this type of genome database yet!';      
      }
    }
    return @APPLESSpeciesNames; # returns array of species names (@[APPLESSpeciesNames])
  } # tell_species_name_from_database #

  method tell_species_names_from_database_hash (HashRef $database_hash_ref) {
    my @APPLESSpeciesNames;
    my @set_of_keys = keys %{$database_hash_ref};
    my %db_hash = %{$database_hash_ref};
    foreach my $database (@set_of_keys) {
      if ($db_hash{$database}->isa('Ensembl_Database_Parameters')) {
	if (defined $db_hash{$database}->alias) {
	  push (@APPLESSpeciesNames, $db_hash{$database}->alias);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no alias has been provided for this Ensembl database ($database->dbname)';# is dying too harsh?
	}
      }
      elsif ($db_hash{$database}->isa('FASTA_Sequence_Database_Parameters')) {
	if (defined $db_hash{$database}->natural_species_name) {
	  push (@APPLESSpeciesNames, $db_hash{$database}->natural_species_name);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no natural_species_name has been provided for this FASTA database ($database->dbname)';# is dying too harsh?
	}
      }
      elsif ($db_hash{$database}->isa('Genbank_Sequence_Database_Parameters')) {
	if (defined $db_hash{$database}->natural_species_name) {
	  push (@APPLESSpeciesNames, $db_hash{$database}->natural_species_name);
	}
	else {
	  die 'cannot determine APPLESSpeciesName - no natural_species_name has been provided for this Genbank database ($database->dbname)';# is dying too harsh?
	}
      }
      else {
	die 'not implemented tell_species_name_from_database method for this type of genome database yet!';      
      }
    }
    return @APPLESSpeciesNames; # returns array of species names (@[APPLESSpeciesNames])
  } # tell_species_names_from_database_hash #

###################################
##### generic private methods #####
###################################

  method private_get_neighbouring_gene_information(Genome_Sequence_Database_Parameters $genome_db,
						   Str $seq_region,
						   Int $position,
						   Boolean $ignore_pseudogenes,
						   Boolean $ignore_RNAgenes, # [HM] introduced in May 2010
						   Str $gene_id_to_ignore) {
    #$GU->user_info(3, "calling  private_get_neighbouring_gene_information\n"); # debugging
    #$GU->user_info( 3, "gene id to ignore: ".$gene_id_to_ignore."\n");
    # computes distances to nearest neighbours to the left and to the right (on positive strand)
    # if $gene_id_to_ignore is a valid gene-ID, then transcripts of this gene are ignored
    # returns an array of length 4 in the format (left neighbour, distance, right neighbour, #distance)
    # if no neighbour gene exists, 'NONE' is returned instead of the distance
    
    my @all_gene_positions;
    
    if ($genome_db->isa('Ensembl_Database_Parameters')) {
      $GU->user_info(3, "getting all ensembl gene positions via job handler and cache\n");
      my $cache = TRUE;
      my $job_handler = Job_Handler->new();
      $job_handler->get_config_settings();

      # a fix for the problem that gene positions are not taken from the cache when the
      # orthology methods change (which do not affect gene positions as of 9/4/10)
      # not a particularly nice fix, though
      my @empty_array = ();
      my $copy_of_genome_db = Ensembl_Database_Parameters->new(
							       orthology_methods => \@empty_array,
							       dbname => $genome_db->dbname,
							       alias => $genome_db->alias,
							       location => $genome_db->location);
      
      push (my @parameters, $copy_of_genome_db);
      my $function = 'private_get_all_ensembl_gene_positions';
      my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
					       wall_time_estimate => 172800,
					       cache_in_memory => TRUE
					      );
      my $GDBU = Genome_DB_Utilities->new();
      my @cache_result = eval { 
	$job_handler->handle_APPLES_function($function, $GDBU, \@parameters, $cache, 365, $job_parameters); 
      };
      
      my $exception_info = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler
      if ($@) {
	my $exception_content = $@;
	if (!UNIVERSAL::can($exception_content, 'isa')) { # need to first test error is an object with method 'isa' (else isa fails with another exception!)
	  die $exception_content; # throw the error string
	} 
	else {
	  if ($exception_content->isa('Job_Information_Exception')) {
	    $exception_info->merge($exception_content);
	    $GU->user_info(1,"this function (private_get_neighbouring_gene_information) can not be completed until the result is available\n");
	    $GU->user_info(1,"throwing an exception!\n");
	    die $exception_info;
	  }
	  else {
	    die $exception_content; # throw the error object
	  }
	}
      }
      # then to handle the returned result (in this case we are not interested in the calling object):
      my $object = shift(@cache_result); # first element of array is the object - shift it off the array
      @all_gene_positions = @cache_result; # the rest of the array is now the result
      #$GU->user_info(3, "all gene positions\n"); # for debugging
      #$GU->user_info(3, Dumper (\@all_gene_positions)); # for debugging
      
    } # endif Ensembl_Database_Parameters
    elsif ($genome_db->isa('FASTA_Sequence_Database_Parameters')) {
      die 'cannot resolve neighbour gene ID from FASTA file!';
    }
    elsif ($genome_db->meta->name eq 'Genbank_Sequence_Database_Parameters') {
      #@all_gene_positions = $self->private_get_all_genbank_gene_positions($genome_db->version); # NOW WRAPPED VIA CALL INTO JOB_HANDLER
      $GU->user_info(3, "getting all genbank gene positions via job handler and cache\n");
      my $cache = TRUE;
      my $job_handler = Job_Handler->new();
      $job_handler->get_config_settings();
      push (my @parameters, $seq_region);
      my $function = 'private_get_all_genbank_gene_positions';
      my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
					       wall_time_estimate => 172800
					      );
      my $GDBU = Genome_DB_Utilities->new();
      my @cache_result = eval { 
	$job_handler->handle_APPLES_function($function, $GDBU, \@parameters, $cache, 365, $job_parameters); 
      };
      
      my $exception_info = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler
      if ($@) {
	my $exception_content = $@;
	if (!UNIVERSAL::can($exception_content, 'isa')) { # check $exception_content has method 'isa'
	  die $exception_content; # throw error string
	} else {
	  if ($exception_content->isa('Job_Information_Exception')) {
	    $exception_info->merge($exception_content);
	    $GU->user_info(1,"this function (private_get_neighbouring_gene_information) can not be completed until result is available\n");
	    $GU->user_info(1,"throwing an exception!\n");
	    die $exception_info;
	  }
	  else {
	    die $exception_content; # throw error object
	  }
	} 
      }
      # then to handle the returned result (in this case we are not interested in the calling object):
      my $object = shift(@cache_result); # first element of array is the object - shift it off the array
      @all_gene_positions = @cache_result; # the rest of the array is now the result
      #$GU->user_info(3, "all gene positions\n"); # for debugging
      #$GU->user_info(3, Dumper (\@all_gene_positions)); # for debugging
    }
    my @typesofpseudogenes = (
			      'snoRNA-pseudogene',  'scRNA-pseudogene',
			      'snRNA_pseudogene',   'scRNA_pseudogene',
			      'snoRNA_pseudogene',  'rRNA_pseudogene',
			      'Mt_tRNA_pseudogene', 'tRNA_pseudogene',
			      'miRNA_pseudogene',   'misc_RNA_pseudogene',
			      'polymorphic_pseudogene', 'unprocessed_pseudogene',
			      'processed_pseudogene'
			     );
    my @typesofRNAgenes = ('snRNA','miRNA','misc_RNA','snoRNA','rRNA');
    my @result = ( 'NONE', 'NONE', 'NONE', 'NONE' );
    foreach my $entry (@all_gene_positions) {
      if ( $entry->{CHROMOSOME} eq $seq_region ) {
	if (   ( $entry->{TYPE} ne 'pseudogene' )
	       || ( !$ignore_pseudogenes ) ) {
	  # check for pseudo-genes
	  my $isakindofpseudogene = FALSE;
	  map {
	    if ( $_ eq $entry->{TYPE} ) {
	      $isakindofpseudogene = TRUE;
	    }
	  } @typesofpseudogenes;
	  
	  # check for RNA-genes
	  my $isakindofRNAgene = FALSE;
	  map {
	    if ( $_ eq $entry->{TYPE} ) {
	      $isakindofRNAgene = TRUE;
	    }
	  } @typesofRNAgenes;
	  
	  if (( !$isakindofpseudogene ) 
	      || ( !$ignore_pseudogenes ) ) {
	    if ((!$isakindofRNAgene)||
		(!$ignore_RNAgenes)) {
	      if ( $entry->{ID} ne $gene_id_to_ignore ) {
		if ( $entry->{START} < $position ) {
		  my $leftdistance = $position - $entry->{END};
		  if (
		      ( $result[1] eq 'NONE' )
		      || (   ( $result[1] ne 'NONE' )
			     && ( $result[1] > $leftdistance ) )
		     ) {
		    $result[0] = $entry->{ID};
		    $result[1] = $leftdistance;
		  }
		}
		if ( $entry->{END} > $position ) {
		  my $rightdistance = $entry->{START} - $position;
		  if (
		      ( $result[3] eq 'NONE' )
		      || (   ( $result[3] ne 'NONE' )
			     && ( $result[3] > $rightdistance ) )
		     ) {
		    $result[2] = $entry->{ID};
		    $result[3] = $rightdistance;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    if ($result[1] ne 'NONE') {
      if ( $result[1] < 0 ) {
	$result[1] = 0;
      }
    }
    if ($result[3] ne 'NONE') {
      if ( $result[3] < 0 ) {
	$result[3] = 0;
      }
    }
    return @result;   
  } # private_get_neighbouring_gene_information #
  
  method private_check_existence_of_neighbouring_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
      my $result;
      eval {
	  $self->get_distance_to_neighbouring_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction);
      };
      if ($@) {
	  my $exception_content = $@;
	  if (!UNIVERSAL::can($exception_content, 'isa')) {
	      die $exception_content; # if exception can't 'isa', re-throw the exception
	  }
	  else {
	      if ($exception_content->isa('No_Neighbouring_Gene_Exception')) {
		  $GU->user_info(3,"encountered gene without a neighbour.\n");
		  $result = FALSE; 
	      }
	      else {
		  die $exception_content;
	      }
	  }
      }
      else {
	  $result = TRUE;
      }
      return $result;
  } # private_check_existence_of_neighbouring_gene # 

  method get_set_of_random_genomic_intervals (Genome_Sequence_Database_Parameters $parameters, ArrayRef $lengths, IntPercentage $max_n_ratio, IntPercentage $max_repeat_ratio, Int $seed) {
      if ($parameters->isa('Ensembl_Database_Parameters')) {
	  return $self->get_set_of_random_genomic_intervals_ensembl($parameters,$lengths,$max_n_ratio,$max_repeat_ratio,$seed);
      }
      elsif ($parameters->isa('FASTA_Sequence_Database_Parameters')) {
	  die 'method not implemented yet'; # 
      }
      elsif ($parameters->isa('Genbank_Sequence_Database_Parameters')) {
	  die 'method not implemented yet'; # 
      }
      else {
	  die 'method not implemented yet'; 
      }
  } # get_set_of_random_genomic_intervals #

############################################
##### Ensembl-specific private methods #####
############################################

  method get_set_of_random_genomic_intervals_ensembl (Ensembl_Database_Parameters $parameters, ArrayRef $lengths, IntPercentage $max_n_ratio, IntPercentage $max_repeat_ratio, Int $seed) {
      
      $GU->user_info(2,"Starting to pick random intervals.\n");
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $slice_adaptor = $self->private_get_ensembl_slice_adaptor($parameters);

      # get genome length
      my $genome_length = 0;
      my @slices = (@{ $slice_adaptor->fetch_all('toplevel') } );
      foreach my $slice (@slices) {
	$GU->user_info( 3, "slice length:\t".$slice->length."\n" );
	$genome_length = $genome_length + $slice->length;
	$GU->user_info( 3, "sum so far = ".$genome_length."\n" );
      }
      $GU->user_info( 3, "total genome length: ".$genome_length."\n" );

      # pick random intervals
      my $genomic_interval;
      my @random_genomic_interval_set;
      srand($seed);
      my $number = @{$lengths};
      my $progress_count = 0;
      for (my $i = 0; $i < $number; $i++) {
	  my $length = ${$lengths}[$i];
	  my $random_position = int (rand ($genome_length));
	  my $selected_slice;
	  my $selected_sequence_unmasked;
	  my $selected_sequence_masked;
	  $GU->user_info( 3, "random position = ".$random_position."\n" );
	  foreach my $slice (@slices) {
	      if ($random_position <= $slice->length) {
		  $GU->user_info( 3, "position is on this slice:".$slice->name."\n" );
		  $selected_slice = $slice;
		  last;
	      }
	      else {
		  $random_position = $random_position - $slice->length;
		  $GU->user_info( 3, "random position = ".$random_position."\n" );
	      }
	  }
	  if ($random_position + $length < $selected_slice->length) {
	      $GU->user_info( 3, $random_position."\t".$length."\t".$selected_slice->length."\tcan get all sequence\n" );
	      my $coord_sys_name = $selected_slice->coord_system()->name();
	      my $region = $selected_slice->seq_region_name;
	      $selected_sequence_unmasked = $selected_slice->subseq($random_position, ($random_position+$length-1));
	      $GU->user_info( 3, "Unmasked: ".$selected_sequence_unmasked."\n" );	      	     
	      my $chosen_part_of_slice = $slice_adaptor->fetch_by_region($coord_sys_name,$region,$random_position,($random_position+$length-1),1); # last parameter is for strand
	      my $hardmasked_slice = $chosen_part_of_slice->get_repeatmasked_seq(); # changes type to Bio::EnsEMBL::RepeatMaskedSlice
	      $selected_sequence_masked = $hardmasked_slice->seq();
	      $GU->user_info( 3, "Masked: ".$selected_sequence_masked."\n" );
	      my $ncountunmasked = ( $selected_sequence_unmasked =~ tr/nN// );
	      my $ncountmasked = ( $selected_sequence_masked =~ tr/nN// );
	      $GU->user_info( 3, "masked N's=\t".$ncountmasked."\n" );
	      $GU->user_info( 3, "percent masked=\t".(($ncountmasked/$length) *100)."\n" );
	      $GU->user_info( 3, "max percent masked=\t".$max_repeat_ratio."\n" );
	      $GU->user_info( 3, "number Ns in unmasked sequence: ".$ncountunmasked."\n" );
	      if (( ($ncountunmasked/$length) * 100 > $max_n_ratio) || ( ($ncountmasked/$length) * 100 > $max_repeat_ratio )) {
		  # reject choice - too many Ns in unmasked or masked sequence
		  $i--;
	      }
	      else {
		  $genomic_interval = Genomic_Interval->new(
		      'genome_db' => $parameters,
		      'coord_sys_name' => $coord_sys_name,
		      'region' => $region,
		      'five_prime_pos' => $random_position,
		      'three_prime_pos' => ($random_position+$length-1),
		      'strand' => 'positive',
		      'gi_sequence' => $selected_sequence_unmasked,
		      'gi_sequence_repeatmasked' => $selected_sequence_masked
		      ); 
		  push (@random_genomic_interval_set, $genomic_interval); 
		  $progress_count++;
		  $GU->user_info(2,"Done ".$progress_count." out of ".$number.".\n");
	      }
	  }
	  else {
	      # reject choice - not enough sequence to pick
	      $i--;
	  }
      }
      return @random_genomic_interval_set;
  } # get_set_of_random_genomic_intervals_ensembl #
      
  method private_connect_to_registry (Str $location, Str $dbname) {
      my $registry = 'Bio::EnsEMBL::Registry';
      

      if (defined $self->registry) {
	  if ($self->current_registry_location ne $location) {
	      $GU->user_info(3,"Reconnecting to registry because location changed from: \'".$self->current_registry_location."\' to \'".$location."\' (database name is: ".$dbname.").\n");
	      $self->private_reconnect_to_registry ($location, $dbname)
	  }
      }

      if (!defined $self->registry) {
	  if ($location eq 'ensembl') { 
	      $registry->load_registry_from_db(
		  -host => 'ensembldb.ensembl.org',
		  -user => 'anonymous'
		  #-verbose => 1 # optional (gives verbose output, useful for checking what DBs are loaded)
		  ); # ONLY loads the databases from ensembl that are available with your installed API software version
	  } 
	  elsif ($location eq 'ensemblgenomes') {
	      $registry->load_registry_from_db(
		  -host => "mysql.ebi.ac.uk",
		  -port => 4157,
		  -user => 'anonymous',
		  -wait_timeout => 86400,
		  ); # testing 24-hour timeout (default is 8 hours)
	       # ONLY loads the databases from ensemblgenomes that are available with your installed API software version
	      $registry->set_disconnect_when_inactive();
	  } 
	  elsif ($location eq 'local') {
      print "\nGenome_DB_Utilities line 1216.\n";
	    $GU->user_info(3, "Requested DB is local\n");#debugging
	    my $ensembl_registry_conf = "/home/grannysmith/webseaweeds/Configuration/ensembl_registry_voland.conf"; #[fixed] - nd
      print "\nGenome_DB_Utilities line 1219.\n";
      $registry->load_all($ensembl_registry_conf, 1); # 1 (optional) gives verbose output

	    # die if api/database versions mismatch
	    $self->private_ensembl_api_and_database_consistent_versions_check($dbname);
	  }	      
	  else {
	    die 'unknown ensembl location\n';
	  }
	  $self->registry($registry); # set the registry attribute
	  $self->current_registry_location($location); # set the current_registry_location attribute
      }
  } # private_connect_to_registry #
  
  method private_reconnect_to_registry (Str $location, Str $dbname) {
    $GU->user_info(3,"Re-connecting to registry.\n");
    if (defined $self->registry){
      $self->registry->clear(); # will clear the registry and disconnect from all databases (but does not seem to be enough to drop registry connection).
    }
    $self->private_clear_registry();# clear registry attribute
    $self->private_clear_location();# clear the current_registry_location attribute
    $self->private_connect_to_registry ($location, $dbname); # reconnect afresh  
  } # private_reconnect_to_registry #
  
# Sascha: "But why is the consistency check only done in one out of three cases within private_connect_to_registry?"
#
# Laura: "The other 2 methods (which use the method load_registry_from_db) have a built in mechanism - they only
# load databases into the $registry object for the software API that you have installed (this is automatically detected).
# It's therefore now up to the user to select the appropriate database version to connect to, considering their installed
# API version, so that the releases match - in all cases I tested where there was a mismatch between the requested
# database version and my installed API, it generated a crash, so I thought the API version checker was redundant in these
# cases."

  method private_get_distance_to_end_of_sequence_region_ensembl(Reg_Loc $reg_loc, GenomicDirection $direction) {
      # $direction interpreted relative to $reg_loc
      
      my $result;
      my $positive_and_5_prime = ($reg_loc->strand eq 'positive')&&($direction eq 'towards_five_prime');
      my $negative_and_3_prime = ($reg_loc->strand eq 'negative')&&($direction eq 'towards_three_prime');
      if ($positive_and_5_prime || $negative_and_3_prime) {
	  $result = $reg_loc->position;
      } else {
	  $self->private_connect_to_registry($reg_loc->genome_db->location, $reg_loc->genome_db->dbname);
	  my $slice_adaptor = $self->private_get_ensembl_slice_adaptor($reg_loc->genome_db);
	  my $slice;
	  if (defined $reg_loc->coord_sys_name) {
	      $slice = $slice_adaptor->fetch_by_region($reg_loc->coord_sys_name, $reg_loc->region);
	  } else {
	      die 'method not implemented for case of missing coordinate system name.';
	  }
	  my $region_length = $slice->length;
	  $result = $region_length - $reg_loc->position;
      }
      return $result;
  } # private_get_distance_to_end_of_sequence_region_ensembl #

 method private_get_distance_to_neighbouring_ensembl_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {

    my $result = $self->private_get_distance_or_id_neighbouring_ensembl_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction, TRUE); 
    return $result;
  } # private_get_distance_to_neighbouring_ensembl_gene #

 method private_get_gene_id_of_neighbouring_ensembl_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {

    my $result = $self->private_get_distance_or_id_neighbouring_ensembl_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction, FALSE);
    return $result;
  } # private_get_gene_id_of_neighbouring_ensembl_gene # 

 method private_get_distance_or_id_neighbouring_ensembl_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction, Boolean $get_distance) {
	
    my $parameters = $reg_loc->genome_db;
    my $position = $reg_loc->position;
    my $gene_id_to_ignore;
    if ($reg_loc->gene_ID) {  
      $gene_id_to_ignore = $reg_loc->gene_ID;
    } 
    else {
      $gene_id_to_ignore = '';
    }

    my @neighbouring_gene_information = $self->private_get_neighbouring_gene_information($reg_loc->genome_db,
											 $reg_loc->region,
											 $position,
											 $ignore_pseudogenes,
											 $ignore_RNAgenes, # [HM] introduced in May 2010
											 $gene_id_to_ignore);
    my $relevant_index = 1;
    if ($direction eq 'towards_three_prime') {
      $relevant_index = 3;
    }
    my $relevant_item = $neighbouring_gene_information[$relevant_index];
    if ($relevant_item eq 'NONE') {
	my $exception = No_Neighbouring_Gene_Exception->new();
	die $exception;
    }
    if ($get_distance) {
	return $relevant_item;
    } 
    else {
	my $id = $neighbouring_gene_information[$relevant_index-1];
	return $id;
    }
  } # private_get_distance_or_id_neighbouring_ensembl_gene # 

method private_get_id_and_distance_of_neighbouring_ensembl_gene (Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes) {
    my $gene_id_to_ignore = $reg_loc->gene_ID;
    my @neighbouring_gene_info = $self->private_get_neighbouring_gene_information($reg_loc->genome_db,
										  $reg_loc->region,
										  $reg_loc->position,
										  $ignore_pseudogenes,
										  $ignore_RNAgenes, # [HM] introduced in May 2010
										  $gene_id_to_ignore);
    return @neighbouring_gene_info;
  } # private_get_id_and_distance_of_neighbouring_ensembl_gene #

 method private_get_all_ensembl_gene_positions (Ensembl_Database_Parameters $parameters) {
      # this method should always be called via the Job_Handler for efficiency, using the handle_APPLES_function method
      
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      
      my @all_gene_positions;
      
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
      my @genelist     = @{ $gene_adaptor->list_stable_ids() };
      $GU->user_info(3,"obtained gene list.\n");
      
      my $progress_count = 0;
      foreach my $gene (@genelist) {
	  $progress_count++;
	  if (($progress_count % 50) == 0) {
	      $GU->user_info(3,$progress_count." genes processed.\n");
	  }
	  my $geneobject;
	  my $seq_region;
	  my $start;
	  my $end;
	  my $typeofgene;
	  my $reconnect = FALSE;
	  while (!$main::global_perseverance->stop_trying) {
	      eval {
		  if ($reconnect) {
		      $self->private_reconnect_to_registry($parameters->location, $parameters->dbname);
		      $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
		      $reconnect = FALSE;
		  }
		  $geneobject = $gene_adaptor->fetch_by_stable_id($gene);
		  $seq_region = $geneobject->slice->seq_region_name();
		  $start      = $geneobject->start();
		  $end        = $geneobject->end();
		  $typeofgene = $geneobject->biotype();
	      };
	      my $error = $@;
	      $main::global_perseverance->decide_on_rerun('Ensembl', TRUE, $error);
	      if ($error) {
		  $reconnect = TRUE;
	      }
	  }
	  $main::global_perseverance->stop_trying(FALSE);
	  my $entry = {
	      ID         => $gene,
	      CHROMOSOME => $seq_region, 
	      # might also be the ID of a scaffold (or other type of sequence region)
	      START      => $start,
	      END        => $end,
	      TYPE       => $typeofgene,
	  };
	  push( @all_gene_positions, $entry );
      }
      return @all_gene_positions; # return type: array
  } # private_get_all_ensembl_gene_positions #

  method private_list_stable_ids_for_ensembl_database (Ensembl_Database_Parameters $parameters) {
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);
    my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
    my @geneIDlist  = @{ $gene_adaptor->list_stable_ids() };
    return @geneIDlist;
  } # private_list_stable_ids_for_ensembl_database #

  method private_list_stable_ids_for_ensembl_database_through_job_handler (Ensembl_Database_Parameters $parameters) {
    ### this (private) method has been added as a quick-fix, it doesn't follow the usual way of starting with a generic method and testing for Genome_sequence_Database_Parameters type etc.
    my $object = Genome_DB_Utilities->new();
    my $function = 'private_list_stable_ids_for_ensembl_database';
    my @parameters = ($parameters);
    my $high_memory = TRUE;
    my @result = $GU->standard_call_to_job_handler($object,$function,\@parameters,$high_memory,TRUE);
    return \@result;
  } # private_list_stable_ids_for_ensembl_database_through_job_handler

  method private_get_genomic_sequence_of_ensembl_gi (Ensembl_Database_Parameters $parameters, Str $coord_sys_name, Str $region, Int $fiveprime, Int $threeprime, StrandType $strand, Str $repeat_masking) {
	
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      #my $alias = $parameters->alias; # redundant
      
      my $startposition = $fiveprime;
      my $endposition = $threeprime;
      
      if ($strand eq 'negative') {
	  $startposition = $threeprime;
	  $endposition = $fiveprime;
      }
      
      my $sequence;
#    my $masked_sequence;
      my $difference = ($endposition - $startposition);
      $GU->user_info(3, "diff is $difference\n");
      if (($endposition - $startposition) < 0) { 
	  # necessary if genomic interval created is empty or non-valid - return empty strings
	  $GU->user_info(3, "genomic interval is empty\n");
	  $sequence = '';
#      $masked_sequence = '';
      }
      else {

	  my $slice_adaptor = $self->private_get_ensembl_slice_adaptor($parameters);
	  my $slice_strand=1;
	  if ($strand eq 'negative') {
	      $slice_strand=-1;
	  }
	  
	  my $slice = $slice_adaptor->fetch_by_region($coord_sys_name, $region, $startposition, $endposition, $slice_strand);
          
	  if ($repeat_masking eq 'repeatmasked') {
	      my $hardmasked_sequence = $slice->get_repeatmasked_seq(); # changes type to Bio::EnsEMBL::RepeatMaskedSlice
	      $sequence = $hardmasked_sequence->seq();
	  }
	  else {
	      $sequence = $slice->seq();
	  }  
      }  
      return $sequence;
  } # private_get_genomic_sequence_of_ensembl_gi #

  method private_get_orthologous_ensembl_gene_IDs (Str $query_gene_id, Ensembl_Database_Parameters $parameters, Generic_Orthology_Finding_Parameters $ortholog_parameters) {
      my $alias = $parameters->alias;
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $db_adaptor = $self->registry->get_DBAdaptor($alias, 'core'); # get db name from alias
      my $species = $db_adaptor->species();
      $GU->user_info(2, "-----------------------------\nStarting to look for $species genes orthologous to ".$query_gene_id.". \n");
      if ( $ortholog_parameters->isa ('Compara_Orthology_Finding_Parameters' )) {
	  $GU->user_info(3, "Using Compara orthology\n");
	  # decide which compara DB to use
	  my $compara_db;
	  
	  if ($parameters->location eq 'ensembl') {
	      $GU->user_info(3, "Using remote ensembl DB\n");
	      $compara_db = 'Multi'; # this is the default for remote ensembl db
	  }
	  else {
	      $compara_db = $ortholog_parameters->compara_db;
	  }
	  #$GU->user_info(3,$compara_db."\n");# debugging
	  ## Get the compara member adaptor
	  my $compara_adaptor_member = $self->registry->get_adaptor($compara_db, 'compara', 'Member');
	  ## Get the compara homology adaptor
	  my $compara_adaptor_homology =  Bio::EnsEMBL::Registry->get_adaptor ($compara_db, 'compara', 'homology');
	  
	  my @targetspecieshomologues;
	  
	  my $qy_member = $compara_adaptor_member->fetch_by_source_stable_id( "ENSEMBLGENE",
									      $query_gene_id );
	  
	  
	  if (!$qy_member) {
	      $GU->user_info( 1, "qy_member stable id is undefined\n" ); # the next part is likely to crash if qy_member stable_id is undefined - have added a line to skip this request
	      $GU->user_info( 1, "skipping ".$query_gene_id."\n");
	      next; # exits the method before it crashes with 'undefined' error
	  }
	  else {
	      $GU->user_info(3,"Gene ENSEMBL ID: ".$qy_member->stable_id.", gene common name: ".$qy_member->display_label.". \n");
	      foreach my $homology (@{$compara_adaptor_homology->fetch_all_by_Member_paired_species( $qy_member, $species )}) { #was $targetspecies
		  if ( $homology->description ne 'between_species_paralog' ) {
		      foreach my $member_attribute (
			  @{ $homology->get_all_Member_Attribute } ) {
			  my ( $member, $attribute ) = @{$member_attribute};
			  next if ( $member->stable_id eq $qy_member->stable_id );
			  push @targetspecieshomologues, $member->stable_id;
			  #$GU->user_info( 3,"stable id: ". $member->stable_id."\n" );
			  $GU->user_info( 2,"Orthologous gene is found: ".$member->stable_id." ( ".$member->display_label." ). \n" );
		      }
		  }
	      }
	  }
	  my $number_of_orthologous_genes_found = scalar(@targetspecieshomologues);
	  if ($number_of_orthologous_genes_found eq 0) {
	      $GU->user_info( 2, "=======> No orthologous to ".$qy_member->stable_id." ( common name ".$qy_member->display_label." ) genes were found in ".$species." species!\n-----------------------------\n");
	  }
	  else {
	      $GU->user_info( 2,"The gene ".$qy_member->stable_id." ( common name ".$qy_member->display_label." ) has ".$number_of_orthologous_genes_found." orthologous gene(s) in ".$species.".\n-----------------------------\n" );
	  }
	  return @targetspecieshomologues; # return type: array
      } 
      elsif ( $ortholog_parameters->meta->name eq 'Synteny_Orthology_Finding_Parameters' ) {
	  die 'cannot find orthologous IDs using synteny method yet';
      }
      elsif ( $ortholog_parameters->meta->name eq 'RBH_Orthology_Finding_Parameters' ) {
	  die 'cannot find orthologous IDs using RBH method yet';
      }
      else {
	  die 'cannot determine orthologs for this method';
      }
  } # private_get_orthologous_ensembl_gene_IDs #
  
  method private_get_orthologous_gene_ID_using_RBH_method (Str $query_gene_id, Genome_Sequence_Database_Parameters $source_parameters, Genome_Sequence_Database_Parameters $target_parameters, Generic_Orthology_Finding_Parameters $ortholog_parameters) {
    # call RBH method via Job Handler to get ortholog mapping object
    $GU->user_info(3, "Getting RBH ortholog mapping between ".$source_parameters->dbname." and ".$target_parameters->dbname."\n");
    my $ortholog_mapping_maker = Ortholog_Mapping_Maker->new();
    my $ortholog_mapping_object = $ortholog_mapping_maker->create_ortholog_mappings_through_job_handler($source_parameters, $target_parameters);
    my $orthologous_gene_id = $ortholog_mapping_object->get_orthologous_ID_for_source_id($query_gene_id);
    $GU->user_info(3, "Reciprocal Best Hit ortholog: ".$orthologous_gene_id."\n");
    
    #return $orthologous_gene_id; # can we return a scalar if the function is expecting an array?
    my @result;

    if ($orthologous_gene_id eq 'not found') {
      $GU->user_info(3, "ortholog not found\n");
      return @result;
    }
    else {
      push (@result, $orthologous_gene_id);
      return @result; # return array containing orthologous gene id
    }
  } # private_get_orthologous_gene_ID_using_RBH_method #

  method private_get_orthologous_gene_ID_using_Random_Assignment_Method (Str $query_gene_id, Genome_Sequence_Database_Parameters $source_parameters, Genome_Sequence_Database_Parameters $target_parameters, Random_Assignment_Orthology_Finding_Parameters $ortholog_parameters) {
    my $seed = $ortholog_parameters->seed; # seed obtained from orthology finding parameters
    srand($seed); # uses fixed seed to enable caching (consistency of result)
    my @result;
    my $orthologous_gene_id;

    my @source_ids = $self->list_stable_ids_for_a_given_genome_database($source_parameters);
    @source_ids = sort (@source_ids);
    my $position;
    for (my $i = 0; $i < scalar(@source_ids); $i++ ) { # get the position of the query gene id in the sorted list of source ids
      if ($source_ids[$i] eq $query_gene_id) {
	$position = $i;
      }
    }
    my @target_ids = $self->list_stable_ids_for_a_given_genome_database($target_parameters);
    @target_ids = sort(@target_ids);
    my $random_number;
    my $range = @target_ids;# length of target_ids array
    # generate $position-th number of random numbers in range
    for (my $j = 0; $j < $position; $j++) {
      $random_number = int(rand($range)); # here
    }
    # use the last random number generated to select a random id from the target source list
    my $random_gene_id = $target_ids[$random_number];
    push (@result, $random_gene_id);
    $GU->user_info(2, $random_gene_id." has been randomly assigned as an ortholog to ".$query_gene_id."\n");
    return @result; # return array containing an 'orthologous' RANDOM gene id
  } # private_get_orthologous_gene_ID_using_Random_Assignment_Method #

  method private_get_all_chromosomes_ensembl (Ensembl_Database_Parameters $parameters) {
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $slice_adaptor = $self->private_get_ensembl_slice_adaptor($parameters);
      my @slices = (@{ $slice_adaptor->fetch_all('toplevel') } );
      my @genomic_intervals;
      my @mitochondrial_or_chloroplast_genome_names = ('Mt','Pt','C','M','MT');
      foreach my $slice (@slices) {
	  my $region = $slice->seq_region_name;
	  my @region_as_array = ($region);
	  if (!$GU->lists_overlap(\@mitochondrial_or_chloroplast_genome_names,\@region_as_array)) {  # exclude mitochondrial and chloroplast genomes	      
	      my $length = $slice->length;
	      $GU->user_info(3,"Region ".$region." has length ".$length.".\n");
	      my $coord_sys_name = $slice->coord_system()->name();
	      my $genomic_interval = Genomic_Interval->new(
		  'genome_db' => $parameters,
		  'coord_sys_name' => $coord_sys_name,
		  'region' => $region,
		  'five_prime_pos' => 1,
		  'three_prime_pos' => $length,
		  'strand' => 'positive',
		  'type' => 'dna',
		  ); 
	      push(@genomic_intervals,$genomic_interval);
	  } else {
	      $GU->user_info(2,"Excluded sequence region: ".$region."\n");
	  }
      }
      my $gis_maker = Genomic_Interval_Set_Maker->new();
      my $gi_set = $gis_maker->make_gi_set_from_gi_array(\@genomic_intervals);
      return $gi_set;
  } # private_get_all_chromosomes_ensembl #

  method private_ensembl_all_transcripts_to_fasta_file (Ensembl_Database_Parameters $parameters) { # was Str $dbname, APPLESSpeciesName $alias
    die 'method not re-implemented yet!';

    #add directory to store in, read from config file

    my $dbname = $parameters->dbname;
    my $outfile = Bio::SeqIO->new(-file => '>' . $dbname . '.all_transcripts.fasta',
				  -format => 'fasta');
    
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);
    my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
    foreach my $geneid (@{ $gene_adaptor->list_stable_ids() }) {
      my $gene = $gene_adaptor->fetch_by_stable_id($geneid);
      
      my @transcripts = @{ $gene->get_all_Transcripts };			
      foreach my $transcript (@transcripts) {
	my $splicedsequence = $transcript->spliced_seq();
	my $bioseq = Bio::Seq->new(-id => $transcript->display_id,
				   -seq => $splicedsequence);	   									   
	$outfile->write_seq($bioseq);
      }
      
    }
    return 1; # file written
  } # private_ensembl_all_transcripts_to_fasta_file #

  method private_ensembl_one_transcript_per_gene_to_fasta_file (Ensembl_Database_Parameters $parameters) { # was (Str $dbname, APPLESSpeciesName $alias)
    die 'method not re-implemented yet!';
    my $dbname = $parameters->dbname;
    my $alias = $parameters->alias;
    my $outfile = Bio::SeqIO->new(-file => '>' . $dbname . '.single_transcripts.fasta',
				  -format => 'fasta');
    
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);

    #my $gene_adaptor = $self->registry->get_adaptor( $dbname, 'core', 'gene' );
    my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
    
    foreach my $geneid (@{ $gene_adaptor->list_stable_ids() }) {
      my $gene = $gene_adaptor->fetch_by_stable_id($geneid);
      
      my @transcripts = @{ $gene->get_all_Transcripts };			
      my $transcript = $transcripts[0];
      my $splicedsequence = $transcript->spliced_seq();
      my $bioseq = Bio::Seq->new(-id => $transcript->display_id,
				 -seq => $splicedsequence);	   									   
      $outfile->write_seq($bioseq);
    }
    return 1;
  } # private_ensembl_one_transcript_per_gene_to_fasta_file #

  method private_ensembl_one_CDS_per_gene_to_fasta_file (Ensembl_Database_Parameters $parameters) { # was (Str $dbname, APPLESSpeciesName $alias)
    die 'method not re-implemented yet!';
    my $dbname = $parameters->dbname;
    my $alias = $parameters->alias;
    my $outfile = Bio::SeqIO->new(-file => '>' . $dbname . '.single_CDS.fasta',
				  -format => 'fasta');
   
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);

    #my $gene_adaptor = $self->registry->get_adaptor( $dbname, 'core', 'gene' );
    my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);

    foreach my $geneid (@{ $gene_adaptor->list_stable_ids() }) {
      my $gene = $gene_adaptor->fetch_by_stable_id($geneid);
      
      my @transcripts = @{ $gene->get_all_Transcripts };			
      my $transcript = $transcripts[0];
      my $cds = $transcript->translateable_seq();
      if ($cds eq '') {
	$GU->user_info( 1, "no translateable sequence for ". $transcript->display_id."!\n" );
      }
      else {
	my $bioseq = Bio::Seq->new(-id => $transcript->display_id,
				   -seq => $cds);	
	$outfile->write_seq($bioseq);
      }
    }
    return 1;
  } # private_ensembl_one_CDS_per_gene_to_fasta_file #

  method private_ensembl_genome_to_fasta_file (Ensembl_Database_Parameters $parameters) { # was (Str $dbname, APPLESSpeciesName $alias)
    die 'method not re-implemented yet!';
    my $dbname = $parameters->dbname;
    my $outfile = Bio::SeqIO->new(-file => '>'. $dbname. '.genome.fasta',
				  -format => 'fasta');
    
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);

    my $cs_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'CoordSystem' );
    
    # List all coordinate systems in this database:
    my @coord_systems = @{ $cs_adaptor->fetch_all() };
    foreach my $cs (@coord_systems) {
      $GU->user_info(3,"Coordinate system: ".$cs->name()." ".$cs->version."\n");
    }
    
    my $slice_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'Slice' );
    
    # Get all slices on the highest coordinate system:
    foreach my $slice (@{ $slice_adaptor->fetch_all('toplevel') } ) {
      $outfile->write_seq($slice);
    }	
    return 1;
  } # private_ensembl_genome_to_fasta_file #

  method private_ensembl_masked_genome_to_fasta_file (Ensembl_Database_Parameters $parameters) { # was (Str $dbname, APPLESSpeciesName $alias)
    die 'method not re-implemented yet!';
    my $dbname = $parameters->dbname;
    my $outfile = Bio::SeqIO->new(-file => '>'. $dbname. '.masked_genome.fasta',
				  -format => 'fasta');
   
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);

    my $cs_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'CoordSystem' );
    
    # List all coordinate systems in this database:
    my @coord_systems = @{ $cs_adaptor->fetch_all() };
    foreach my $cs (@coord_systems) {
      $GU->user_info(3,"Coordinate system: ".$cs->name()." ".$cs->version."\n");
    }
    
    my $slice_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'Slice' );
    
    # Get all slices on the highest coordinate system:
    foreach my $slice (@{ $slice_adaptor->fetch_all('toplevel') } ) {
      my $maskedseq = $slice->get_repeatmasked_seq;
      $outfile->write_seq($maskedseq);
    }	
  } # private_ensembl_masked_genome_to_fasta_file #	

  method private_ensembl_database_stats (Ensembl_Database_Parameters $parameters) {
      my $dbname = $parameters->dbname;
      $GU->user_info( 2, "\nEnsembl Genome Statistics for ". $dbname."\n" );
      
      $self->private_connect_to_registry($parameters->location, $dbname);
      
      my $cs_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'CoordSystem' );
      
      my @coord_systems = @{ $cs_adaptor->fetch_all() };
      foreach my $cs (@coord_systems) {
	  $GU->user_info(2,"Coordinate system: ".$cs->name()." ".$cs->version."\n");
      }
      
      my $slice_adaptor = $self->registry->get_adaptor( $dbname, 'Core', 'Slice' );
      my $counter;
      foreach my $slice (@{ $slice_adaptor->fetch_all('toplevel') } ) {
	  $counter++;
      }
      $GU->user_info( 2, "Toplevel slices: ". $counter."\n" );
      
      #my $gene_adaptor = $self->registry->get_adaptor( $dbname, 'core', 'gene' );
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);

      my $genecounter;			
      foreach my $geneid (@{ $gene_adaptor->list_stable_ids() }) {
	  $genecounter++;
      }
      $GU->user_info( 2, "Genes: ".$genecounter."\n" );
  } # private_ensembl_database_stats #

  method private_ensembl_api_and_database_consistent_versions_check (Str $dbname) {# (APPLESSpeciesName $alias) {
    # does this work with dbname and/or alias?
    # private method, only to be used with Ensembl databases
    # uses ensembl's version_check method to check if database version and api version match
    # a mismatch results in 'die'; a match returns 1.
    my $registry = 'Bio::EnsEMBL::Registry';
    my $ensembl_registry_conf = "/home/grannysmith/webseaweeds/Configuration/ensembl_registry_voland.conf"; # environmental variable (set in $GU->load_includes, from value specified in APPLES.dat config file)
    $registry->load_all($ensembl_registry_conf);
    
    # die if api/database versions mismatch

    my $API_version = $registry->software_version();
    my $check;
    if ($API_version < 55) { # at which API version did get_all_DBAdaptors_by_dbname method not exist?
	my $db_API_version = $self->private_get_api_version_from_db_name($dbname);
        

	$check = ($db_API_version == $API_version);
    } else {
	my @db_adaptors = @{ $registry->get_all_DBAdaptors_by_dbname($dbname) };
	$check = $registry->version_check($db_adaptors[0]);
    }
    if (!$check) {
      die 'mismatch between API version and ensembl database version';
    }
    else {
      return 1;
    }  
  } # private_ensembl_api_and_database_consistent_versions_check #

  method private_get_api_version_from_db_name (Str $dbname) {
      # assumptions: a) "_core_" is part of the database name and precedes the API version number
      #              b) API versions have two digits
    
      my $result;
      my $index = index($dbname, "_core_");
            
      my $start_position = $index+6;
      $result = substr($dbname,$start_position,2);
      return $result;
  } # private_get_api_version_from_db_name #

  method private_ensembl_proteins_to_fasta_file (Ensembl_Database_Parameters $parameters, Str $filename, TranscriptChoice $transcript_choice) {
    $GU->user_info(1,"XXX making fasta file from ".$parameters->dbname."\n");
    my $outfile = Bio::SeqIO->new(-file =>  '>'.$filename,
				  -format => 'fasta');
    my $dbname = $parameters->dbname;
    my $alias = $parameters->alias;
    
    $self->private_connect_to_registry($parameters->location, $parameters->dbname);

    my $gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
    my @gene_ids = $self->list_stable_ids_for_a_given_genome_database ($parameters); 
    #@gene_ids = @gene_ids[0..20]; # VOLATILE
    foreach my $geneid (@gene_ids) { 
      
      # start - testing of Running_Perseverance here (around "my $gene = $gene_adaptor->fetch_by_stable_id($geneid);")
      my $gene;
      my $reconnect = FALSE;
      while (!$main::global_perseverance->stop_trying) {
	eval {
	    if ($reconnect) {
		# attempt to refresh the connection with ensembl
		$GU->user_info(3, "attempting to reconnect to ensembl\n");
		$self->private_reconnect_to_registry($parameters->location, $parameters->dbname);
		$gene_adaptor = $self->private_get_ensembl_gene_adaptor ($parameters);
		$GU->user_info(3, "successful reconnect?\n");
		$reconnect = FALSE;
	    }
	    $gene = $gene_adaptor->fetch_by_stable_id($geneid);
	    my @transcripts = @{ $gene->get_all_Transcripts };
	    my @transcripts_to_translate;
	    my $gene_external_id;
	    if ($transcript_choice eq 'one') {
		$gene_external_id = $gene->external_name();
		if (!defined $gene_external_id) {
		    $gene_external_id = '';
		}
		$GU->user_info(3, "gene external id: ".$gene_external_id."\n");
		my $transcript = $self->private_select_one_ensembl_transcript (\@transcripts, $gene_external_id);
		$GU->user_info (3, "transcript: ".$transcript."\n");
		push (@transcripts_to_translate, $transcript);
	    }
	    
	    elsif ($transcript_choice eq 'all') {
		@transcripts_to_translate = @transcripts
	    }
	    
	    else {
		die 'method to deal with this TranscriptChoice not implemented yet!';
	    }	
	    foreach my $transcript_to_translate (@transcripts_to_translate) {
		my $cds = $transcript_to_translate->translateable_seq();
		if ($cds eq '') {
		    $GU->user_info( 3, "no translateable sequence for ". $transcript_to_translate->display_id."!\n" );
		}
		else {
		    $GU->user_info( 3, "translateable sequence exists for ". $transcript_to_translate->display_id."\n" );
		    my $protein = $transcript_to_translate->translate();
		    my $bioseq = Bio::Seq->new(-id => $geneid,
					       -seq => $protein->seq());
		    $outfile->write_seq($bioseq); # not using '-id => $transcript_to_translate->display_id' here any more, because we want to retain the original gene id for the purposes of ortholog mapping - what about in other cases though? add another parameter to use gene id/display id?
		}
	    }
	}; # end of eval statement
	my $error = $@;
	$main::global_perseverance->decide_on_rerun('Ensembl', FALSE, $error);
	# specific-case error handling code here: (optional)
	if ($error) { # something has gone wrong...
	    $reconnect = TRUE;	   
	}
      } # end of while loop - if there was an error we (may) try again
      # and reset the try counter and stop_trying attribute
      $main::global_perseverance->stop_trying(FALSE);
      # end - testing of Running_Perseverance
    } # end of foreach geneid loop
    $GU->user_info(3,$filename." file written\n");
    chmod (0777, $filename);
    return 1;
  } # private_ensembl_proteins_to_fasta_file #
  
  method private_select_one_ensembl_transcript (ArrayRef $transcripts_ref, Str $gene_external_name) {
    my $selected_transcript = $$transcripts_ref[0];
    foreach my $transcript (@{$transcripts_ref}) {
      if ($gene_external_name eq $transcript->external_name) {
	$selected_transcript = $transcript;
	$GU->user_info(3, "selected transcript: $selected_transcript.\n");
      }
    }
    return $selected_transcript;
  } # private_select_one_ensembl_transcript #

  method private_get_ensembl_gene_adaptor (Ensembl_Database_Parameters $parameters) {
    my $API_version = $self->registry->software_version();
    #$GU->user_info(1,"API version: ". $API_version."\n"); # debugging
    my $gene_adaptor;
    if ($API_version < 55) { # at which API version did get_all_DBAdaptors_by_dbname method not exist?
      $gene_adaptor = $self->registry->get_adaptor( $parameters->alias, 'core', 'gene' );
    }
    else {
      my @db_adaptors = @{ $self->registry->get_all_DBAdaptors_by_dbname($parameters->dbname) };
   
      my $species = $db_adaptors[0]->species(); 
      $gene_adaptor = $self->registry->get_adaptor( $species, 'core', 'gene' );
    }
    if (!defined $gene_adaptor) { 
      die 'gene adaptor is undefined: possibly there is a mismatch between the release versions of your installed API and the requested database. Please ensure they match\n';
    }
    return $gene_adaptor;
  } # private_get_ensembl_gene_adaptor #

  method private_get_ensembl_slice_adaptor (Ensembl_Database_Parameters $parameters) {

    my $API_version = $self->registry->software_version();
    my $slice_adaptor;
    if ($API_version < 55) { # at which API version did get_all_DBAdaptors_by_dbname method not exist?
      $slice_adaptor = $self->registry->get_adaptor($parameters->alias, 'core', 'slice');
    }
    else {
      my @db_adaptors = @{ $self->registry->get_all_DBAdaptors_by_dbname($parameters->dbname) };
      my $species = $db_adaptors[0]->species(); 
      $slice_adaptor = $self->registry->get_adaptor( $species, 'core', 'slice' );
        
    }
    if (!defined $slice_adaptor) { 
      die 'slice adaptor is undefined: possibly there is a mismatch between the release versions of your installed API and the requested database. Please ensure they match\n';
    }
    return $slice_adaptor;
  } # private_get_ensembl_slice_adaptor #

  method private_get_coord_system_name_from_ensembl_gene_id (Ensembl_Database_Parameters $parameters, Str $geneid) {  	
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $slice_adaptor = $self->private_get_ensembl_slice_adaptor($parameters);
      my $slice = $slice_adaptor->fetch_by_gene_stable_id($geneid);
      my $coord_sys  = $slice->coord_system()->name();
      return $coord_sys;
  } # private_get_coord_system_name_from_ensembl_gene_id #

  method private_get_seq_region_name_from_ensembl_gene_id (Ensembl_Database_Parameters $parameters, Str $geneid) {
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor($parameters);
      my $geneobject = $gene_adaptor->fetch_by_stable_id($geneid);
      my $region = $geneobject->seq_region_name;
      return $region;
  } # private_get_seq_region_name_from_ensembl_gene_id #

  method private_get_ensembl_gene_start_position (Ensembl_Database_Parameters $parameters, Str $geneid) {
      my $strand = $self->private_get_strand_of_ensembl_gene($parameters,$geneid);
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor($parameters);
      my $geneobject = $gene_adaptor->fetch_by_stable_id($geneid);
      my $position = $geneobject->start;
      if ($strand eq 'negative') {
	  $position = $geneobject->end;
      }
      return $position;
  } # private_get_ensembl_gene_start_position #

  method private_get_ensembl_gene_end_position (Ensembl_Database_Parameters $parameters, Str $geneid) {
      my $strand = $self->private_get_strand_of_ensembl_gene($parameters,$geneid);
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor($parameters);
      my $geneobject = $gene_adaptor->fetch_by_stable_id($geneid);
      my $position = $geneobject->end;
      if ($strand eq 'negative') {
	  $position = $geneobject->start;
      }
      return $position;
  } # private_get_ensembl_gene_end_position #

  method private_get_strand_of_ensembl_gene (Ensembl_Database_Parameters $parameters, Str $geneid) {
      $self->private_connect_to_registry($parameters->location, $parameters->dbname);
      my $gene_adaptor = $self->private_get_ensembl_gene_adaptor($parameters);
      my $geneobject = $gene_adaptor->fetch_by_stable_id($geneid);
      my $strand;
      if ($geneobject->strand == 1) {
	  $strand = 'positive';
      }
      elsif ($geneobject->strand == -1) {
	  $strand = 'negative';
      }
      else {
	  die 'ensembl strand has an incorrect value!';
      }
      return $strand;
  } # private_get_strand_of_ensembl_gene #

  method private_ensembl_protein_intervals_to_fasta_file (Genome_Sequence_Database_Parameters $parameters, Str $filename, TranscriptChoice $transcript_choice) {
    die 'method not implemented yet!';
  } # private_ensembl_protein_intervals_to_fasta_file #

############################################
##### GenBank-specific private methods #####
############################################

method private_get_distance_to_neighbouring_genbank_gene(Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
    my $result = $self->private_get_distance_or_id_neighbouring_genbank_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction, TRUE); 
    return $result;
  } # private_get_distance_to_neighbouring_genbank_gene #

method private_get_gene_id_of_neighbouring_genbank_gene (Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction) {
    my $result = $self->private_get_distance_or_id_neighbouring_genbank_gene($reg_loc, $ignore_pseudogenes, $ignore_RNAgenes, $direction, FALSE);
    return $result;
  } # private_get_gene_id_of_neighbouring_genbank_gene # 

method private_get_distance_or_id_neighbouring_genbank_gene (Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes, GenomicDirection $direction, Boolean $get_distance) {
	
    my $parameters = $reg_loc->genome_db;
    my $position = $reg_loc->position;
    if ($reg_loc->gene_ID) {
	if ($reg_loc->strand eq 'negative') {
	    my $geneinf = $self->private_get_genbank_geneinfo($reg_loc->gene_ID);
	    
	    $position = $geneinf->{END};
	}
    }

    my $geneinf = $self->private_get_genbank_geneinfo($reg_loc->gene_ID);
 
    my $gene_id_to_ignore;
    if ($reg_loc->gene_ID) {  
      $gene_id_to_ignore = $reg_loc->gene_ID;
    } 
    else {
      $gene_id_to_ignore = '';
    }

    my @neighbouring_gene_information = $self->private_get_neighbouring_gene_information($reg_loc->genome_db,
											 $reg_loc->region,
											 $position,
											 $ignore_pseudogenes,
											 $ignore_RNAgenes, # [HM] introduced in May 2010
											 $gene_id_to_ignore);
    my $relevant_index = 1;
    if ($direction eq 'towards_three_prime') {
      $relevant_index = 3;
    }
    my $relevant_item = $neighbouring_gene_information[$relevant_index];
    if ($relevant_item eq 'NONE') {
	my $exception = No_Neighbouring_Gene_Exception->new();
	die $exception;
    }
    if ($get_distance) {
      return $relevant_item;
    } 
    else {
      my $id = $neighbouring_gene_information[$relevant_index-1];
      return $id;
    }
  } # private_get_distance_or_id_neighbouring_genbank_gene #

method private_get_id_and_distance_of_neighbouring_genbank_gene (Reg_Loc $reg_loc, Boolean $ignore_pseudogenes, Boolean $ignore_RNAgenes) {
    my $gene_id_to_ignore = $reg_loc->gene_ID;

    my @neighbouring_gene_info = $self->private_get_neighbouring_gene_information($reg_loc->genome_db,
										  $reg_loc->region,
										  $reg_loc->position,
										  $ignore_pseudogenes,
										  $ignore_RNAgenes, # [HM] introduced in May 2010
										  $gene_id_to_ignore);
    return @neighbouring_gene_info;
  } #  private_get_id_and_distance_of_neighbouring_genbank_gene # 

  method private_fetch_genbank_file (Str $version) {
    #	A method for fetching a single genbank file from ncbi
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    my $rawfile = $data_dumps_dir.'raw/'.$version . '.raw';
    if (!-e $rawfile) { # if file doesn't already exist, fetch it
    	mkdir ($data_dumps_dir.'raw');
    	chmod (0777,$data_dumps_dir.'raw');
		my $OK = FALSE;
		while ($OK == FALSE) {
			eval{
				#	In fact there 
				my $genObj = $gb->get_Stream_by_version($version);
  				my $seq = $genObj->next_seq();
  				if (!defined($seq)) {
  				  # We have tried downloading a non-existant file
  				  open (my $f, '>'.$rawfile);
  				  close $f;
  				}
  				else {
        	                   my $seqOut = new Bio::SeqIO(-file => '>'.$rawfile, -format => 'genbank');
        	                   $seqOut->write_seq($seq);
  				}
				$OK = TRUE;
			};
			if ($@) {
			    $gb = new Bio::DB::GenBank;
			    $GU->user_info( 3, "Error obtaining genbank file.  Trying again with new i/f\n" );
			    sleep(3);
			}
		}
    }
    if (!-e $rawfile) { # if file doesn't exist now, die
      die 'could not write file\n';
    }
    return $rawfile;
  } # private_fetch_genbank_file # 

  method private_fetch_genbank_files (ArrayRef $files) {
    #	A method for fecthing a set of genbank files from ncbi
    #	More efficient than fetching one at a time
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    mkdir($data_dumps_dir.'raw');
    chmod(0777,$data_dumps_dir.'raw');
    my @fetchFiles;
    my @fetchFilenames;
    my @rawfiles;
    #	See which rawfiles already exist, and create a list of files that don't and need 
    #	to be collected
    foreach my $file (@{$files}) {
	my $rawfile = $data_dumps_dir.'raw/'.$file . '.raw';
	push(@rawfiles,$rawfile);
    	if (!-e $rawfile) {
	    push(@fetchFiles,$file);
	    push (@fetchFilenames,$rawfile);
    	}
    }
    if (scalar(@fetchFiles) > 0) {
	my $OK = 0;
	while ($OK == 0) {
	    eval {
		my $count = 0;
		#	Oddly documentation says that get_Stream_by_version method does not exist, but it does
		#	and skips the delay that is used in get_Seq_by_Version.  The processing of the
		#	data within Perl adds sufficient delay that the additional delay is largely unecessary
		my $genObj = $gb->get_Stream_by_version(\@fetchFiles);
		while( my $seq = $genObj->next_seq() ) {
		    my $seqOut = new Bio::SeqIO(-file => '>'.$fetchFilenames[$count], -format => 'genbank');
		    $seqOut->write_seq($seq);
		    if (!-e $fetchFilenames[$count]) { # if file doesn't exist now, die
			die 'could not write file\n';
		    }
		    $count++;
		}
		$OK = 1;
	    };
	    if ($@) {
		$gb = new Bio::DB::GenBank;
		$GU->user_info( 3, "Error obtaining genbank file.  Trying again with new i/f\n" );
		sleep(3);
	    }
	}
    }
    return @rawfiles;
  } # private_fetch_genbank_files # 

  method private_get_all_genbank_gene_positions (Str $version) {
    # this method should always be called via the Job_Handler for efficiency, using the handle_APPLES_function method
    # have split into retrieval of raw data (private_fetch_genbank_file) and creation of geneinfdump
    my @gene_position_data;
    my $region_name = $version;
#    my @translations;
    my @proteins;
    my $chromosome = '0';
    my $coord_system  = '';
    my $start;
    my $end;
    my $strand;
    my $typeofgene;
    my $gene = '';
    my $locus = '';
    
  	sub get_db_xref {
  		#	Dont try stepping through this code using Eclipse, 
  		#	Perl just exists with exit code 22.
		my ($self, $featObj) = @_;
		my $xref = '';
		my @values = $$featObj->get_tag_values('db_xref');
		foreach my $value (@values) {
		  if (index($value,'TAIR:') == 0) {
		    $xref = $value;
		    last; 
		  }
		  if (index($value,'EcoGene') == 0) {
		    $xref = $value;
		    last; 
		  }
		  if (index($value,'UniProt') == 0) {
		    $xref = $value;
		    last; 
		  }
		}
		if ($xref eq '') {
			$xref = $values[0];
			}
		return $xref;
	}
	
	sub make_name_usable(){
		my ($self, $gene) = @_;
		#	This is going to be used as a filename, s make sure it has no filename incompatible 
		#	characters, and also does not match any of the reserved Windows filenames
		$gene =~ s/[\:,\\,\/,\<,\>,\|]/-/g;
		#
		#	And convert spaces into underscores to make it easier for unix, and so
		#	that it can be used as an identifier in Fasta files
		$gene =~ s/ /_/g;
		#	\A and \Z match the start and end, so must be the whole string
		#	\d is a single digit (deals with lpt1 through 9)
		#	$& is the matched string
		#	/i at the end makes it case insensitive
		$gene =~ s/\A(con|prn|nul|aux|com\d|lpt\d)\Z/$&_1/i;
		
		return $gene;
	}
	
  
    
    my $gene_position_stored = TRUE;	#So we dont store intial garbage
    my %locusMap;		#A map from locus to genes
    my %geneMap;
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    
    mkdir($data_dumps_dir.'inf');
    chmod(0777,$data_dumps_dir.'inf');
    mkdir($data_dumps_dir.'gene');
    chmod(0777,$data_dumps_dir.'gene');
    mkdir($data_dumps_dir.'gene/'.$region_name);
    chmod(0777,$data_dumps_dir.'gene/'.$region_name);

    my $rawfile = $self->private_fetch_genbank_file($version);    
    my $seqio_object = Bio::SeqIO->new(-file => $rawfile,-format => 'genbank');
    my $seqObj = $seqio_object->next_seq();
    my $first = TRUE; 
    foreach my $featObj ($seqObj->get_SeqFeatures()) {
      if ($featObj->primary_tag() eq 'source') {
      	if ($featObj ->has_tag('chromosome')) {
      		$coord_system = 'chromosome';
	      	$chromosome = ($featObj->get_tag_values('chromosome'))[0];
      	}
      	elsif ($featObj ->has_tag('mol_type')) {
      		my @values = $featObj->get_tag_values('mol_type');
	      	if ($values[0] eq 'genomic DNA') {
	      		$coord_system = 'chromosome';
	      	}
	      	else{
	      		$coord_system = $values[0];
	      	}
      	}
      }
    ###########################
    #	Check the 'gene' primary region
    elsif ($featObj->primary_tag() eq 'gene') {
		if ($locus ne '') {
			if (not $gene_position_stored) {
	      			#	Prepare the data that will be cached
					my $entry = {
					     ID         => $gene,
					     CHROMOSOME => $region_name, 
					     # might also be the ID of a scaffold (or other type of sequence region)
					     START      => $start,
					     END        => $end,
					     TYPE       => $typeofgene,
					    };
					if ($first == TRUE) {
						$entry -> {VERSION} = DBVERSION;
						$first = FALSE;
					}
					push( @gene_position_data, $entry );
      			}
  			$gene_position_stored = FALSE;

   			#	And the data that will be stored in files
			my $geneInf = {
			     ID           => $gene,
			     COORD_SYSTEM => $coord_system,
			     REGION       => $region_name, 
			     CHROMOSOME   => $chromosome, 
			     START        => int($start),
			     END          => int($end),
			     STRAND       => $strand,
			     TYPE         => $typeofgene,
			     PROTEINS         => \@proteins,
				};
	        my $genfile = $data_dumps_dir.'gene/'.$gene.'.geneInf';
			nstore (\$geneInf,$genfile);
			chmod(0664,$genfile);

			$locusMap{$locus} = $gene;
			$geneMap{$gene} = 1;
			
			@proteins = ();
      		}

		$gene = '';      	
		$start      = $featObj->start();
		$end        = $featObj->end();
		$strand 	   = ($featObj->strand() == 1)?'positive':'negative';
		$typeofgene = '';
		if ($featObj->has_tag('locus_tag')) {
		    $locus = ($featObj->get_tag_values('locus_tag'))[0];
			}
		#	First preference for the gene name is the gene tag
		if ($featObj->has_tag('gene')) {
		    $gene = ($featObj->get_tag_values('gene'))[0];
			}
		#	Otherwise it is the cross reference to other identifiers, 
		#	with some being preferred
		elsif ($featObj ->has_tag('db_xref')) {
			$gene = $self -> get_db_xref(\$featObj);
		}
		#	Pseudo gene identifier
		elsif ($featObj->has_tag('pseudo'))	{
		    $typeofgene = 'pseudo';
		}
			
		#	If all else fails use the locus, which is really just internal to the genbank file
		#	and not always present
		if ($gene eq '') {
		    $gene = $locus;
		}
		
		$gene = $region_name.'/'.$self -> make_name_usable($gene);
		
		#	Now make sure we have not seen it before	
		if (defined($geneMap{$gene})) {
		    my $i = 2;
		    while (defined($geneMap{$gene.'_'.$i})){$i++;};
		    $gene = $gene.'_'.$i;			
		}
			
		if ($locus eq '') {
		    $locus = $gene;
		}
      }
	###########################
    #	Check the 'ncRNA' primary region
    elsif ($featObj->primary_tag() eq 'ncRNA') {
		$typeofgene = ($featObj->get_tag_values('ncRNA_class'))[0];
      	}
	###########################
    #	Check the 'CDS' primary region
	elsif ($featObj->primary_tag() eq 'CDS') {
		my $CDSlocus;
	  	if ($featObj ->has_tag('locus_tag')) {
	    	$CDSlocus = ($featObj->get_tag_values('locus_tag'))[0];
	  	}
	  	elsif ($featObj ->has_tag('gene'))	{
			$CDSlocus = $region_name.'/'.($featObj->get_tag_values('gene'))[0];
	  	}
		else {
	    	$GU->user_info(2, "CDS without a gene or a locus tag\n");
		}
		my $CDSprotein;
	    if ($locus ne $CDSlocus) {
			if (not $gene_position_stored) {
      			
		  #	Store current gene data, as we are about to replace it
		  my $entry = {
		      ID         => $gene,
		      CHROMOSOME => $region_name, 
		      # might also be the ID of a scaffold (or other type of sequence region)
		      START      => $start,
		      END        => $end,
		      TYPE       => $typeofgene,
		  };
		  if ($first == TRUE){
		      $entry -> {VERSION} = DBVERSION;
		      $first = FALSE;
	    }
		  push( @gene_position_data, $entry );
	      }
	      
	      #	First of all check to see whether we might have previously seen the protein
		my $oldgene = $locusMap{$CDSlocus};
	      
	    if (defined($oldgene)) {
			#	Store away the current gene info, because we are going to replace it with
			#	the previously stored data.
			my $geneInf = {
				ID           => $gene,
				COORD_SYSTEM => $coord_system,
				REGION       => $region_name, 
				CHROMOSOME   => $chromosome, 
				START        => int($start),
				END          => int($end),
				STRAND       => $strand,
				TYPE         => $typeofgene,
				PROTEINS     => \@proteins,
				};
			my $genfile = $data_dumps_dir.'gene/'.$gene.'.geneInf';
			nstore (\$geneInf,$genfile);
			chmod(0664,$genfile);
	  
			$locusMap{$locus} = $gene;

			#	Get the previous data								
			my $geneinf	= $self -> private_get_genbank_geneinfo($gene);
			$start = $geneinf ->{START};
			$end = $geneinf ->{END};
			$strand = $geneinf ->{STRAND};
			$typeofgene = $geneinf ->{TYPE};
			@proteins = @{$geneinf ->{PROTEINS}};
				
			#	We are reverting to the previously stored data, so dont store it again
			$gene_position_stored = TRUE;
			  
			$locus = $CDSlocus;
	      }
	      else {
			  #	And the data that will be stored in files
			  if ($gene ne '') {
			      my $geneInf = {
				  ID           => $gene,
				  COORD_SYSTEM => $coord_system,
				  REGION       => $region_name, 
				  CHROMOSOME   => $chromosome, 
				  START        => int($start),
				  END          => int($end),
				  STRAND       => $strand,
				  TYPE         => $typeofgene,
				  PROTEINS         => \@proteins,
			      };
		      my $genfile = $data_dumps_dir.'gene/'.$gene.'.geneInf';
		      nstore (\$geneInf,$genfile);
		      chmod(0664,$genfile);
		      
		      $locusMap{$locus} = $gene;
		      $geneMap{$gene} = 1;
			  }
				
		  $GU->user_info(2, 'Orphan CDS: Gene locus = '.$locus.' CDS locus = :'. $CDSlocus."\n");
		  $locus = $CDSlocus;
		  $gene = '';
		  $gene_position_stored = FALSE;
		  @proteins = ();
	      }
	  }
			
	  my $CDSstart      = $featObj->start();
	  my $CDSend        = $featObj->end();
	  my $CDSstrand 	   = ($featObj->strand() == 1)?'positive':'negative';
	  my $trans = '';	# zero length string indicates that the translation is not available
	  
	  if ($gene eq '') {
	      #	We have a locus without an associated gene bit, eg in the files that
	      #	have been converted
	      $gene = $locus;
	      $gene =~ s/\s.*//g;			#  Remove trailing garbage
	      
	      $gene = $region_name.'/'.$self -> make_name_usable($gene);
	      
	      #	Get gene location from the CDS
	      $start = $CDSstart;
	      $end = $CDSend;
	      $strand = $CDSstrand;
	  }
	  if ($featObj->has_tag('pseudo'))	{
	      $typeofgene = 'pseudo';
	  }

	  if ($featObj->has_tag('protein_id')) {
	      $typeofgene = 'protein_coding';
	      $CDSprotein = ($featObj->get_tag_values('protein_id'))[0];
	      $trans = undef;	# We can get the translation from the genbank file later
	  }
	  elsif ($featObj->has_tag('db_xref')) {
	      if ($typeofgene eq '') {
			  $typeofgene = 'unknown';
	      }
	      $CDSprotein = $self ->get_db_xref(\$featObj);
	  }
	  else {
	      if ($typeofgene eq '') {
		  $typeofgene = 'unknown';
	      }
	      $CDSprotein = "unknown";
	  }

	  $CDSprotein =~ s/\s.*//g;	#Get rid of trailling garbage that might be put in by local genbank creator
	  
	  my $CDSproteinFilename = $self -> make_name_usable($CDSprotein);

	  push(@proteins,$CDSproteinFilename);	#So we can retreive it from the filename
	
	  if ($featObj->has_tag('translation')) {
	      $trans = ($featObj->get_tag_values('translation'))[0];
	      $trans =~ s/\s//g;
		  }

	  my $cdsInf = {
	      ID           => $CDSprotein,	#The actual CDS name, rather than what it had to be changed to
	      GENE		  => $gene,
	      COORD_SYSTEM => $coord_system,
	      REGION       => $region_name, 
	      CHROMOSOME   => $chromosome, 
	      START        => int($CDSstart),
	      END          => int($CDSend),
	      STRAND       => $CDSstrand,
	      TYPE         => $typeofgene,
	      TRANSLATION => $trans,
	  };
	  my $cdsInfName = $data_dumps_dir.'inf/'.$CDSproteinFilename.'.cdsInf';
	  nstore (\$cdsInf,$cdsInfName);
	  chmod(0664,$cdsInfName);
      }
    }
	
    #	Finshed the file, store any data that is still current
    if ($gene ne '') {
		my $geneInf = {
		    ID           => $gene,
		    COORD_SYSTEM => $coord_system,
		    REGION       => $region_name, 
		    CHROMOSOME   => $chromosome, 
		    START        => int($start),
		    END          => int($end),
		    STRAND       => $strand,
		    TYPE         => $typeofgene,
		    PROTEINS     => \@proteins,
		};
	    
	    my $genfile = $data_dumps_dir.'gene/'.$gene.'.geneInf';
		nstore (\$geneInf,$genfile);
		chmod(0664,$genfile);
		
		if (!$gene_position_stored){
			my $entry = {
			    ID         => $gene,
			    CHROMOSOME => $region_name, 
			    # might also be the ID of a scaffold (or other type of sequence region)
			    START      => $start,
			    END        => $end,
			    TYPE       => $typeofgene,
				};
			push( @gene_position_data, $entry );
		}
    }
    $GU->user_info(2, scalar(@gene_position_data).' genes found\n');
    return @gene_position_data;
  } # private_get_all_genbank_gene_positions #

  method private_list_stable_ids_for_genbank_database (Genbank_Sequence_Database_Parameters $parameters) {
    my @all_gene_positions;
    
    $GU->user_info(3, "getting all genbank gene positions via job handler and cache\n");
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    my $cache = TRUE;
    my $job_handler = Job_Handler->new();
    my @geneIDlist;
    
    foreach my $file (@{ $parameters->filenames}) {
#    	my $recache = TRUE;
    	my $recache = FORCE_REFRESH;
    	my $dataOK = FALSE;
    	
    	do {
		    $job_handler->get_config_settings();
		    push (my @parameters, $file);
		    my $function = 'private_get_all_genbank_gene_positions';
		    my $job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
							     wall_time_estimate => 172800,
							    );
			if ($recache == TRUE) {
				$job_parameters ->{recache} = TRUE;
			}							    
							    
		    my $GDBU = Genome_DB_Utilities->new();
		    my @cache_result = eval { 
		      $job_handler->handle_APPLES_function($function, $GDBU, \@parameters, $cache, 365, $job_parameters); 
		    };
		    
		    my $exception_info = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler
		    if ($@) {
			my $exception_content = $@;
			if (!UNIVERSAL::can($exception_content, 'isa')) {   
			    die $exception_content; # throw error string
			} else {  
			    if ($exception_content->isa('Job_Information_Exception')) {
				$exception_info->merge($exception_content);
				$GU->user_info(1,"this function (private_list_stable_ids_for_genbank_database) can not be completed until result is available\n");
				$GU->user_info(1,"throwing an exception!\n");
				die $exception_info;
			    }
			    else {
				die $exception_content; # throw error object
			    }
			} 
		    }
		    # otherwise (no errors caught), handle the returned result (in this case we are not interested in the calling object):
		    my $object = shift(@cache_result); # first element of array is the object - must shift it off the array
		    @all_gene_positions = @cache_result; # the rest of the array is now the result
		    #$GU->user_info(3, "all gene positions\n"); # for debugging
		    #$GU->user_info(3, Dumper (\@all_gene_positions)); # for debugging
		
			#	Need to check that the genes are in the local directory
			#	They may not be because someone else may have put the results into the cache.  If so then we have to put the
			#	data into the directory
			
			if ((scalar(@all_gene_positions) == 0) || 
			    (($all_gene_positions[0] ->{CHROMOSOME} eq $file) && 
			     (defined($all_gene_positions[0] ->{VERSION})) && 
			     ($all_gene_positions[0] ->{VERSION} == DBVERSION))) {
			    my $genfile = $data_dumps_dir.'gene/'.$all_gene_positions[0]->{ID}.'.geneInf';
				if (!-e $genfile) {
					# Call the generator function but throw away the result
					$self -> private_get_all_genbank_gene_positions($file);
				}
				$dataOK = TRUE;
			}
			else {
				#	We have cached data with the old format
				if ($recache == TRUE) {
					die "Attempted to resave all genbank gene positions in the new format, and it still failed"
				}
				$recache = TRUE;
			}
	    }
	    until ($dataOK == TRUE);
		 
	foreach my $gene (@all_gene_positions) {
	    push(@geneIDlist,$gene->{ID});
	}
    }
    return @geneIDlist;
  } # private_list_stable_ids_for_genbank_database #

  method private_get_internal_genbank_fasta_file (Str $region){
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    my $fastafile = $data_dumps_dir.'fa/'.$region.'.fa';

    if (!-e $fastafile) { # if file doesn't already exist, fetch it
    	mkdir ($data_dumps_dir.'fa');
	    my $rawfile = $self->private_fetch_genbank_file($region);    
    	my $seqio_object = Bio::SeqIO->new(-file => $rawfile,-format => 'genbank');
    	my $seqObj = $seqio_object->next_seq();
    
    	my $fastafile = $data_dumps_dir.'fa/'.$region.'.fa';
    
    	my $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => '>'.$fastafile);
    	$seqout->write_seq($seqObj);
    }
	return $fastafile  	
  }
  
  method private_get_genomic_unmasked_sequence_of_genbank_gi (Genome_Sequence_Database_Parameters $parameters, Str $coord_sys_name, Str $region, Int $fiveprime, Int $threeprime, StrandType $strand) {

      my $fastafile = $self -> private_get_internal_genbank_fasta_file($region);

      my $startposition = $fiveprime;
      my $endposition = $threeprime;

      if ($strand eq 'negative') {
	  $startposition = $threeprime;
	  $endposition = $fiveprime;
      }
      my $sequence;

      my $number_of_Ns_start = 0;
      my $number_of_Ns_end = 0;

#    my $masked_sequence;
      my $difference = ($endposition - $startposition);
      $GU->user_info(3, "diff is $difference\n");
      if ($difference < 0) { 
	  die 'ill-defined genomic positions.';
      }
      else {		  
	  my $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => '<'.$fastafile);
	  my $slice = $seqin -> next_seq();
	  if ($startposition<1) {
	      if ($endposition<1) {
		  $number_of_Ns_start = 1+$difference;
		  $startposition = 1;
		  $endposition = 1;
	      } else {
		  $number_of_Ns_start = 1-$startposition;
		  $startposition = 1;
	      }
	  }
	  my $slice_length = $slice->length();
	  if ($endposition>$slice_length) {
	      $number_of_Ns_end = $endposition-$slice_length;
	      $endposition = $slice_length;
	  }
	  $sequence = $slice->subseq($startposition,$endposition);
	  $sequence = uc($sequence);
	  my $start_extension = $GU->repeat_character('N',$number_of_Ns_start);
	  my $end_extension = $GU->repeat_character('N',$number_of_Ns_end);
	  $sequence = $start_extension.$sequence.$end_extension;
	  if ($strand eq 'negative') { 
	      #make from strings in script
	      my $seqobj = Bio::Seq->new(-seq => $sequence);
	      my $rev = $seqobj->revcom();
	      $sequence = $rev->seq();
	  }
##    my $hardmasked_sequence = $slice->get_repeatmasked_seq(); # changes type to Bio::EnsEMBL::RepeatMaskedSlice
#    $masked_sequence = $slice->subseq($startposition,$endposition);
##    $masked_sequence = $hardmasked_sequence->seq();
		  
#    if ($strand eq 'negative') { 
#	# make from strings in script
#	  my $seqobj = Bio::Seq->new(-seq => $masked_sequence);
#	  my $rev = $seqobj->revcom();
#	  $masked_sequence = $rev->seq();
#      }
      }
#    my @result;
#    push (@result, $sequence);
#    push (@result, $masked_sequence);
		 
      return $sequence;
  } # private_get_genomic_unmasked_sequence_of_genbank_gi #

  method private_genbank_proteins_to_fasta_file (Genbank_Sequence_Database_Parameters $parameters, Str $filename, TranscriptChoice $transcript_choice) {
    $GU->user_info(1,"XXX making fasta file from ". $parameters->dbname."\n");
    my $outfile = Bio::SeqIO->new(-file =>  '>'.$filename, -format => 'fasta');

    my $dbname = $parameters->dbname;
    my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
    mkdir($data_dumps_dir.'inf');

	my @geneids = $self->list_stable_ids_for_a_given_genome_database($parameters);
	my $Ngenes = @geneids;
	my $curGene = 0;
	while ($curGene < $Ngenes) {
	    my @transcripts;
	    if (not(defined($geneids[$curGene]))) {
	    	$curGene++;
	    	next;
	    }
        my $geneinf = $self->private_get_genbank_geneinfo($geneids[$curGene]);
        my @proteins = @{$geneinf ->{PROTEINS}};
        if (scalar(@proteins) && ($proteins[0] ne 'unknown')) {
	        my $protein = $proteins[0];
	        my $cdsinf = $self->private_get_genbank_cdsinfo($protein);
	        my $sequence;
	        my $count = $Ngenes - $curGene;
	        if ($count > 10){$count = 10;};
	        my @cdsInfsToFetch;
	        my @cdsFilesToFetch;
	        
			my $fetchGene = $curGene;
			
			my $id = $geneinf ->{ID};
	        
	        while (defined($cdsinf) && !(defined($cdsinf->{TRANSLATION})) && ($count > 0)) {
			    for (my $i = 0;$i < scalar(@proteins);$i++) {
					$cdsinf = $self->private_get_genbank_cdsinfo($proteins[$i]);
					if (!(defined($cdsinf->{TRANSLATION}))){
						push (@cdsInfsToFetch,$cdsinf);
						push (@cdsFilesToFetch,$cdsinf -> {ID});
					}
		    	}
			    $fetchGene++;
			    $count--;
			    if ($count > 0) {
					$geneinf = $self->private_get_genbank_geneinfo($geneids[$fetchGene]);
					@proteins = @{$geneinf ->{PROTEINS}};
					if (scalar(@proteins)) {
					    $cdsinf = $self->private_get_genbank_cdsinfo($proteins[0]);
					}
					else {
					    undef($cdsinf);
					}
		    	}
	        }
	        if (scalar(@cdsInfsToFetch) > 0) {
	        	my $count = 0;
	        	# There may be a number like this
		      	my @rawfiles = $self->private_fetch_genbank_files(\@cdsFilesToFetch);
		      	foreach my $rawfile(@rawfiles) {
			      	my $seqio_object = Bio::SeqIO->new(-file => $rawfile,-format => 'genbank');
			    	while (my $seqObj = $seqio_object->next_seq()) {
				    	my $translation = $seqObj -> seq();
				    	my $cdsinf = shift(@cdsInfsToFetch);
				    	$cdsinf->{TRANSLATION} = $translation;
        				my $cdsfileName = $data_dumps_dir.'inf/'.$cdsinf -> {ID}.'.cdsInf';
				        nstore(\$cdsinf,$cdsfileName);
					    chmod (0664, $cdsfileName);
		    		}
		    	}
		    	$geneinf = $self->private_get_genbank_geneinfo($geneids[$curGene]);
				@proteins = @{$geneinf ->{PROTEINS}};
		        $protein = $proteins[0];
	        	$cdsinf = $self->private_get_genbank_cdsinfo($protein);
		   	}

			if ($cdsinf -> {TRANSLATION} eq '')
			{
				$GU->user_info( 3, "No sequence for ". $cdsinf -> {ID}." in ".$id."\n" );
			}
			else
			{
				$GU->user_info( 3, "sequence exists for ".$cdsinf -> {ID}." in ".$id."\n" );
				
				my $bioseq = Bio::Seq->new(-id => $id, -seq => $cdsinf -> {TRANSLATION});
				$outfile->write_seq($bioseq); # not using '-id => $transcript_to_translate->display_id' here any more, because we want to retain the original gene id for the purposes of ortholog mapping - what about in other cases though? add another parameter to use gene id/display id?
			}
      		
			if ($transcript_choice eq 'all') {
			    for (my $i = 1;$i < scalar(@proteins);$i++) {
			        $cdsinf = $self->private_get_genbank_cdsinfo($proteins[$i]);
					if ($cdsinf -> {TRANSLATION} eq '')
					{
						$GU->user_info( 3, "No sequence for ".$cdsinf -> {ID}." in ".$id."\n" );
					}
					else
					{
						$GU->user_info( 3, "sequence exists for ".$cdsinf -> {ID}." in ".$id."\n");
						my $bioseq = Bio::Seq->new(-id => $id."_".$i+1, -seq => $cdsinf -> {TRANSLATION});
						$outfile->write_seq($bioseq); # not using '-id => $transcript_to_translate->display_id' here any more, because we want to retain the original gene id for the purposes of ortholog mapping - what about in other cases though? add another parameter to use gene id/display id?
					}
			    }
			}
	        $curGene++;
        }
	    else {
			$GU->user_info( 3, "No sequence for ". $geneinf->{ID}."\n" );
			$curGene++;
	    }
	}
    $GU->user_info(3,$filename." file written\n");
    return 1;
  } # private_genbank_proteins_to_fasta_file #

  method private_genbank_genome_to_fasta_file (Genbank_Sequence_Database_Parameters $parameters) {
  #
  #   Fetches  
    my $dbname = $parameters->dbname;
    my $outfile = Bio::SeqIO->new(-file => '>'. $dbname. '.genome.fasta',
				  -format => 'fasta');

    foreach my $file (@{ $parameters->filenames}) {
    	my $fastafile = $self ->private_get_internal_genbank_fasta_file($file);
    	my $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => '<'.$fastafile);
        my $slice = $seqin -> next_seq();
        $outfile->write_seq($slice)
    }
    return 1;
  } # private_ensembl_genome_to_fasta_file #

  method private_get_all_chromosomes_genbank(Genbank_Sequence_Database_Parameters $parameters) {
      my @genomic_intervals;
      foreach my $file (@{ $parameters->filenames}) {
	  my $fastafile = $self->private_get_internal_genbank_fasta_file($file);
	  my $seqin = Bio::SeqIO->new( -format => 'Fasta', -file => '<'.$fastafile);
	  my $slice = $seqin -> next_seq();
	  my $slice_length = $slice->length;
	  $GU->user_info(3,"Region ".$file." has length ".$slice_length.".\n");
	  my $genomic_interval = Genomic_Interval->new(
	      'genome_db' => $parameters,	 
	      'region' => $file,
	      'five_prime_pos' => 1,
	      'three_prime_pos' => $slice_length,
	      'strand' => 'positive',
	      'type' => 'dna',
	      ); 
	  push(@genomic_intervals,$genomic_interval);
      }
      my $gis_maker = Genomic_Interval_Set_Maker->new();
      my $gi_set = $gis_maker->make_gi_set_from_gi_array(\@genomic_intervals);
      return $gi_set;
  } # private_get_all_chromosomes_genbank #

  method private_genbank_all_transcripts_to_fasta_file (Genbank_Sequence_Database_Parameters $parameters) { # 
#
# For genbank genomes we need to get the complete genome, which we cache locally (dataDumps/fa/*.fa)
    my %fastaFiles;
    foreach my $file (@{ $parameters->filenames}) {
    	my $fastafile = $self ->private_get_internal_genbank_fasta_file($file);
    	my $seqin= Bio::SeqIO->new( -format => 'Fasta', -file => '<'.$fastafile);
        my $slice = $seqin -> next_seq();
    	$fastaFiles{$file} = $slice;
    }

# And then we get the list of the names for all of the genes 
    my $dbname = $parameters->dbname;
    my $outfile = Bio::SeqIO->new(-file => '>' . $dbname . '.all_transcripts.fasta',
				  -format => 'fasta');
    
    my @geneids = $self->list_stable_ids_for_a_given_genome_database($parameters);
    
    my $Ngenes = @geneids;
    my $curGene = 0;
    
#   And for each of the genes we need the start and end positions so that we can extract the DNA sequence
#   for each gene
    while ($curGene < $Ngenes) {
      
      my $geneinf = $self->private_get_genbank_geneinfo($geneids[$curGene]);
      $curGene++;
      
      my $startposition = $geneinf -> {START};
      my $endposition = $geneinf -> {END};
      my $region = $geneinf -> {REGION};
      my $strand = $geneinf -> {STRAND};
      my $id = $geneinf -> {ID};
      
      #   The region information tells use which chromosome to use.   Particularly important 
      #   when we have genbank data for an organisim with multiple chromosomes.
      my $slice = $fastaFiles{$region};
      my $seqobj;

      if ($strand eq 'positive') {
        my $sequence = $slice->subseq($startposition,$endposition);
        $seqobj = Bio::Seq->new(-seq => $sequence,-id  => $id);
      }
      else {
        my $sequence = $slice->subseq($startposition,$endposition);
        my $posseq = Bio::Seq->new(-seq => $sequence,-id  => $id);
        $seqobj = $posseq->revcom();
        }
      $outfile->write_seq($seqobj);
      }
      
    return 1; # file written
  } # private_genbank_all_transcripts_to_fasta_file #


  method private_get_genbank_geneinfo(Str $geneid){
      my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
      my $genfile = $data_dumps_dir.'gene/'.$geneid.'.geneInf';
      my $geneInf = ${retrieve ($genfile)};
      return $geneInf;  	
  } # private_get_genbank_geneinfo

  method private_get_genbank_cdsinfo(Str $cdsid){
      my $data_dumps_dir = $ENV{'data_dumps'}; # download data dumps to here
      my $genfile = $data_dumps_dir.'inf/'.$cdsid.'.cdsInf';
      my $cdsInf = ${retrieve ($genfile)};
      return $cdsInf;  	
  } # private_get_genbank_cdsinfo #

  method private_get_genbank_gene_start_position (Str $geneid) {
  	  my $geneinf = $self->private_get_genbank_geneinfo($geneid);
      my $position;
      if ($geneinf->{STRAND} eq 'positive') {
		  $position = $geneinf->{START};
      }
      else{
		  $position = $geneinf->{END};
      }
      return $position;
  } # private_get_genbank_gene_start_position #
  
  
  method private_get_genbank_gene_parameters (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
      my $geneinf;
      if ($geneid eq 'not found'){
      	$geneinf ->{REGION} = '';
      	$geneinf ->{STRAND} = 'positive';
      	$geneinf ->{START} = 0;
      	$geneinf ->{END} = 0;
      }
      else {
        $geneinf = $self->private_get_genbank_geneinfo($geneid);
      }
      my $params = Genomic_Interval_Parameters -> new(
        geneid => $geneid,
        database => $parameters,
        genomic_region_type => 'gene_space',
        region => $geneinf->{REGION},
        strand => $geneinf->{STRAND},
        five_prime_pos => $geneinf->{START},
        three_prime_pos => $geneinf->{END}
        );
      return $params;
  } # private_get_genbank_gene_parameters #

  method private_get_genbank_go_terms( Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
    my $geneinf;
    my %go_terms;
    my $data_dumps_dir = $ENV{'data_dumps'};          # download data dumps to here

    if ( $geneid ne 'not found' ) {
      $geneinf = $self->private_get_genbank_geneinfo($geneid);
      my @proteins = @{ $geneinf->{PROTEINS} };
      if ( scalar(@proteins) && ( $proteins[0] ne 'unknown' ) ) {
	foreach my $protein ( @proteins ) {
	  my $cdsinf = $self->private_get_genbank_cdsinfo($protein);
	  if ( !defined( $cdsinf->{GO_PROCESS} ) ) {
	    my @go_component;
	    my @go_process;
	    my $rawfile = $self->private_fetch_genbank_file($protein);

	    my $seqio_object = Bio::SeqIO->new( -file => $rawfile, -format => 'genbank' );
	    while ( my $seqObj = $seqio_object->next_seq() ) {
	      foreach my $featObj ( $seqObj->get_SeqFeatures() ) {
		if ( $featObj->primary_tag() eq 'CDS' ) {
		  if ( $featObj->has_tag('note') ) {
		    my $note = ( $featObj->get_tag_values('note') )[0];
		    my @sections = split( ';', $note );
		    foreach my $section (@sections) {
		      $section =~ s/^\s+//;    #Strip leading spaces
		      if ( substr( $section, 0, 10 ) eq 'GO_process' ) {
			my $process = substr( $section, 11 );
			$process =~ s/^\s+//;    #Strip leading spaces
			push( @go_process, $process );
		      }
		      elsif ( substr( $section, 0, 12 ) eq 'GO_component' ) {
			my $component = substr( $section, 13 );
			$component =~ s/^\s+//;    #Strip leading spaces
			push( @go_component, $component );
		      }
		    }
		  }
		}
	      }
	    }
            $cdsinf->{GO_PROCESS}   = \@go_process;
	    $cdsinf->{GO_COMPONENT} = \@go_component;

            my $cdsfile = $data_dumps_dir.'inf/'.$protein.'.cdsInf';
            nstore( \$cdsinf, $cdsfile );
            chmod( 0664, $cdsfile );
	  }
          my $goData = Genome_CDS_GOterm_Data -> new(
            go_process => $cdsinf->{GO_PROCESS},
	    go_component => $cdsinf->{GO_COMPONENT});
	  $go_terms{$protein} = $goData;
	}
      }
    }
    return \%go_terms;
  }    # private_get_genbank_gene_parameters #

  method private_get_genbank_proteins( Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
    my $geneinf;
    my @proteins;
    my $data_dumps_dir = $ENV{'data_dumps'};          # download data dumps to here

    if ( $geneid ne 'not found' ) {
      $geneinf = $self->private_get_genbank_geneinfo($geneid);
      @proteins = @{$geneinf->{PROTEINS}};
    }
    return @proteins;
  }    # private_get_genbank_gene_parameters #
  

#############################################
##### FASTA-DB-specific private methods #####
#############################################

  method private_list_stable_ids_for_fasta_database (FASTA_Sequence_Database_Parameters $parameters) {
      my @result;

      my $inputstream = Bio::SeqIO->new( -file => $parameters->filename, -format => 'Fasta' );
      while ( my $sequence = $inputstream->next_seq() ) {
	  push(@result,$sequence->id);
      }
      return @result;
  } # private_list_stable_ids_for_fasta_database #

  method private_retrieve_all_sequences_from_fasta_DB (FASTA_Sequence_Database_Parameters $parameters) {
      my $seqio_object = Bio::SeqIO->new(-file => $parameters->filename);
      my @sequence_array;
      my $rec;
      while (my $seq = $seqio_object->next_seq) {
	  $rec = {
	      NAME   => $seq->display_id,
	      SEQUENCE => $seq->seq
	  };
	  push (@sequence_array, $rec);
      }
      return @sequence_array;
  } # private_retrieve_all_sequences_from_fasta_DB #

  method private_get_genomic_unmasked_sequence_of_fasta_gi (FASTA_Sequence_Database_Parameters $parameters, Str $coord_sys_name, Str $region, Int $fiveprime, Int $threeprime, StrandType $strand) {

      my $inputstream = Bio::SeqIO->new( -file => $parameters->filename, -format => 'FASTA' );
      my $sequence = $inputstream->next_seq();
      while ($sequence->display_id() ne $region) {
	  if (!($sequence = $inputstream->next_seq())) {
	      die 'sequence ('.$region.') does not exist in file('.$parameters->filename.').';
	  }
      }
      my $string = $sequence->seq();
      my @bothcoordinates = ($fiveprime,$threeprime);
      my $min = $GU->minimum(\@bothcoordinates);
      my $max = $GU->maximum(\@bothcoordinates);
      if ($min<1) {
	  die 'invalid coordinate';
      }
      if ($max>length($string)) {
	  die 'invalid coordinate';
      }
      my $substring = substr($string,$min,($max-$min)+1);
      if ($strand eq 'negative') { 
	  my $bioseq = Bio::Seq->new( -seq => $substring );
	  $substring = $bioseq->revcom()->seq();
      }  
      return $substring;
  } # private_get_genomic_unmasked_sequence_of_fasta_gi #

} # Genome_DB_Utilities #
