=pod

=head1 Class ReMo_Maker

ReMo_Maker - methods to create ReMo(s).

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<nmake_remos_from_reglocs>

Returns a new My::Module object.

=item C<make_remo_from_core_promoter>

Returns a new My::Module object.

=item C<make_remos_from_orthologous_core_promoters>

Returns a new My::Module object.

=item C<make_remos_from_chip_seq_fragment_set>

Returns a new My::Module object.

=item C<make_remos_from_known_bisi>

Returns a new My::Module object.

=back

=head1 LICENSE

This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package and distributed 
under Academic Non-Commercial Use Licence.

=head1 COPYRIGHT

(c) Copyright University of Warwick 2009-2010

=head1 AUTHOR

=cut

use MooseX::Declare;

class ReMo_Maker {
  use Parameters;
  use ReMo;
  use BiSi;
  use Genome_DB_Utilities;
  use Data::Dumper;
  use Job_Handler;
  use APPLES_Datatypes qw (APPLESSpeciesName Boolean PositiveInt);
  use constant {FALSE => 0,
		TRUE	=> 1};	
  use ChIP_Seq_Fragment_Set;

=pod
    
=head2 make_remos_from_reglocs( ReMo_Constructor_Parameters $parameters, Reg_Loc $reg_loc )

This method makes a ReMo from the specified RegLoc using specified parameters.

Return: Returns a ReMo object.

=cut
 
  method make_remos_from_reglocs(ReMo_Constructor_Parameters $parameters, Reg_Loc $reg_loc) {
    
    my $neighbour_search_direction;
  	my $five_prime_pos;
    my $three_prime_pos;
    my $gdbu = Genome_DB_Utilities->new();
 
    if ($reg_loc->strand eq 'positive') {
       $neighbour_search_direction = 'towards_five_prime';
    } else {
       $neighbour_search_direction = 'towards_three_prime';
    }
    
    my $distance_to_neighbouring_gene = $gdbu->get_distance_to_neighbouring_gene($reg_loc, 
                                              $parameters->ignore_pseudogenes,
                                              $parameters->ignore_RNAgenes,
                                              $neighbour_search_direction);                                
                                 
  	if ($reg_loc->strand eq 'positive' ) {
        $five_prime_pos = $reg_loc->position - $distance_to_neighbouring_gene;
        $three_prime_pos = $reg_loc->position - 1;
    } else {
        $five_prime_pos = $reg_loc->position + $distance_to_neighbouring_gene; 
        $three_prime_pos = $reg_loc->position + 1;
    }
    
    my $remo = ReMo->new('genome_db' => $reg_loc->genome_db,
             'coord_sys_name' => $reg_loc->coord_sys_name,
             'region' => $reg_loc->region,
             'five_prime_pos' => $five_prime_pos,
             'three_prime_pos' => $three_prime_pos,
             'strand' => $reg_loc->strand,
             'label' => $reg_loc->gene_ID
              );
    
    
    $remo->get_working_sequence();
    
    return $remo;
  }

  method make_remo_from_core_promoter(ReMo_Core_Promoter_Constructor_Parameters $parameters,
				      Reg_Loc $reg_loc,
                      Boolean $obtain_sequences) {
    
    my $gdbu = Genome_DB_Utilities->new();
    my $GU = General_Utilities->new();
    
    my $length = $parameters->length;
    my $minimum_length = $parameters->minimum_length;
    my $distance_to_neighbouring_gene;
    my $five_prime_pos;
    my $three_prime_pos;
	
    if ($length<=0) {
	die 'core promoter must have length greater than zero.';
    }
    if ($parameters->stop_at_neighbouring_gene) {
	   my $neighbour_search_direction;
	   
	   if ($reg_loc->strand eq 'positive') {
	       $neighbour_search_direction = 'towards_five_prime';
	   } else {
	       $neighbour_search_direction = 'towards_three_prime';
	   }
       
       my $neighbour_gene_exists = $gdbu->check_existence_of_neighbouring_gene($reg_loc, 
										$parameters->ignore_pseudogenes,
										$parameters->ignore_RNAgenes,
										$neighbour_search_direction);
	if ($neighbour_gene_exists) {
	    $distance_to_neighbouring_gene = $gdbu->get_distance_to_neighbouring_gene($reg_loc, 
										      $parameters->ignore_pseudogenes,
										      $parameters->ignore_RNAgenes,
										      $neighbour_search_direction);
	    if ($distance_to_neighbouring_gene < $length) {
		$length = $distance_to_neighbouring_gene;
	    }
	}
    }

    if ($length<$minimum_length) {
      $length = $minimum_length;
    }
    			
    if ($reg_loc->strand eq 'positive' ) {
	$five_prime_pos = $reg_loc->position - $length;
	$three_prime_pos = $reg_loc->position - 1;
    }
    else {
	$five_prime_pos = $reg_loc->position + $length; 
	$three_prime_pos = $reg_loc->position + 1;
    }
    
    my $remo = ReMo->new('genome_db' => $reg_loc->genome_db,
			 'coord_sys_name' => $reg_loc->coord_sys_name,
			 'region' => $reg_loc->region,
			 'five_prime_pos' => $five_prime_pos,
			 'three_prime_pos' => $three_prime_pos,
			 'strand' => $reg_loc->strand,
			 'label' => $reg_loc->gene_ID
			  );
    $GU->user_info(3,"DEBUG1\n");
    if ($obtain_sequences) {
	$remo->get_working_sequence();
    }
    $GU->user_info(3,"DEBUG2\n");
    return $remo;
  } # make_remo_from_core_promoter #
  
  method make_remos_from_orthologous_core_promoters(ReMo_Core_Promoter_Constructor_Parameters $parameters,  Genome_Sequence_Database_Parameters $genome_db, Ref $reg_loc_ref, Orthology_Method_And_Species_Restriction $orthology_method_and_species_restriction, Str $source_dbname) { # APPLESSpeciesName $source_alias
    # check $source_dbname is in the list of allowed_source_dbnames
    my @allowed_source_dbnames = @{$orthology_method_and_species_restriction->allowed_source_dbnames};
    my $orthology_parameters = $orthology_method_and_species_restriction->generic_orthology_finding_parameters;
    my $exists = 'false';
    foreach my $species (@allowed_source_dbnames) {
      if ($source_dbname eq $species) {
	$exists = 'true';
      }
    }
    if ($exists eq 'false') {
      die 'this orthology finding method is not valid for this source species $source_dbname (check the orthology_method_and_species_restriction parameters)';
    }
    my $gdbu = Genome_DB_Utilities->new();
    
    my @orthologous_gene_ids = $gdbu->get_orthologous_gene_IDs( ${$reg_loc_ref}->gene_ID, ${$reg_loc_ref}->genome_db, $genome_db, $orthology_parameters); 
    
    my @remos_one_species;
    # make into Reg_Locs first
    my $reg_loc_maker = Reg_Loc_Maker->new();
    my @orthologous_reg_locs = $reg_loc_maker->make_reg_locs_from_list( $genome_db,  \@orthologous_gene_ids);
    
    foreach my $orthologous_reg_loc (@orthologous_reg_locs) {
      # for each orthologous_reg_log, create a remo object
      
      my $remo_one_species = $self->make_remo_from_core_promoter ( $parameters, $orthologous_reg_loc,TRUE);
      
      push (@remos_one_species, $remo_one_species);
    }
    
    # to return a set (array) of ReMos 
    return @remos_one_species;	   
  } # make_remos_from_orthologous_core_promoters #		                   				   
  
  method make_remos_from_chip_seq_fragment_set(ChIP_Seq_Fragment_Set $fragments) {
    
    die 'method not implemented yet!';
    
    #my $GI = $ChIP_seq_fragment_set->genomic_interval;
    #my $ReMo_from_CSFS = ReMo->new(genomic_interval => $GI);
    
  } # make_remos_from_chip_seq_fragment_set #
  
  method make_remos_from_known_bisi(BiSi $bisi, Int $number_of_flanking_bases_upstream, Int $number_of_flanking_bases_downstream) {
    
    die 'method not implemented yet!';
    
    my $gi_of_bisi = $bisi->Genomic_Interval; # get genomic interval of BiSi   
  } # make_remos_from_known_bisi #	
  
} # ReMo_Maker #
  
	
