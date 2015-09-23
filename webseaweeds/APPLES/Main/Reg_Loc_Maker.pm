### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Reg_Loc_Maker class ###
# methods for creating sets of regulated loci
use MooseX::Declare;

class Reg_Loc_Maker {
    use Reg_Loc;
    use APPLES_Datatypes qw(APPLESSpeciesName LocusType Boolean);
    use constant {FALSE => 0,
		  TRUE	=> 1};
    use Genome_DB_Utilities;
    use Data::Dumper;
    use General_Utilities;
    use Parameters;
    
    my $GDBU = Genome_DB_Utilities->new();
    my $GU = General_Utilities->new();

    method make_all_reg_locs_for_a_genome(Genome_Sequence_Database_Parameters $parameters) {
	my $genome_seq_handler = Genome_DB_Utilities->new();
	
	my @geneIDlist = $genome_seq_handler->list_stable_ids_for_a_given_genome_database($parameters);
	my @reg_locs = $self->make_reg_locs_from_list_through_job_handler($parameters,\@geneIDlist);
	return @reg_locs;
    } # make_all_reg_locs_for_a_genome #
  
    method make_reg_locs_from_list (Genome_Sequence_Database_Parameters $parameters, ArrayRef $geneIDlist) {
	my @reg_locs;
	my $progress_count = 0;
	my $list_length = @{$geneIDlist};
	foreach my $geneID (@{$geneIDlist}) {
	    $progress_count++;
	    if (($progress_count % 10) == 0) {
		$GU->user_info(2,$progress_count."/".$list_length." reg_locs processed.\n");
	    }
	    if (($progress_count % 1000) == 0) {
		$GDBU = Genome_DB_Utilities->new();
		$GU->user_info(3,"re-newed GDBU.\n");
	    }
	    my $reg_loc = $self->make_reg_loc_from_gene_id($parameters,$geneID);
	    push(@reg_locs,$reg_loc);
	}
	return @reg_locs;
    } # make_reg_locs_from_list #
    
    method make_reg_locs_from_list_through_job_handler (Genome_Sequence_Database_Parameters $parameters, ArrayRef $geneIDlist) {
	my $object = Reg_Loc_Maker->new();
	my $function = 'make_reg_locs_from_list';
	my @all_parameters = ($parameters,$geneIDlist);
	my $high_memory = TRUE;
	my @result = $GU->standard_call_to_job_handler($object,$function,\@all_parameters,$high_memory,FALSE);
	return @result;
    } # make_reg_locs_from_list_through_job_handler #
  
    method make_reg_locs_from_chip_locations(Str $filename, Genome_Sequence_Database_Parameters $EDP) {
    # Sascha: may be more natural to produce ReMo objects rather than Reg_Locs for ChIP-fragments - nobody seems to call this function anyway
	my @reg_locs;
	# read file in line at a time
	open (FILE, $filename);
	while (<FILE>) {	
	    chomp;	
	    my ($encode, $position, $region) = split("\t");
	    my $strand = 'positive';
	    my $coord_sys = 'CHROMOSOME';
	    
	    my $reg_loc = Reg_Loc->new(genome_db=>$EDP,
				       coord_sys_name=>$coord_sys,
				       region=>$region,
				       position=>$position,
				       strand=>$strand,					
				       locus_type=>'chip_fragment',
				       gene_ID=>'',
				       transcript_ID=>'');
	    push (@reg_locs, $reg_loc);									  
	}
	close (FILE);
	return @reg_locs;
    } # make_reg_locs_from_chip_locations #
    
    method make_reg_loc_from_gene_id (Genome_Sequence_Database_Parameters $parameters, Str $geneid) {
    	
    my $coord_sys_name;	
    my $seq_region_name;
    my $position;
    my $strand;
    
    while (!$main::global_perseverance->stop_trying){
	    eval {	
			$coord_sys_name = $GDBU->get_coord_system_name_from_gene_id($parameters, $geneid);
			$seq_region_name = $GDBU->get_seq_region_name_from_gene_id($parameters, $geneid);
			$position = $GDBU->get_gene_start_position($parameters, $geneid);
			$strand = $GDBU->get_strand_of_gene($parameters, $geneid);			
	    };
	    my $error = $@;
	    $main::global_perseverance->decide_on_rerun('Ensembl',TRUE,$error);
	    if ($error) {
	    	$GDBU = Genome_DB_Utilities->new();
    	}		
    }
    $main::global_perseverance->stop_trying(FALSE);        
    my $reg_loc = Reg_Loc->new(genome_db => $parameters,
			       coord_sys_name => $coord_sys_name,
			       region => $seq_region_name,
			       position => $position,
			       strand => $strand,					
			       locus_type => 'gene',
			       gene_ID => $geneid,
			       transcript_ID => '');
    return $reg_loc;
    } # make_reg_loc_from_gene_id #

} # Reg_Loc_Maker #
