### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### BiSi_Set Class: set of binding sites ###
use MooseX::Declare;
use Serializable;


class BiSi_Set extends Serializable {
	use Parameters;
	use BiSi;
	use Data::Dumper;
	use APPLES_Datatypes qw (Boolean PositiveInt);
	use constant {FALSE => 0,
		      TRUE	=> 1};	

	has 'bisi_set' =>	(is => 'ro', isa => 'ArrayRef[BiSi]', required => 1);
	
	my $GU = General_Utilities->new();
	
	# Returns BiSi Set consisting of a single Bisi
	method expand_bisi_set( PositiveInt $five_prime_flank_length, PositiveInt $three_prime_flank_length, Boolean $mask_original_bisi){
		
		# For each BiSi in the BiSi_Set, create a new genomic interval that has a modified
		# sequence based on expanded flanking regions
		
		# expand binding site by flanking bases
		my @genomic_interval_set;
		my $gi_maker = Genomic_Interval_Maker->new();
		my @genomic_intervals_to_merge;
		foreach my $BISI ( @{$self->bisi_set()} ){
			foreach my $genomic_interval ( @{ $BISI->genomic_interval_set()->genomic_interval_set() } ){
				my $modified_gi = Genomic_Interval_Maker->copy_genomic_interval($genomic_interval);
				$modified_gi->expand_interval($five_prime_flank_length,$three_prime_flank_length);
				push(@genomic_intervals_to_merge, $modified_gi);
			}
		}
		# mask binding sites
		my $gi_set_maker = Genomic_Interval_Set_Maker->new();
		my $gi_set = $gi_set_maker->make_gi_set_from_gi_array(\@genomic_intervals_to_merge);
		my $overlap_free_gi_set = $gi_set->merge(TRUE);
		if ($mask_original_bisi) {
			my $number = $overlap_free_gi_set->return_number_of_intervals();
			for (my $index = 0; $index < $number; $index++) {
				# @{ $BISI->genomic_interval_set()->genomic_interval_set() }
				foreach my $bisi_part ( @{${$self->bisi_set}[$index]->genomic_interval_set()->genomic_interval_set()} ){
					${$overlap_free_gi_set->genomic_interval_set}[$index]->mask_genomic_interval( $bisi_part );
				}
			}
		}
		return $overlap_free_gi_set;
	} # expand_bisi_set #
	
	method get_top_n_pvalues(Int $n){
		
		# If n pvalues requested is greater than the total number of  bisi's in set then crash
		
		if ( $n > scalar(@{$self->bisi_set}) ) {
				
			die "Number of p-values requested greater than the number available (# of BiSi's)\n";
		}
		
		my @top_n_pvalues;
		
		if ( $self->private_has_pvalues() ) {
			
			# Sort bisi set by pvalue
			my @sorted_bisis = @{$self->private_sort_bisi_set_by_pvalue()};

			my $index = 0;
			
			for (my $i = 0; $i < $n; $i++) {
				
				my $bisi = $sorted_bisis[$i];
				
				# Identify the index of this BiSi in the BiSi_Set objects array of BiSi objects
				my $index = 0;
				++$index until ${$self->bisi_set}[$index] == $bisi;

				my $pvalue = {
					PVALUE => $bisi->single_site_pvalue,
					BISI_SET_INDEX => $index
				};
				
				push(@top_n_pvalues, $pvalue);
				
			}
			
			return \@top_n_pvalues;
		}
		else {		
			die "Not all of the BiSi's in this BiSi_Set contain pvalue scores\n";
		}
	} # get_top_n_pvalues # 

	method get_number_of_bisis() {
	    my $result = @{$self->bisi_set};
	    return $result;
	} # get_number_of_bisis #
	
	method private_sort_bisi_set_by_pvalue(){
		
		my @sorted_bisis = sort{ ($a->single_site_pvalue) <=> ($b->single_site_pvalue) } @{$self->bisi_set};
		
		return \@sorted_bisis;
	} # private_sort_bisi_set_by_pvalue # 
	
	method private_has_pvalues(){
	
		foreach my $bisi (@{$self->bisi_set}){
			if(!defined ($bisi->single_site_pvalue)){
				$GU->user_info(1,"undefined bisi pvalue\n");
				return FALSE;
			}
		}
		
		return TRUE;
	} # private_has_pvalues #
	
#	override to_hash() {
#        
#        my @bisi_array;
#        
#        foreach my $bisi (@{$self->bisi_set}) {
#        	push(@bisi_array, $bisi->to_array());
#        }
#        
#        my $hash = {'bisi_set' => {'bisi' => [@bisi_array]}};
#        
#        return $hash;     
#    }
	
} # BiSi_Set #
