### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### BiFa_Raw_Result Class ###

use MooseX::Declare;

class BiFa_Raw_Result {
	use Data::Dumper;
	use Parameters; 
	use General_Utilities;

	has 'raw_result' => (is => 'ro', isa => 'ArrayRef', required => 1);

	my $GU = General_Utilities->new();

	method process_result {
		my @ordered_sites_one_result;
		my @sites_all_results;
		my $gi_counter = 0;
		
		my @all_raw_scores = @{$self->raw_result};
		
		for my $i (0..@all_raw_scores-1) {
			
			if (exists($all_raw_scores[$i])) {
				my $s = $all_raw_scores[$i];
				my @one_result_raw_scores = @{$s};
				$GU->user_info( 2, "going to get ordered sites...\n" );
				@ordered_sites_one_result = $self->all_sites_ordered_by_score(\@one_result_raw_scores, $gi_counter);
				$GU->user_info( 3, Dumper(@ordered_sites_one_result) );
				push (@sites_all_results, @ordered_sites_one_result);
			}
			else {
				$GU->user_info( 3, "no scores for this sequence!\n" );
			}
			
			$gi_counter++;

		}
			$GU->user_info( 3, Dumper (@sites_all_results) ); 
			
		my @ordered_sites_all_results = $self->ordered_sites_all_results(\@sites_all_results);
			
			$GU->user_info( 3, Dumper (@ordered_sites_all_results) );
	} # process_result #
	
	method all_sites_ordered_by_score (ArrayRef $one_result_ref, Int $gi_counter) {
	
		my @raw_result = @{$one_result_ref};
		
		my @positive_strand_scores;
		my @negative_strand_scores;
		my $strand = 'NONE';
		my $rec;
		my @recs1;
		my @recs2;
		my @recs;
		
		for (my $i = 0; $i < scalar(@raw_result); $i += 2) {
			push(@positive_strand_scores, $raw_result[$i]);
			push(@negative_strand_scores, $raw_result[$i+1]);
		}
		for ( my $i = 0; $i < scalar(@positive_strand_scores); $i++ ) {
			$strand = 'positive';
			$rec = {
				INTERVALINDEX => $gi_counter,
				SITE   => $i,
				STRAND => $strand,
				PVALUE => $positive_strand_scores[$i]
			};
			push (@recs1, $rec);
		}

		for ( my $i = 0; $i <  scalar(@negative_strand_scores); $i++ ) {
			$strand = 'negative';
			$rec = {	
				INTERVALINDEX => $gi_counter,
				SITE   => $i,
				STRAND => $strand,
				PVALUE => $negative_strand_scores[$i]
			};
			push (@recs2, $rec);
		}
		
		@recs = (@recs1, @recs2);
		
		my @sorted_recs = (sort { $a->{PVALUE} <=> $b->{PVALUE} } @recs); # smallest->biggest

		return @sorted_recs;
	} # all_sites_ordered_by_score #
	
	method ordered_sites_all_results (ArrayRef $sites_all_results_ref) {

		my @sites_all_results = @{$sites_all_results_ref};
		
		my @sorted_result = (sort { $a->{PVALUE} <=> $b->{PVALUE} } @sites_all_results);
		
		return @sorted_result;
		
	} # ordered_sites_all_results #

#### these methods below are not used in code yet - might not be needed ####
	method significant_sites (ArrayRef $ordered_sites, BiFa_Parameters_Core $parameters) {

		my $cutoff = $parameters->significance_threshold;
		my @validscores;
		foreach my $rec (@{$ordered_sites}) {
			$GU->user_info( 3, Dumper ($rec) );
			if ($rec->{PVALUE} > $cutoff) {
				$GU->user_info( 3, $rec->{PVALUE}."\n" );
				$GU->user_info( 3, "significant hit!\n" );
				push (@validscores, $rec);
			}
		}
		$GU->user_info( 3, "valid scores array:\n" );
		$GU->iser_info( 3, Dumper (@validscores) );	
		return @validscores;
	} # significant_sites #

	method find_best_site (ArrayRef $raw_result_ref) {

		my @raw_result = @{$raw_result_ref};

		my @positive_strand_scores;
		my @negative_strand_scores;
		my $best_score = 0;
		my $strand = 'NONE';
		my $rec;
		for (my $i = 0; $i < scalar(@raw_result); $i += 2) {
			push(@positive_strand_scores, $raw_result[$i]);
			push(@negative_strand_scores, $raw_result[$i+1]);
		}
		
		for ( my $i = 0; $i < scalar(@positive_strand_scores); $i++ ) {
			$strand = 'positive';
			if ($positive_strand_scores[$i] > $best_score) {
				$rec = {
					SITE   => $i,
					STRAND => $strand,
					SCORE => $positive_strand_scores[$i]
				};
				$best_score = $positive_strand_scores[$i];
			}
		}

		for ( my $i = 0; $i <  scalar(@negative_strand_scores); $i++ ) {
			$strand = 'negative';
			if ($negative_strand_scores[$i] > $best_score) {
				$rec = {
					SITE   => $i,
					STRAND => $strand,
					SCORE => $negative_strand_scores[$i]
				};
				$best_score = $negative_strand_scores[$i];
			}
		}

		$GU->user_info (3, $best_score."\n");
		my $best_site = $rec; 
		return $best_site;
	} # find_best_site #
	
} # BiFa_Raw_Result #
