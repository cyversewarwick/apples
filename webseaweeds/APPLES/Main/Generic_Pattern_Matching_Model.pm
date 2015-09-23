### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Generic_Pattern_Matching_Model Class ###

use MooseX::Declare;

class Generic_Pattern_Matching_Model {
        use Generic_Sequence_Pattern;
	use Genomic_Interval;
	use Genomic_Interval_Set;
	use General_Utilities;
	use Data::Dumper;
	use ReRe_Set;
	use APPLES_Datatypes qw (PositiveInt Boolean);
	use constant {FALSE => 0,
		      TRUE	=> 1};	

	has 'pattern' => (is => 'ro', isa => 'Generic_Sequence_Pattern', required => 1);
	
	my $GU = General_Utilities->new();

	method is_pattern_in_sequence(Genomic_Interval_Set $gi_set) {
	    # returns a Boolean indicating presence or absence of pattern in $gi_set

	    die 'method not implemented yet';
	} # is_pattern_in_sequence #

	method is_pattern_in_sequence_rere_set(ReRe_Set $rere_set) {
	    my @result;
	    my $progress_count = 0;
	    my $number_of_reres = @{$rere_set->rere_set_members};
	    foreach my $rere (@{$rere_set->rere_set_members}) {
		$progress_count++;
		$GU->user_info(2,"processing ReRe ".$progress_count." out of ".$number_of_reres.".\n");
		my $gi_set = $rere->make_gi_set_from_rere(TRUE);
		my $p_value = $self->is_pattern_in_sequence($gi_set);
		push(@result,$p_value);
	    }
	    return @result;
	} # is_pattern_in_sequence_rere_set #
	
	method find_top_n_matches_in_sequence(Genomic_Interval_Set $gi_set, Int $n) {
	    # returns a BiSi_Set for top non-overlapping sites

	    die 'method not implemented yet';			
	} # find_top_n_matches_in_sequence #
	
	method number_of_possible_sites(Genomic_Interval_Set $gi_set) {
	    # returns Int indicating at how many positions pattern could be matched within $gi_set
 
	    die 'method not implemented yet';
	} # number_of_possible_sites #

	method binomial_overrepresentation_test(Genomic_Interval_Set $gi_set, PositiveInt $max_number_of_sites) {
          # uses binomial statistic to test if pattern is overrepresented looking at up to
	  # $max_number_of_sites best matches in $gi_set
	    my $top_bisis = $self->find_top_n_matches_in_sequence($gi_set,$max_number_of_sites);

	    my @bisi_bifa_assignments = ();
	    if ($self->pattern->isa('Generic_Weight_Matrix') ){
                push(@bisi_bifa_assignments,  $self->pattern->wm_identifier() );
	    }

	    my $number_of_sites = $top_bisis->get_number_of_bisis();
	    my $p_values = $top_bisis->get_top_n_pvalues($number_of_sites);
	    #$GU->user_info(1,Dumper($p_values));
	    my $number_of_sites_tested = $self->number_of_possible_sites($gi_set);
	    my $stats = Statistics->new();
	    $GU->user_info(3,"number of possible binding sites in binomial overrepresentation test: ".$number_of_sites_tested."\n");
	    my $result = $stats->best_binomial_pvalue($p_values, $number_of_sites_tested);
	    my $index = $$result{INDEX};
	    my @best_bisis;
	    for my $i ( 0 .. $index) {
		push(@best_bisis, ${$top_bisis->bisi_set}[$i]);
	    }

	    foreach my $best_bisi (@best_bisis) {
                 $best_bisi->bifa_assignment(\@bisi_bifa_assignments);
            }

	    $$result{N_BEST_BISIS} = \@best_bisis;
	    return $result;
	} # binomial_overrepresentation_test #
	
	method hypergeometric_overrepresentation_test(Boolean $test_for_significantly_large_overlap, ReRe_Set $rere_set, ArrayRef[ArrayRef[Str]] $gene_ID_clusters, Str $threshold, Boolean $use_job_handler) {
		
		# Use when there is one pattern model; returns one p-value.
		
		my $die_on_exception = TRUE;
		return $self -> private_hypergeometric_overrepresentation_test($test_for_significantly_large_overlap, $rere_set, $gene_ID_clusters, $threshold, $use_job_handler,$die_on_exception);
	
	} # hypergeometric_overrepresentation_test #
		
	method hypergeometric_overrepresentation_test_in_parallel(ArrayRef[Generic_Pattern_Matching_Model] $pattern_models,Boolean $test_for_significantly_large_overlap, ReRe_Set $rere_set, ArrayRef[ArrayRef[Str]] $gene_ID_clusters, Str $threshold, Str $progress_file = "") {
	    
	    # Use when there is an arrayref of pattern models. Returns an array of p-values, one for each pattern model.
	    # When this method is called, self just functions as a shell (there is no need for the attributes). 
	    	    
	    my @all_p_values;
	    my $use_job_handler = TRUE;
	    my $die_on_exception = FALSE; # die_on_exception is set to false, allowing continuation of the pattern_model loop after a job (is_pattern_in_rere_set) has been initiated.
	    my $aggregate_exception = Aggregate_Exception->new();
	    my $stats_required = FALSE;
	    my $number_of_patterns = @{$pattern_models};
	    my $pattern_count = 1;
	    foreach my $pattern_model (@{$pattern_models}) {
		$GU->user_info(2,"Doing hypergeometric test for matrix ".$pattern_model->pattern->wm_identifier." (".$pattern_count." out of ".$number_of_patterns.")\n"); 
		
		#last param defaults to empty string, but can be the name of a file used to monitor progress.
		if ($progress_file ne "") {
			#write proportion complete to file
			open(FP, ">$progress_file");
			my $prog = $pattern_count/$number_of_patterns;
			print FP "$prog\n";
			close FP;
		}
		
		
		$pattern_count++;
		my @p_values;
			eval {			    
			    @p_values = $pattern_model->private_hypergeometric_overrepresentation_test($test_for_significantly_large_overlap, $rere_set, $gene_ID_clusters, $threshold, $use_job_handler,$die_on_exception);
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
			  #my $ID = $pattern_model->pattern->wm_identifier; # VOLATILE
			  #$GU->user_info(2,"\'".$ID."\',\n");
			  push (@all_p_values, \@p_values);
			}
		      }
	    if ($stats_required) {
			$aggregate_exception->print_statistics_and_die;
	    }
	    return @all_p_values;
	} # hypergeometric_overrepresentation_test_in_parallel #	
	
	method private_hypergeometric_overrepresentation_test(Boolean $test_for_significantly_large_overlap, ReRe_Set $rere_set, ArrayRef[ArrayRef[Str]] $gene_ID_clusters, Str $threshold, Boolean $use_job_handler, Boolean $die_on_exception) {
	    # made decision on usage of Job_Handler a parameter here as
	    # a) there are Job_Handler-calls underneath this function, so at times it may be necessary to run this one plainly so
            #    as to see any problems more easily
	    # b) this function itself should never be called through the Job_Handler, so a redundant parameter cannot make a problem
            #    here
	    # returns a p-value
	    
	    # the variable die_on_exception allows us to handle errors from the job handler differently depending on whether we call 
	    # hypergeometric function for one pattern model or for multiple pattern models.
  		# die_on_exception = TRUE: no parallelisation; print statistics and die after a job is submitted
  		# die_on_exception = FALSE: allows parallelisation

	    my @all_gene_IDs = $rere_set->all_gene_ids();
	    if ($GU->list_has_redundancy(\@all_gene_IDs)) {
		die 'there may be a problem with your ReRe_Set, some gene-IDs occur more than once.';
	    }
	    
	    my @score_p_values;
	    if ($use_job_handler) {
		@score_p_values = $self->private_is_pattern_in_sequence_rere_set_through_job_handler($rere_set,$die_on_exception); # NOW P VALUES
	    } else {
		@score_p_values = $self->is_pattern_in_sequence_rere_set($rere_set); #NOW P VALUES
	    }
	    my @all_calls;
	    foreach my $p_value(@score_p_values){
	    	if ($p_value<$threshold) {
	    	my $result = TRUE;
	    	push (@all_calls,$result);
	    	}
	    	else{
	    		my $result = FALSE;
	    		push (@all_calls,$result);
	    	}
	    }
	#if ($p_value<$self->threshold) {
	#    $result = TRUE; # CHANGE THIS TO p-value to avoid thresholding issues.
	#}
	    my @positive_gene_IDs = $GU->subset_array_by_booleans(\@all_gene_IDs,\@all_calls);     
            # remove cluster-members that are not in universe (which is @all_gene_IDs)
        my @p_values;
        foreach my $gene_ID_cluster (@{$gene_ID_clusters}) {
			my @non_redundant_cluster = $GU->remove_duplicates_from_list($gene_ID_cluster);
			my @relevant_cluster_members = $GU->list_intersection(\@all_gene_IDs,\@non_redundant_cluster);
			# work out relevant set sizes
			my @cluster_members_with_pattern = $GU->list_intersection(\@positive_gene_IDs,\@relevant_cluster_members);
			my $universe_size = @all_gene_IDs;
			my $cluster_size = @relevant_cluster_members;
			my $positive_genes = @positive_gene_IDs;
			my $intersection_size = @cluster_members_with_pattern;
			# "translate" this into R's parameterisation of the hypergeometric test
			my $q = $intersection_size;
			my $m = $positive_genes;
			my $n = $universe_size-$positive_genes;
			my $k = $cluster_size;
			my $lower_tail;
			if (!$test_for_significantly_large_overlap) {
			$lower_tail = TRUE;
			} else {
			$lower_tail = FALSE;
			}
			my $stats = Statistics->new();
			my $p_value = $stats->phyper($q,$m,$n,$k,$lower_tail);
			if (($cluster_size == 0)||($positive_genes == 0)) {
			$p_value = 1;
			}
			push(@p_values,$p_value);
		}
	    return @p_values;
	} # private_hypergeometric_overrepresentation_test #

	method private_is_pattern_in_sequence_rere_set_through_job_handler(ReRe_Set $rere_set,Boolean $die_on_exception) {
	    my $cache = TRUE;
	    my $cache_duration = 200;
	    my $job_handler = Job_Handler->new();
	    $job_handler->get_config_settings();
	    my $function = 'is_pattern_in_sequence_rere_set';
	    push (my @parameters, $rere_set);
	    my $job_parameters = Job_Parameters->new(memory_requirement_high => TRUE,
						     wall_time_estimate => 172800
		);
	    my @cache_result = eval {
		$job_handler->handle_APPLES_function($function, $self, \@parameters, $cache, $cache_duration, $job_parameters);
	    };
	    my $aggregate_exception = Aggregate_Exception->new();
	    if ($@) {
		my $exception_content = $@;
		my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
		# no parallelism here, though
		$aggregate_exception = $outcome[1];
		if ($outcome[0]) {
			if ($die_on_exception){
			    $aggregate_exception->print_statistics_and_die;			    
			} else {
				die $aggregate_exception; ## Not sure about this line! ##
			}
		} else {
		    die 'did not expect this could happen at design time of this code.';
		}
	    }
	    shift @cache_result;
	    return @cache_result;
	} # private_is_pattern_in_sequence_rere_set_through_job_handler #

} # Generic_Pattern_Matching_Model #

class WM_Pattern_Model extends Generic_Pattern_Matching_Model {
    use Parameters;
    use General_Utilities;
    use Generic_Sequence_Pattern;
    use Genomic_Interval;
    use Genomic_Interval_Set;
    use Data::Dumper;
    use WM_Utilities;
    
    has +'pattern' => (is => 'ro', isa => 'Generic_Weight_Matrix', required => 1) ;
    has 'scoring_model' => (is => 'ro', isa => 'Generic_WM_Scoring_Parameters', required => 1);

    
    my $GU = General_Utilities->new();
    my $wm_util = WM_Utilities->new();
    
    override number_of_possible_sites(Genomic_Interval_Set $gi_set) {		
		
		my $wm_length = $self->pattern->get_pssm_length();
		my $possible_sites_to_score = 0;
		$GU->user_info(3,"LENGTH: $wm_length"."\n");
		foreach my $gi ( @{$gi_set->genomic_interval_set} ) {			
			my $gi_length = $gi->get_length();		

			$GU->user_info(3,"GI-LENGTH: $gi_length\n");

			my $max_sites = 2*($gi_length-$wm_length+1);		
			$possible_sites_to_score += $max_sites;
		}
		
		return $possible_sites_to_score;
		
    } # number_of_possible_sites #
    
    override is_pattern_in_sequence(Genomic_Interval_Set $gi_set) {	
	# need to add an attribute for threshold to support this method in this class
	
	die 'method not implemented yet';
    } # is_pattern_in_sequence #
    
    override find_top_n_matches_in_sequence(Genomic_Interval_Set $gi_set, Int $n) {
		
		my $raw_result = $wm_util->score_wm_on_sequences($self, $gi_set);
		my $wm_scores = WM_Scores->new(raw_result => $raw_result, gi_set => $gi_set);
		my $wm_length = $self->pattern->get_pssm_length();
		my $result = $wm_scores->get_n_nonoverlapping_sites($n, $wm_length);
		
		return $result;
		
    } # find_top_n_matches_in_sequence #
    
} # WM_Pattern_Model #

class Multi_Occurrence_WM_Pattern_Model extends Generic_Pattern_Matching_Model {
    # uses binomial overrepresentation test to test 1 or more occurrences of the same 
    # weight matrix for significance taking total length of sequence tested into account
    
    use Parameters;
    use APPLES_Datatypes qw (Boolean NonNegativeNum PositiveInt);
    use constant {FALSE => 0,
		  TRUE	=> 1};	
    use General_Utilities;
    use Generic_Sequence_Pattern;
    use Genomic_Interval_Set;

    has +'pattern' => (is => 'ro', isa => 'Generic_Weight_Matrix', required => 1);
    has 'scoring_model' => (is => 'ro', isa => 'Generic_WM_Scoring_Parameters', required => 1);
    has 'maximum_number_of_sites' => (is => 'ro', isa => PositiveInt, required => 1);
    #has 'threshold' => (is => 'ro', isa => 'NonNegativeNum', required => 1);

    method is_pattern_in_sequence(Genomic_Interval_Set $gi_set) {
	# returns a Boolean indicating presence or absence of pattern in $gi_set
	
	my $wm_pattern_model = WM_Pattern_Model->new(pattern => $self->pattern, scoring_model => $self->scoring_model);
	my $overrepresentation_test_result = $wm_pattern_model->binomial_overrepresentation_test($gi_set,
	    $self->maximum_number_of_sites);
	my $p_value = $overrepresentation_test_result->{PVALUE};
	#my $GU->user_info(3,"P value: $p_value\n");
	return $p_value; 
    } # is_pattern_in_sequence #

} # Multi_Occurrence_WM_Pattern_Model #

class Pair_WM_Pattern_Model extends Generic_Pattern_Matching_Model {
	use Genomic_Interval;
	use Genomic_Interval_Set;
	use BiFa_Raw_Result;
	use Data::Dumper;
	use General_Utilities;
	use WM_Utilities;
	use Generic_Sequence_Pattern;
	#use Generic_Pattern_Matching_Model;
	use APPLES_Datatypes qw(ScoringModel);
		
	has +'pattern' => (is => 'ro', isa => 'Generic_Pair_Weight_Matrices', required => 1);
	has 'scoring_model' => (is => 'rw', isa => ScoringModel, required => 1);
	
	my $GU = General_Utilities->new();

	method score (Genomic_Interval_Set $gi_set){
		$| = 1;# Buffering
		my $wm_util = WM_Utilities->new();
		
		# Get scores for the two individual weight matrices
		# Weight matrix A
		my $wmpmA = WM_Pattern_Model->new(pattern => $self->pattern->wm_a, scoring_model => $self->scoring_model);
		my $raw_result_A = $wm_util->score_wm_on_sequences($wmpmA, $gi_set);
		my $wm_scores_A = WM_Scores->new(raw_result => $raw_result_A, gi_set => $gi_set);
		my $positive_strand_scores_A = $wm_scores_A->get_positive_strand_scores();
		my $negative_strand_scores_A = $wm_scores_A->get_negative_strand_scores();
		# Weight matrix B
		my $wmpmB = WM_Pattern_Model->new(pattern => $self->pattern->wm_b, scoring_model => $self->scoring_model);
		my $raw_result_B = $wm_util->score_wm_on_sequences($wmpmB, $gi_set);
		my $wm_scores_B = WM_Scores->new(raw_result => $raw_result_B, gi_set => $gi_set);
		my $positive_strand_scores_B = $wm_scores_B->get_positive_strand_scores();
		my $negative_strand_scores_B = $wm_scores_B->get_negative_strand_scores();
		#$GU->user_info(3,Dumper($positive_strand_scores_B));
		my $L_wmA = $wm_util->get_wm_length($self->pattern->wm_a);
		my $L_wmB = $wm_util->get_wm_length($self->pattern->wm_b);
		#$GU->user_info(3,"$L_wmA:$L_wmB\n");
		#$GU->user_info(3,Dumper($negative_strand_scores_A));
		my $i=0;
		foreach (@$positive_strand_scores_A) {
			foreach(@$_) {
				$GU->user_info(3,"$i: $_\n");
				$i++;
			}
		}
		#$GU->user_info(3,"\n");
		#$GU->user_info(3,Dumper($negative_strand_scores_B));
		
		# Cycle through each genomic interval in the set
		#foreach my $gi ( @{$gi_set->genomic_interval_set} ){
		my @all_scores_for_gi_set;
		for (my $x = 0; $x < scalar(@{$gi_set->genomic_interval_set}); $x++) {
			
			my $gi = ${$gi_set->genomic_interval_set}[$x];
			#$GU->user_info(3,$gi."\n");
			#$GU->user_info(3,"Testing: ".$gi->label."\n");
			my $gi_counter = $x;
			#$GU->user_info(3,$$positive_strand_scores_A[$x]);
			#$GU->user_info(3,"\n");
			#my $total_possible_pair_scores = ( length($gi->get_sequence()) - ($L_wmA + $self->pattern->max_spacing + $L_wmB) );
			#$GU->user_info(3,@{$$positive_strand_scores_A[$x]}."\n");
			#$GU->user_info(3,@{$$positive_strand_scores_B[$x]}."\n");
			# Score pair weight matrices on the +ve strand with wm A being the central pattern
			#my $LF_scores_S = $self->private_score_pair_wm_left_flank($$positive_strand_scores_A[$x], $$positive_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()), 1, $gi_counter, $gi->strand);
			my $RF_scores_S = $self->private_score_pair_wm_right_flank($$positive_strand_scores_A[$x], $$positive_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()), 1, $gi_counter, $gi->strand);
			#$GU->user_info(3,Dumper($RF_scores_S));
			# Score pair weight matrices on the -ve strand with wm A being the central pattern
			#my $LF_scores_S = $self->private_score_pair_wm_left_flank($$negative_strand_scores_A[$x], $$negative_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()), -1, $gi_counter);
			#my $RF_scores_AS = $self->private_score_pair_wm_right_flank_REVCOM($$negative_strand_scores_A[$x], $$negative_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()), -1, $gi_counter, $gi->strand);
			#my $LF_scores_AS = $self->private_score_pair_wm_left_flank_REVCOM($$negative_strand_scores_A[$x], $$negative_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()), -1, $gi_counter, $gi->strand);
			
			#$GU->user_info(3,Dumper($LF_scores_AS));
			#my @all_scores = [$RF_scores_S];#, $RF_scores_AS];#, $RF_scores_S];
			my @all_scores = [$RF_scores_S];#, $LF_scores_AS];#, $LF_scores_AS];
			#$GU->user_info(3,Dumper($RF_scores_AS));
			
			# Score pair weight matrices on the +ve strand with wm B being the central pattern
			#my ($LF_scores2, $RF_scores2) = $self->private_score_pair_wm($$positive_strand_scores_B[$x], $$positive_strand_scores_A[$x], $L_wmB, $L_wmA, length($gi->get_sequence()));
			# Score pair weight matrices on the -ve strand with wm A being the central pattern
			#my ($LF_scores3, $RF_scores3) = $self->private_score_pair_wm($$negative_strand_scores_A[$x], $$negative_strand_scores_B[$x], $L_wmA, $L_wmB, length($gi->get_sequence()));
			# Score pair weight matrices on the -ve strand with wm B being the central pattern
			#my ($LF_scores4, $RF_scores4) = $self->private_score_pair_wm($$negative_strand_scores_B[$x], $$negative_strand_scores_A[$x], $L_wmB, $L_wmA, length($gi->get_sequence()));
			#$GU->user_info(3,Dumper ($LF_scores1));
			push(@all_scores_for_gi_set, \@all_scores);
			#$GU->user_info(3,Dumper($LF_scores_AS));
		}
		#$GU->user_info(3,@all_scores_for_gi_set);
		#$GU->user_info(3,"\n");
		my $pair_wm_scores = Pair_WM_Scores->new(raw_result => \@all_scores_for_gi_set, gi_set => $gi_set);
		$pair_wm_scores->process_raw_result();
		
		my $result = $pair_wm_scores->get_n_nonoverlapping_sites(1);
		$GU->user_info(3,"**********************************************************\n");
		$GU->user_info(3,Dumper($result));
		#my (@LF_scores, @RF_scores) = private_score_pair_wm(@$positive_strand_scores_A, @$positive_strand_scores_B);
			# Score pair weight matrices on the +ve strand with wm B being the central pattern
			
			# Score pair weight matrices on the -ve strand with wm A being the central pattern
			
			# Score pair weight matrices on the -ve strand with wm B being the central pattern
			
		#}
		
		#my $num_intervals = @{$gi_set->genomic_interval_set};	
	} # score #
	
	override is_pattern_in_sequence(Genomic_Interval_Set $gi_set) {
		die "Not implemented yet\n";
	} # is_pattern_in_sequence #
	
	override find_top_n_matches_in_sequence(Genomic_Interval_Set $gi_set, Int $n) {
		die "Not implemented yet\n";
	} # find_top_n_matches_in_sequence #
	
	override number_of_possible_sites(Genomic_Interval_Set $gi_set) {
		die "Not implemented yet\n";
	} # number_of_possible_sites #

	method private_score_pair_wm_left_flank_REVCOM (ArrayRef $centre, ArrayRef $flank, Int $centre_pattern_length, Int $flanking_pattern_length, Int $sequence_length, Int $strand, Int $gi_counter, Str $source_strand){
		
		my @L_SCORES;
		#my @R_SCORES;
		
		my @reverse_centre = reverse(@$centre);
		@{$centre} = ();
		@{$centre} = @reverse_centre;
		my @reverse_flank = reverse(@$flank);
		@{$flank} = ();
		@{$flank} = @reverse_flank;
		
		my $strand_label;
		if ($strand == 1) {
			$strand_label = 'positive';	
		}
		elsif ($strand == -1) {
			$strand_label = 'negative';
		}
		else {
			die "Strand label not defined correctly\n";
		}
		my $total_possible_pair_scores = $sequence_length;
		my $count = 0;
		for (my $i = 0; $i < $total_possible_pair_scores - $centre_pattern_length + 1; $i++) {
			#$GU->user_info(3,("\r $i"));
			my $centre_i_score = $$centre[$i];
			my $offset = ($i - ( $flanking_pattern_length + $self->pattern->max_spacing ));
			my $length = ($self->pattern->max_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			#$GU->user_info(3,"$i > $centre_i_score\n");
			if( $i < ( $flanking_pattern_length + $self->pattern->max_spacing ) ) {
				#$score = { SCORE => 1, FLANKING_WM_INDEX => 'NA', STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				#FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i };
				#$count++;
				#$GU->user_info(3,"Index: $i\tNA\n");
			}
			else {
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				my $flanking_wm_index = $offset + $score_position;
				#$GU->user_info(3,"index: $i\tscore -> $score\tf_index -> $score_position\n");
				if ($flanking_wm_index < 0){
					die "oops\n";
				}
				
				$score = { SCORE => $score, FLANKING_WM_INDEX => ( $sequence_length - ($flanking_wm_index + $flanking_pattern_length) ) , STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => ( $sequence_length - ( $i + $centre_pattern_length ) ), SOURCESTRAND => $source_strand};
				push (@L_SCORES, $score);
			}
			
			#if( !defined($score) ){
			#	die;
			#}
			#else{
			#	push (@L_SCORES, $score);
			#}
			
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		return \@L_SCORES;
	} # private_score_pair_wm_left_flank_REVCOM #

	method private_score_pair_wm_right_flank_REVCOM (ArrayRef $centre, ArrayRef $flank, Int $centre_pattern_length, Int $flanking_pattern_length, Int $sequence_length, Int $strand, Int $gi_counter,  Str $source_strand){
		
		#my @L_SCORES;
		my @R_SCORES;
		
		my @reverse_centre = reverse(@$centre);
		@{$centre} = ();
		@{$centre} = @reverse_centre;
		my @reverse_flank = reverse(@$flank);
		@{$flank} = ();
		@{$flank} = @reverse_flank;
		
		
		my $strand_label;
		if ($strand == 1) {
			$strand_label = 'positive';	
		}
		elsif ($strand == -1) {
			$strand_label = 'negative';
		}
		else {
			die "Strand label not defined correctly\n";
		}
		# Compute combined scores for central matrix and right hand flank
		#my $total_possible_pair_scores = ( $sequence_length - ($centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) );
		my $total_possible_pair_scores = $sequence_length;
		my $count = 0;
		for (my $i = 0; $i < $total_possible_pair_scores; $i++) {
			#$GU->user_info(3,("\r $i"));
			my $centre_i_score = $$centre[$i];
			my $offset = ($i + $centre_pattern_length);
			my $length = ($self->pattern->max_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			if ( ($i + $centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) > $total_possible_pair_scores) {
				#$score = { SCORE => 1, FLANKING_WM_INDEX => 'NA', STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				#FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i};
				#$count++;
				#$GU->user_info(3,"index: $i\tNA\n");
			}
			else {
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				#$GU->user_info(3,"index: $i\tscore -> $score\tf_index -> $score_position\n");
				# Need to define the position of the flanking wm that gave the pair with the best combined score
				
				my $flanking_wm_index = $i + $centre_pattern_length + $score_position;
				$score = { SCORE => $score, FLANKING_WM_INDEX => ( $sequence_length - ($flanking_wm_index + $flanking_pattern_length) ) , STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => ($sequence_length - ($i + $centre_pattern_length) ), SOURCESTRAND => $source_strand };
				
				push (@R_SCORES, $score);
				
			}
			
			#if(!defined($score)){
				#	die;
				#}
			#else{
				#	push (@R_SCORES, $score);
				#}
			
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		return \@R_SCORES;
	} # private_score_pair_wm_right_flank_REVCOM #
	
	method private_score_pair_wm (ArrayRef $centre, ArrayRef $flank, Int $centre_pattern_length, Int $flanking_pattern_length, Int $sequence_length, Int $strand){
	
		my @L_SCORES;
		my @R_SCORES;
		
		my $strand_label;
		if ($strand == 1) {
			$strand_label = '+';	
		}
		elsif ($strand == -1) {
			$strand_label = '-';
		}
		else {
			die "Strand label not defined correctly\n";
		}
		# Compute combined scores for central matrix and right hand flank
		#my $total_possible_pair_scores = ( $sequence_length - ($centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) );
		my $total_possible_pair_scores = $sequence_length;
		my $count = 0;
		for (my $i = 0; $i < $total_possible_pair_scores; $i++){
			#$GU->user_info(3,"\r $i");
			my $centre_i_score = $$centre[$i];
			my $offset = ($i + $centre_pattern_length);
			my $length = ($self->pattern->max_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			if ( ($i + $centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) > $total_possible_pair_scores) {
				$score = { SCORE => 'NA', FLANKING_WM_INDEX => 'NA', STRAND => $strand_label };
				$count++;
			}
			else{
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				
				# Need to define the position of the flanking wm that gave the pair with the best combined score
				
				my $flanking_wm_index = $i + $centre_pattern_length + $score_position;
				$score = { SCORE => $score, FLANKING_WM_INDEX => $flanking_wm_index, STRAND => $strand_label };
				
			}
			
			if (!defined($score)) {
				die;
			}
			else{
				push (@R_SCORES, $score);
			}
		
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		# Compute combined scores for central matrix and left hand flank
		
		for (my $i = 0; $i < $total_possible_pair_scores - $centre_pattern_length + 1; $i++){
			#$GU->user_info(3,"\r $i");
			my $centre_i_score = $$centre[$i];
			my $offset = ($i - ( $flanking_pattern_length + $self->pattern->max_spacing ));
			my $length = ($self->pattern->max_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			#$GU->user_info(3,$centre_i_score." | ");
			if ( $i < ( $flanking_pattern_length + $self->pattern->max_spacing ) ) {
				$score = { SCORE => 'NA', FLANKING_WM_INDEX => 'NA', STRAND => $strand_label };
				$count++;
			}
			else {
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				my $flanking_wm_index = $offset + $score_position;
				
				
				$score = { SCORE => $score, FLANKING_WM_INDEX => $flanking_wm_index , STRAND => $strand_label};
			}
			
			if ( !defined($score) ) { 
				die;
			}
			else {
				push (@L_SCORES, $score);
			}
			
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		#foreach(@L_SCORES){
		#$GU->user_info(3,substr($_,0,6)."\t");
		#}
		#$GU->user_info(3,"\n");
		#$GU->user_info(3,$count."\n");
		return (\@L_SCORES, \@R_SCORES);
		
	} # private_score_pair_wm #
	
	method private_score_pair_wm_left_flank (ArrayRef $centre, ArrayRef $flank, Int $centre_pattern_length, Int $flanking_pattern_length, Int $sequence_length, Int $strand, Int $gi_counter, Str $source_strand){
		
		my @L_SCORES;
				
		my $strand_label;
		if ($strand == 1){
			$strand_label = 'positive';	
		}
		elsif ($strand == -1){
			$strand_label = 'negative';
		}
		else{
			die "Strand label not defined correctly\n";
		}
		my $total_possible_pair_scores = $sequence_length;
		my $count = 0;
		for (my $i = 0; $i < $total_possible_pair_scores - $centre_pattern_length + 1; $i++) {
			#$GU->user_info(3,("\r $i"));
			my $centre_i_score = $$centre[$i];
			my $offset = ($i - ( $flanking_pattern_length + $self->pattern->max_spacing ));
			my $length = ($self->pattern->max_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			#$GU->user_info(3,$centre_i_score." | ");
			if ( $i < ( $flanking_pattern_length + $self->pattern->max_spacing ) ) {
				#$score = { SCORE => 1, FLANKING_WM_INDEX => 'NA', STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				#FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i };
				#$count++;
				#$GU->user_info(3,"Index: $i\tNA\n");
			}
			else {
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				my $flanking_wm_index = $offset + $score_position;
				#$GU->user_info(3,"index: $i\tscore -> $score\tf_index -> $score_position\n");
				if ($flanking_wm_index < 0){
						die "oops\n";
				}
				$score = { SCORE => $score, FLANKING_WM_INDEX => $flanking_wm_index , STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
							FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i, SOURCESTRAND => $source_strand};
				push (@L_SCORES, $score);
			}
			
			#if( !defined($score) ){
			#	die;
			#}
			#else{
			#	push (@L_SCORES, $score);
			#}
			
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		return \@L_SCORES;
	} # private_score_pair_wm_left_flank #
	
	method private_score_pair_wm_right_flank (ArrayRef $centre, ArrayRef $flank, Int $centre_pattern_length, Int $flanking_pattern_length, Int $sequence_length, Int $strand, Int $gi_counter, Str $source_strand){
		
		my @L_SCORES;
		my @R_SCORES;
		
		my $strand_label;
		if ($strand == 1) {
			$strand_label = 'positive';	
		}
		elsif ($strand == -1) {
			$strand_label = 'negative';
		}
		else {
			die "Strand label not defined correctly\n";
		}
		# Compute combined scores for central matrix and right hand flank
		#my $total_possible_pair_scores = ( $sequence_length - ($centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) );
		my $total_possible_pair_scores = $sequence_length;
		my $count = 0;
		for (my $i = 0; $i < $total_possible_pair_scores; $i++) {
			#$GU->user_info(3,"\r $i");
			my $centre_i_score = $$centre[$i];
			my $offset = ($i + $centre_pattern_length + $self->pattern->min_spacing);
			my $length = ($self->pattern->max_spacing - $self->pattern->min_spacing);
			#$GU->user_info(3,"OFFSET $offset, Length = $length\n");
			my $score;
			my $score_position;
			if ( ($i + $centre_pattern_length + $self->pattern->max_spacing + $flanking_pattern_length) > $total_possible_pair_scores) {
				#$score = { SCORE => 1, FLANKING_WM_INDEX => 'NA', STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				#FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i};
				#$count++;
				#$GU->user_info(3,"index: $i\tNA\n");
			}
			else {
				my @flank_possible_scores = @$flank[ $offset .. ($offset + $length) ];
				#$GU->user_info(3,@flank_possible_scores);
				#$GU->user_info(3,"\n");
				($score, $score_position) = $self->private_min_score($centre_i_score, \@flank_possible_scores);
				#$GU->user_info(3,"index: $i\tscore -> $score\tf_index -> $score_position\n");
				# Need to define the position of the flanking wm that gave the pair with the best combined score
				
				my $flanking_wm_index = $i + $centre_pattern_length + $score_position;
				$score = { SCORE => $score, FLANKING_WM_INDEX => $flanking_wm_index, STRAND => $strand_label, CENTRAL_WM_LEN => $centre_pattern_length,
				FLANKING_WM_LEN => $flanking_pattern_length, GI_INDEX => $gi_counter, CENTRAL_WM_INDEX => $i, SOURCESTRAND => $source_strand};
				
				push (@R_SCORES, $score);
				
			}
			
			#if(!defined($score)){
			#die;
			#}
			#else{
			#	push (@R_SCORES, $score);
			#}
			
			#my $score = $self->private_min_score($centre_i_score, \@flank_possible_scores);
			#$GU->user_info(3,$score.",");
		}
		
		return \@R_SCORES;
	} # private_score_pair_wm_right_flank #
		
	method private_min_score (Num $wma_score, ArrayRef $wmb) {
		
		my $min_score;
		my $min_score_position = 0;
		#$GU->user_info(3,"->".$#{$wmb}."\n");
		#$GU->user_info(3,"c_score => $wma_score\n");
		for (my $m = 0; $m <= $#{$wmb}; $m++ ) {
			#$GU->user_info(3,"\t$$wmb[$m]\n");
			my $cs = $self->private_combined_score($wma_score, $$wmb[$m], $m+1);
			
			if (!defined($min_score)) {
				$min_score = $cs;
				$min_score_position = $m;
			}
			elsif ( $cs < $min_score) {
				$min_score = $cs;
				$min_score_position = $m;
			}
			else {
				# The min score stays the same
			}
			
		}
		#$GU->user_info(3,"ms ".$$wmb[$min_score_position]."\n");
		return ($min_score, $min_score_position);
	} # private_min_score #
	
	method private_combined_score (Num $p1, Num $p2, Int $m) {
		
		# Compute both and return minimum
		#$GU->user_info(3,"p1 = $p1 | $p2 | $m\n");
	
		my $cs1 = $p1 * ( (1 - (1 - $p2)**$m) );	
		my $cs2 = $p2 * ( (1 - (1 - $p1)**$m) );	
		if ($cs1 <= $cs2) {
		    return $cs1;
		}
		else {
		    return $cs2;
		}
	} # private_combined_score #
	
} # Pair_WM_Pattern_Model #


