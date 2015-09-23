### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### WM_Scores Class ###

use MooseX::Declare;

class WM_Scores {	
	use Data::Dumper;
	use General_Utilities;
	use APPLES_Datatypes qw (Probability);
	use Genomic_Interval;
	use BiSi_Set;
	use constant {FALSE => 0,
		      TRUE  => 1};
	
	has 'raw_result' => (is => 'ro', isa => 'ArrayRef', required => 1, trigger => \&create_processed_result);
	has 'gi_set' => (is => 'ro', isa => 'Genomic_Interval_Set', required => 1);
	has 'processed_result' => (is => 'rw', isa => 'ArrayRef');
	
	my $GU = General_Utilities->new();

	sub create_processed_result(){	
		my ($self, $arg) = @_;	
		$self->processed_result( $self->process_raw_result() );
		#$GU->user_info(3,Dumper($self->gi_set));	
	} # create_processed_result #
	
	method get_positive_strand_scores() {	
		my @all_positive_strand_scores;
		my $l = scalar( @{$self->raw_result} );
		
		for (my $i = 0; $i < scalar(@{$self->raw_result}); $i ++) {
			
			my @raw_scores = @{ ${$self->raw_result}[$i] };
			my @positive_strand_scores;
			for (my $j = 0; $j < scalar(@raw_scores); $j += 2) {
				#$GU->user_info(3,${$self->raw_result}[$i]);
				#$GU->user_info(3,"\n");
				push(@positive_strand_scores, $raw_scores[$j]);
			
			}
			push(@all_positive_strand_scores, \@positive_strand_scores);
		}
		#$GU->user_info(3,"elements = " .$#all_positive_strand_scores."\n");
		return \@all_positive_strand_scores;
		
	} # get_positive_strand_scores #
	
	method get_negative_strand_scores() {	
		my @all_negative_strand_scores;
		my $l = scalar( @{$self->raw_result} );
		
		for (my $i = 0; $i < scalar(@{$self->raw_result}); $i ++) {
			
			my @raw_scores = @{ ${$self->raw_result}[$i] };
			my @negative_strand_scores;
			for (my $j = 0; $j < scalar(@raw_scores); $j += 2) {
				#$GU->user_info(3,${$self->raw_result}[$i]);
				#$GU->user_info(3,"\n");
				push(@negative_strand_scores, $raw_scores[$j+1]);
				
			}
			push (@all_negative_strand_scores, \@negative_strand_scores);
		}
			
		return \@all_negative_strand_scores;
	} # get_negative_strand_scores #
	
	method process_raw_result() {
		my @ordered_sites_one_result;
		my @sites_all_results;
		my $gi_counter = 0;
		
		my @all_raw_scores = @{$self->raw_result};
		
		for my $i (0..@all_raw_scores-1) {
			
			if (exists($all_raw_scores[$i])) {
				my $s = $all_raw_scores[$i];
				my @one_result_raw_scores = @{$s};
				
				#$GU->user_info( 2, "going to get ordered sites...\n" );
				#$GU->user_info(3,"Going to get ordered sites.....\n");
				#$GU->user_info(3,"GENOMIC_INTERVAL: ". ${$self->gi_set->genomic_interval_set}[$i]{strand}."\n");
				my $genomic_interval_strand_in_db = ${$self->gi_set->genomic_interval_set}[$i]{strand};
				@ordered_sites_one_result = $self->all_sites_ordered_by_score(\@one_result_raw_scores, $gi_counter, $genomic_interval_strand_in_db);
				#$GU->user_info( 3, Dumper(@ordered_sites_one_result) );
				#$GU->user_info(3,Dumper(@ordered_sites_one_result));
				push (@sites_all_results, @ordered_sites_one_result);
			}
			else {
				#$GU->user_info( 3, "no scores for this sequence!\n" );
				$GU->user_info(3,"No Scores for this sequence!\n");
			}
			
			$gi_counter++;
			
		}
		#$GU->user_info( 3, Dumper (@sites_all_results));
		my @ordered_sites_all_results = $self->ordered_sites_all_results(\@sites_all_results);
		return \@ordered_sites_all_results;
		#$GU->user_info( 3, Dumper (@ordered_sites_all_results) );
	} # process_raw_result #
	
	method process_raw_result {
		my @ordered_sites_one_result;
		my @sites_all_results;
		my $gi_counter = 0;
		
		my @all_raw_scores = @{$self->raw_result};
		
		for my $i (0..@all_raw_scores-1) {
			
			if (exists($all_raw_scores[$i])) {
				my $s = $all_raw_scores[$i];
				my @one_result_raw_scores = @{$s};
				
				#$GU->user_info( 2, "going to get ordered sites...\n" );
				#$GU->user_info(3,"Going to get ordered sites.....\n");
				#$GU->user_info(3,"GENOMIC_INTERVAL: ". ${$self->gi_set->genomic_interval_set}[$i]{strand}."\n");
				my $genomic_interval_strand_in_db = ${$self->gi_set->genomic_interval_set}[$i]{strand};
				@ordered_sites_one_result = $self->all_sites_ordered_by_score(\@one_result_raw_scores, $gi_counter, $genomic_interval_strand_in_db);
				#$GU->user_info( 3, Dumper(@ordered_sites_one_result) );
				#$GU->user_info(3,Dumper(@ordered_sites_one_result));
				push (@sites_all_results, @ordered_sites_one_result);
			}
			else {
				$GU->user_info(3,"No Scores for this sequence!\n");
			}
			
			$gi_counter++;
			
			
		}
		#$GU->user_info( 3, Dumper (@sites_all_results) ); 
		my @ordered_sites_all_results = $self->ordered_sites_all_results(\@sites_all_results);
		return \@ordered_sites_all_results;
		#$GU->user_info( 3, Dumper (@ordered_sites_all_results) );
	} # process_raw_result #
	
	method all_sites_ordered_by_score (ArrayRef $one_result_ref, Int $gi_counter, Str $source_strand) {
		
		my @raw_result = @{$one_result_ref};
		
		my @positive_strand_scores;
		my @negative_strand_scores;
		my $strand = 'NONE';
		my $rec;
		my @recs1;
		my @recs2;
		my @recs;
		
		#$GU->user_info(3,Dumper($one_result_ref));
		for (my $i = 0; $i < scalar(@raw_result); $i += 2) {
			push(@positive_strand_scores, $raw_result[$i]);
			push(@negative_strand_scores, $raw_result[$i+1]);
		}
		for ( my $i = 0; $i < scalar(@positive_strand_scores); $i++ ) {
			#$GU->user_info(3,"DEBUG: $i\n");
			#if($strand eq 'positive'){
			#	$strand = 'positive';
			#}
			#elsif($strand eq 'negative'){
			#	$GU->user_info(3,"hello\n");
			#	$strand = 'negative';
			#}
			#else{
			#	die "strand type error!\n";
			#}
			
			$strand = 'positive';
			$rec = {
				INTERVALINDEX => $gi_counter,
				SITE   => $i,
				STRAND => $strand,
				SOURCESTRAND => $source_strand,
				PVALUE => $positive_strand_scores[$i]
			};
			
			if ( defined($positive_strand_scores[$i]) ){
				
				
			}
			else{
				
			}
			
			push (@recs1, $rec);
		}
		
		for ( my $i = 0; $i <  scalar(@negative_strand_scores); $i++ ) {
			
			#if($strand eq 'positive'){
			#	$strand = 'positive';
			#}
			#elsif($strand eq 'negative'){
			#	$strand = 'negative';
			#}
			#else{
			#	die "strand type error!\n";
			#}
			
			$strand = 'negative';
			$rec = {	
				INTERVALINDEX => $gi_counter,
				SITE   => $i,
				STRAND => $strand,
				SOURCESTRAND => $source_strand,
				PVALUE => $negative_strand_scores[$i]
			};
			push (@recs2, $rec);
			if ( defined($negative_strand_scores[$i]) ){
				
				
			}
			else{
				
			}
		}
		
		@recs = (@recs1, @recs2);
		
		my @sorted_recs = (sort { $a->{PVALUE} <=> $b->{PVALUE} } @recs); # smallest->biggest
		
		return @sorted_recs;
	} #all_sites_ordered_by_score
	
	method ordered_sites_all_results (ArrayRef $sites_all_results_ref) {
		
		my @sites_all_results = @{$sites_all_results_ref};
		
		my @sorted_result = (sort { $a->{PVALUE} <=> $b->{PVALUE} } @sites_all_results);
		
		return @sorted_result;
		
	} #ordered_sites_all_results
	
	method private_site_to_gi(Any $site, Int $length) {
		
		my $five_pos;
        my $three_pos;
        my $STRAND;
		
		my $genomic_interval = (@{$self->gi_set->genomic_interval_set}[$site->{INTERVALINDEX}]);
            
            if ($site->{STRAND} eq 'positive') {
                
                if($site->{SOURCESTRAND} eq 'positive'){
                    $five_pos = $genomic_interval->five_prime_pos+$site->{SITE};
                    $three_pos =  $genomic_interval->five_prime_pos+$site->{SITE}+$length -1;
                    $STRAND = 'positive';
                }
                elsif($site->{SOURCESTRAND} eq 'negative'){
                    $five_pos = $genomic_interval->five_prime_pos-$site->{SITE};
                    $three_pos =  $genomic_interval->five_prime_pos-$site->{SITE} - $length +1;
                    $STRAND = 'negative';
                }
                else{
                    die;
                }
                
            }
            elsif ($site->{STRAND} eq 'negative') {
                
                if($site->{SOURCESTRAND} eq 'positive'){
                    $three_pos = $genomic_interval->five_prime_pos+$site->{SITE};
                    $five_pos =  $genomic_interval->five_prime_pos+$site->{SITE}+$length -1;
                    $STRAND = 'negative';
                }
                elsif($site->{SOURCESTRAND} eq 'negative'){
                    $three_pos = $genomic_interval->five_prime_pos-$site->{SITE};
                    $five_pos =  $genomic_interval->five_prime_pos-$site->{SITE} - $length +1;
                    $STRAND = 'positive';
                }
                else{
                    die;
                }
                
            }     

            my $gi = Genomic_Interval->new(
                genome_db => $genomic_interval->genome_db,
                coord_sys_name => $genomic_interval->coord_sys_name,
                region => $genomic_interval->region,
                five_prime_pos => $five_pos,
                three_prime_pos => $three_pos,
                strand => $STRAND,
                working_sequence => $genomic_interval->working_sequence,
                rel_to_tss => $site->{SITE} );
                
            if (defined $genomic_interval->label) {
                $gi->label($genomic_interval->label);
            }
            
            $gi->get_sequence();
            
            return $gi;
		
	}
	
	method get_n_nonoverlapping_sites (Int $n_sites, Int $length) {
		
		my @ordered_sites = @{$self->processed_result};
		
		my @non_overlapping_bisis;
		if ($#ordered_sites<0) {
		    my $bisi_set = BiSi_Set->new(bisi_set => \@non_overlapping_bisis);
		    return $bisi_set;
		}

		# make best site into gi->bisi, push into array.
		my @non_overlapping_gi_set;
			
		my @gi_array;
		my $gi = $self->private_site_to_gi($ordered_sites[0], $length);
        push(@gi_array, $gi);
            
        my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
			
			
		my $best_bisi = BiSi->new(genomic_interval_set=>$gi_set, single_site_pvalue=>$ordered_sites[0]->{PVALUE});
		push (@non_overlapping_bisis, $best_bisi);
			
		push (@non_overlapping_gi_set, $gi);
		
		#$GU->user_info(3,$gi->{five_prime_pos}."\t".$gi->{three_prime_pos}."\t".$gi->{gi_sequence}."\n");
			
		
		my $counter = 1;
		
		for (my $i = 1; $i < scalar(@ordered_sites); $i ++) {		
			if ($counter < $n_sites) {		
				
				my @gi_array;
				my $gi = $self->private_site_to_gi($ordered_sites[$i], $length);
        		push(@gi_array, $gi);
            
                my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
				
				my $make_bisi = TRUE;
				
				foreach my $binding_site (@non_overlapping_bisis) {

					my $overlap = $gi_set->gi_set_overlap($binding_site->genomic_interval_set);
					if ($overlap) {
						$make_bisi = FALSE;
						#$GU->user_info(3,$gi->{five_prime_pos}."\t".$gi->{three_prime_pos}."\t".$gi->{gi_sequence}."\tOVERLAPS with ".
						#${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{five_prime_pos}."\t".${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{three_prime_pos}.
						#"\t".${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{gi_sequence}."\n");
					}
					
				}
				
				if ($make_bisi) {
					
					my $bisi = BiSi->new(genomic_interval_set=>$gi_set, single_site_pvalue=>$ordered_sites[$i]->{PVALUE});
					
					push (@non_overlapping_bisis, $bisi);
					$counter++;
				}
				
			}
			else{
				last;
			}
			
			
		}
		if ($counter == $n_sites) {
			#$GU->user_info( 3, "\n\nn sites found!\n" );
		}
		
		# returns array of n best nonoverlapping sites, as an array of BiSis
		my $bisi_set = BiSi_Set->new(bisi_set => \@non_overlapping_bisis);
		return $bisi_set;

	} # get_n_nonoverlapping_sites #
	
	method get_nonoverlapping_sites (Probability $threshold, Int $length) {
        
            my @ordered_sites = @{$self->processed_result};
        
        my @non_overlapping_bisis;
        if ($#ordered_sites<0) {
            my $bisi_set = BiSi_Set->new(bisi_set => \@non_overlapping_bisis);
            return $bisi_set;
        }

        # make best site into gi->bisi, push into array.
        my @non_overlapping_gi_set;
            
        my @gi_array;
        my $gi = $self->private_site_to_gi($ordered_sites[0], $length);
        push(@gi_array, $gi);
            
        my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
               
        my $best_bisi = BiSi->new(genomic_interval_set=>$gi_set, single_site_pvalue=>$ordered_sites[0]->{PVALUE});
        push (@non_overlapping_bisis, $best_bisi);
            
        push (@non_overlapping_gi_set, $gi);
        
        #$GU->user_info(3,$gi->{five_prime_pos}."\t".$gi->{three_prime_pos}."\t".$gi->{gi_sequence}."\n");
            
        
        for (my $i = 1; $i < scalar(@ordered_sites); $i ++) {  
        	
        	if($ordered_sites[$i]->{PVALUE} <= $threshold) {
     
                my @gi_array;
                my $gi = $self->private_site_to_gi($ordered_sites[$i], $length);
                push(@gi_array, $gi);
            
                my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_array);
                
                my $make_bisi = TRUE;
                foreach my $binding_site (@non_overlapping_bisis) {
                    
                    my $overlap = $gi_set->gi_set_overlap($binding_site->genomic_interval_set);
                    if ($overlap) {
                        $make_bisi = FALSE;
                        #$GU->user_info(3,$gi->{five_prime_pos}."\t".$gi->{three_prime_pos}."\t".$gi->{gi_sequence}."\tOVERLAPS with ".
                        #${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{five_prime_pos}."\t".${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{three_prime_pos}.
                        #"\t".${$binding_site->genomic_interval_set->genomic_interval_set}[0]->{gi_sequence}."\n");
                    }
                    
                }
                
                if ($make_bisi) {
                    
                    my $bisi = BiSi->new(genomic_interval_set=>$gi_set, single_site_pvalue=>$ordered_sites[$i]->{PVALUE});
                    
                    push (@non_overlapping_bisis, $bisi);
                    #$GU->user_info(3,$gi->{five_prime_pos}."\t".$gi->{three_prime_pos}."\t".$gi->{gi_sequence}."\n");
                }
            }              
        }

        # returns array of n best nonoverlapping sites, as an array of BiSis
        my $bisi_set = BiSi_Set->new(bisi_set => \@non_overlapping_bisis);
        return $bisi_set;

    } # get_n_nonoverlapping_sites #
	
} # WM_Scores #

class Pair_WM_Scores {	
	use Data::Dumper;
	use General_Utilities;
	use BiSi_Set;
	use constant {FALSE => 0,
		      TRUE  => 1};	

	has 'raw_result' => (is => 'ro', isa => 'ArrayRef', required => 1);# trigger => \&create_processed_result);
	has 'gi_set' => (is => 'ro', isa => 'Genomic_Interval_Set', required => 1);
	has 'processed_result' => (is => 'rw', isa => 'ArrayRef');

	my $GU = General_Utilities->new();
	
	method process_raw_result() {		
		my @ordered_sites_one_result;
		my @sites_all_results;
		my $gi_counter = 0;
		
		my @all_raw_scores = @{$self->raw_result};
		my @all_raw_scores_combined;
		
		# Put all scores for this gi set in a single array
		for my $i (0 .. @all_raw_scores-1) {
			
			if (exists($all_raw_scores[$i])) {
				
				# Get all scores for gi with index i in set
				my $gi_scores = $all_raw_scores[$i];
				#$GU->user_info(3,"->".$scores);
				#$GU->user_info(3,"\n");
				my $scores = $$gi_scores[0];
				#$GU->user_info(3,"\t-->".$s."\n");
				
				foreach my $s (@{$scores}){
						
					push (@all_raw_scores_combined, @{$s});
				}
				
			}
			else {
				#$GU->user_info( 3, "no scores for this sequence!\n" );
				die "No Scores for this sequence!\n";
			}
			
			$gi_counter++;
			
		}
		
		# Sort array by score
		#$GU->user_info(3,$all_raw_scores_combined[0]{POSITION_IN_GI});
		my @sorted_all_scores = (sort { $a->{SCORE} <=> $b->{SCORE} } @all_raw_scores_combined); # smallest->biggest
		
		$self->{processed_result} = \@sorted_all_scores;
	} # process_raw_result #
	
	method get_n_nonoverlapping_sites(Int $n_sites){
		
		my @ordered_sites = @{$self->processed_result};
		
		my @non_overlapping_bisis;
		my $gi_set = $self->gi_set;
		my $five_pos_centre;
		my $three_pos_centre;
		my $five_pos_flank;
		my $three_pos_flank;
		my $STRAND;
		# make best site into gi->bisi, push into array.
		my $best_site = $ordered_sites[0];
			
		# Define the genomic interval set for this, top scoring pair.

		my $genomic_interval = ${$self->gi_set->genomic_interval_set}[ $$best_site{GI_INDEX} ];
			
		#$GU->user_info(3,ref($best_site));	
		# Define gi for the central WM
		#$GU->user_info(3,Dumper($best_site));
		$GU->user_info(3,"-->".$$best_site{CENTRAL_WM_INDEX}."\n");
		#$GU->user_info(3,$genomic_interval->five_prime_pos."\n");
		#if ($$best_site{STRAND} eq 'positive') {
		#	$five_pos_centre = $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX};
		#	$three_pos_centre =  $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX}+$$best_site{CENTRAL_WM_LEN} -1;
		#}
		#elsif ($$best_site{STRAND} eq 'negative') {
		#	$five_pos_centre = $genomic_interval->five_prime_pos-$$best_site{CENTRAL_WM_INDEX} - 1;
		#	$three_pos_centre =  $genomic_interval->five_prime_pos-($$best_site{CENTRAL_WM_INDEX}+$$best_site{CENTRAL_WM_LEN});
		#}
		###############
		if ($$best_site{STRAND} eq 'positive') {
			#$five_pos = $genomic_interval->five_prime_pos+$best_site->{SITE};
			#$three_pos =  $genomic_interval->five_prime_pos+$best_site->{SITE}+$length;
			
			if($$best_site{SOURCESTRAND} eq 'positive'){
				#$GU->user_info(3,"MARKER1\n");
				$five_pos_centre = $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX};
				$three_pos_centre =  $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX}+$$best_site{CENTRAL_WM_LEN}-1;
				$STRAND = 'positive';
			}
			elsif($$best_site{SOURCESTRAND} eq 'negative'){
				#$GU->user_info(3,"MARKER2\n");
				$five_pos_centre = $genomic_interval->five_prime_pos-$$best_site{CENTRAL_WM_INDEX};
				$three_pos_centre =  $genomic_interval->five_prime_pos-$$best_site{CENTRAL_WM_INDEX} - $$best_site{CENTRAL_WM_LEN} +1;
				$STRAND = 'negative';
			}
			else{
				die;
			}
			
		}
		elsif ($$best_site{STRAND} eq 'negative') {
			
			#$three_pos = $genomic_interval->five_prime_pos-$best_site->{SITE};
			#$five_pos =  $genomic_interval->five_prime_pos-$best_site->{SITE}+$length;
			
			if($$best_site{SOURCESTRAND} eq 'positive'){
				#$GU->user_info(3,"MARKER3\n");
				#$five_pos = $genomic_interval->five_prime_pos+$best_site->{SITE};
				#$three_pos =  $genomic_interval->five_prime_pos+$best_site->{SITE}+$length;
				#$STRAND = 'positive'; The above 3 lines get the gi for the revcom of hit on the positive strand
				$three_pos_centre = $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX};
				$five_pos_centre =  $genomic_interval->five_prime_pos+$$best_site{CENTRAL_WM_INDEX} + $$best_site{CENTRAL_WM_LEN} -1;
				$STRAND = 'negative';
			}
			elsif($$best_site{SOURCESTRAND} eq 'negative'){
				#$GU->user_info(3,"MARKER4\n");
				$three_pos_centre = $genomic_interval->five_prime_pos-$$best_site{CENTRAL_WM_INDEX};
				$five_pos_centre =  $genomic_interval->five_prime_pos-$$best_site{CENTRAL_WM_INDEX} - $$best_site{CENTRAL_WM_LEN} +1;
				$STRAND = 'positive';
				#$GU->user_info(3,"5' = $five_pos, 3' = $three_pos\n");
			}
			else{
				die;
			}
			
		}
		#############
			
		my $gi_centre = Genomic_Interval->new(
							
							genome_db => $genomic_interval->genome_db,
							coord_sys_name => $genomic_interval->coord_sys_name,
							region => $genomic_interval->region,
							five_prime_pos => $five_pos_centre,
							three_prime_pos => $three_pos_centre,
							strand => $STRAND,#$$best_site{STRAND},
							label => $genomic_interval->label
		);

		$gi_centre->get_sequence();	
		# Define gi for the flanking WM
			
		#if ($$best_site{STRAND} eq 'positive') {
		#$five_pos_flank = $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX};
		#$three_pos_flank =  $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX}+$$best_site{FLANKING_WM_LEN} -1;
		#}
		#elsif ($$best_site{STRAND} eq 'negative') {
		#	$three_pos_flank = $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX}-1;
		#	$five_pos_flank =  $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX}+$$best_site{FLANKING_WM_LEN} -1;
		#}
		
		if ($$best_site{STRAND} eq 'positive') {
			#$five_pos = $genomic_interval->five_prime_pos+$best_site->{SITE};
			#$three_pos =  $genomic_interval->five_prime_pos+$best_site->{SITE}+$length;
			
			if($$best_site{SOURCESTRAND} eq 'positive'){
				#$GU->user_info(3,"MARKER1\n");
				$five_pos_flank = $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX};
				$three_pos_flank =  $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX}+$$best_site{FLANKING_WM_LEN}-1;
				$STRAND = 'positive';
			}
			elsif($$best_site{SOURCESTRAND} eq 'negative'){
				#$GU->user_info(3,"MARKER2\n");
				$five_pos_flank = $genomic_interval->five_prime_pos-$$best_site{FLANKING_WM_INDEX};
				$three_pos_flank =  $genomic_interval->five_prime_pos-$$best_site{FLANKING_WM_INDEX} - $$best_site{FLANKING_WM_LEN} +1;
				$STRAND = 'negative';
			}
			else{
				die;
			}
			
		}
		elsif ($$best_site{STRAND} eq 'negative') {
			
			#$three_pos = $genomic_interval->five_prime_pos-$best_site->{SITE};
			#$five_pos =  $genomic_interval->five_prime_pos-$best_site->{SITE}+$length;
			
			if($$best_site{SOURCESTRAND} eq 'positive'){
				#$GU->user_info(3,"MARKER3\n");
				#$five_pos = $genomic_interval->five_prime_pos+$best_site->{SITE};
				#$three_pos =  $genomic_interval->five_prime_pos+$best_site->{SITE}+$length;
				#$STRAND = 'positive'; The above 3 lines get the gi for the revcom of hit on the positive strand
				$three_pos_flank = $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX};
				$five_pos_flank =  $genomic_interval->five_prime_pos+$$best_site{FLANKING_WM_INDEX} + $$best_site{FLANKING_WM_LEN} -1;
				$STRAND = 'negative';
			}
			elsif($$best_site{SOURCESTRAND} eq 'negative'){
				#$GU->user_info(3,"MARKER4\n");
				$three_pos_flank = $genomic_interval->five_prime_pos-$$best_site{FLANKING_WM_INDEX};
				$five_pos_flank =  $genomic_interval->five_prime_pos-$$best_site{FLANKING_WM_INDEX} - $$best_site{FLANKING_WM_LEN} +1;
				$STRAND = 'positive';
				#$GU->user_info(3,"5' = $five_pos, 3' = $three_pos\n");
			}
			else{
				die;
			}
			
		}

		
		my $gi_flank = Genomic_Interval->new(
			
							genome_db => $genomic_interval->genome_db,
							coord_sys_name => $genomic_interval->coord_sys_name,
							region => $genomic_interval->region,
							five_prime_pos => $five_pos_flank,
							three_prime_pos => $three_pos_flank,
							strand => $STRAND,
							label => $genomic_interval->label
			);
		$gi_flank->get_sequence();	
		my @gi_set;
		push(@gi_set, $gi_centre);
		push(@gi_set, $gi_flank);
			
		my $best_hit_genomic_interval_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_set);
			
		my $best_bisi = BiSi->new(genomic_interval_set => $best_hit_genomic_interval_set, bisi_score => $$best_site{SCORE});
			
		push (@non_overlapping_bisis, $best_bisi);
		
		#$GU->user_info(3,"BEST SITE:\n");
		#$GU->user_info(3,Dumper($best_bisi));
		# ---------------------------------------------------------------------------------
		# Now add more non-overlapping BiSi's until the list is n long
		# ---------------------------------------------------------------------------------
		my $counter = 1;
			
		for (my $i = 1; $i < scalar(@ordered_sites); $i ++) {		
				
			# If the non-overlapping list of BiSi's is of length n then stop
			if ($counter < $n_sites) {
					
				# for remaining pairs of sites, turn into gi set.  If overlap, next; else add to array.
				# continue until array is length n, or all sites have been checked.
				# Define the genomic interval set for this, top scoring pair.
					
				my $next_site = $ordered_sites[$i];
				
				# ***Important*** If flanking wm index is 'NA' then this means that this is an invalid score and should be totally ignored by this procedure
				if ($$next_site{FLANKING_WM_INDEX} eq 'NA'){
					next;
				}
				
				my $genomic_interval = ${$self->gi_set->genomic_interval_set}[ $$next_site{GI_INDEX} ];
					
					
				# Define gi for the central WM
					
				if ($$next_site{STRAND} eq 'positive') {
					$five_pos_centre = $genomic_interval->five_prime_pos+$$next_site{CENTRAL_WM_INDEX};
					$three_pos_centre =  $genomic_interval->five_prime_pos+$$next_site{CENTRAL_WM_INDEX}+$$next_site{CENTRAL_WM_LEN} -1;
				}
				elsif ($$next_site{STRAND} eq 'negative') {
					$three_pos_centre = $genomic_interval->five_prime_pos+$$next_site{CENTRAL_WM_INDEX};
					$five_pos_centre =  $genomic_interval->five_prime_pos+$$next_site{CENTRAL_WM_INDEX}+$$next_site{CENTRAL_WM_LEN} -1;
				}
					
				
				my $gi_centre = Genomic_Interval->new(
					
								genome_db => $genomic_interval->genome_db,
								coord_sys_name => $genomic_interval->coord_sys_name,
								region => $genomic_interval->region,
								five_prime_pos => $five_pos_centre,
								three_prime_pos => $three_pos_centre,
								strand => $$next_site{STRAND},
								label => $genomic_interval->label
					);
					
					
					
				# Define gi for the flanking WM
					
				if ($$next_site{STRAND} eq 'positive') {
					
					$five_pos_flank = $genomic_interval->five_prime_pos+$$next_site{FLANKING_WM_INDEX};
					$three_pos_flank =  $genomic_interval->five_prime_pos+$$next_site{FLANKING_WM_INDEX}+$$next_site{FLANKING_WM_LEN} -1;
				}
				elsif ($$next_site{STRAND} eq 'negative') {
					
					$three_pos_flank = $genomic_interval->five_prime_pos+$$next_site{FLANKING_WM_INDEX};
					$five_pos_flank =  $genomic_interval->five_prime_pos+$$next_site{FLANKING_WM_INDEX} + $$next_site{FLANKING_WM_LEN} -1;
				}
				
					
				my $gi_flank = Genomic_Interval->new(
					
								genome_db => $genomic_interval->genome_db,
								coord_sys_name => $genomic_interval->coord_sys_name,
								region => $genomic_interval->region,
								five_prime_pos => $five_pos_flank,
								three_prime_pos => $three_pos_flank,
								strand => $$next_site{STRAND},
								label => $genomic_interval->label
					);
					
				my @gi_set;
				push(@gi_set, $gi_centre);
				push(@gi_set, $gi_flank);
				#$GU->user_info(3,$gi_centre->five_prime_pos."\n");	
				my $next_hit_genomic_interval_set = Genomic_Interval_Set->new(genomic_interval_set => \@gi_set);
				###########
				#Test for overlap of this gi set against each existing site in bisi array
					
				my $make_bisi = TRUE;
				
				# Make BiSi anyway but delete following two lines
				my $bisi = BiSi->new(genomic_interval_set =>$next_hit_genomic_interval_set, bisi_score => $$next_site{SCORE});
				#$GU->user_info(3,Dumper($bisi));
				foreach my $binding_site (@non_overlapping_bisis) {
						#$GU->user_info(3,"-->".$binding_site."\n");
					my $overlap = $next_hit_genomic_interval_set->gi_set_overlap($binding_site->genomic_interval_set);
						
					if ($overlap) {
						$make_bisi = FALSE;
					}
				}
				#$GU->user_info(3,"IS OVERLAP = $make_bisi\n");
				if ($make_bisi) {
					my $bisi = BiSi->new(genomic_interval_set =>$next_hit_genomic_interval_set, bisi_score => $$next_site{SCORE});
							
					push (@non_overlapping_bisis, $bisi);
					
					$counter++;
				}
				
					
					#########
				}
				
			if ($counter == $n_sites) {
				#$GU->user_info( 3, "\n\nn sites found!\n" );
				#$GU->user_info(3,"\n\nn sites found!\n");
				last;
			}
			
		}
		
		# returns array of n best nonoverlapping sites, as an array of BiSis
		my $bisi_set = BiSi_Set->new(bisi_set => \@non_overlapping_bisis);
		#$GU->user_info(3,Dumper($bisi_set));
		#$GU->user_info(3,@non_overlapping_bisis."\n");
		return $bisi_set;
	} # get_n_nonoverlapping_sites #

} # Pair_WM_Scores #

