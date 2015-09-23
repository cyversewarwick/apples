### (c) copyright University of Warwick 2013 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### WPBundler Class ###
#Notes

use MooseX::Declare;

class WPTools::WPBundler
{
    use Data::Dumper;
    use WPTools::Conserved_Region;
    use WPTools::Region_Match;
    use List::Util qw(first max maxstr min minstr reduce shuffle sum);

    
    has 'cur_ID' => (is => 'rw', isa => 'Int');
    #has 'conserved_regions' => (is => 'rw', isa => 'HashRef');
    
    #<== add_to_conserved_regions_from_seaweed_results ==>#
    #This method takes in:
    #1. A result object from the Seaweed_Job class when run (Datatypes::Results::Alignment_Plot)
    #2. A conservation threshold. This is the score at which you judge a window pair to be significantly conserved
    #And it outputs:
    #1. A set of conserved regions in species 1 and species 2, and the mappings betweem them
    #i.e. which regions of sequence 1 are conserved in sequence 2, and where are they conserved? (and vice-versa)
    
    #<== Start method ==>#
    method add_to_conserved_regions_from_seaweed_results(Datatypes::Results::Alignment_Plot $conservation_result, Num $conservation_threshold, Str $species_1_gene_acc, Str $species_2_gene_acc) {
        
        my $start_id = $self->cur_ID();
        my $end_id;
        
        #<=== Get data ===>#
        #Get the two raw sequences
        my $seq_1 = $conservation_result->sequence_1->seq || die "\nFatal error: No sequence 1 defined";
        my $seq_2 = $conservation_result->sequence_2->seq || die "\nFatal error: No sequence 2 defined";
        
        #Get the two species we're working with
        my $species_1 = $conservation_result->sequence_1->species;
        my $species_2 = $conservation_result->sequence_2->species;

        #Get the five prime positions
        my $fiveprime_1 = $conservation_result->sequence_1->{"five_prime_pos"};
        my $fiveprime_2 = $conservation_result->sequence_2->{"five_prime_pos"};
        
        if($species_1 =~ m/unknown species/ || $species_2 =~ m/unknown species/)
        {
            die "\nWPBundler error: one or both of the species in the alignment is unknown";
        }
        
        #Get the window size
        my $window_size = $conservation_result->job->windowsize;
        
        #Retrieve the conservation scores for each window pair
        my %plot = %{$conservation_result->plot};
        
        #We keep track of which windows have high conservation scores in each sequence
        #And which window(s) they correspond to in the other species
        my %species_1_conserved_windows = ();
        my %species_2_conserved_windows = ();
        
        #<=== Find conserved regions ===>#
        #Go through all of the window pairs
        foreach my $window_pair (keys %plot)
        {
            #If a window pair has a good enough conservation score, as defined by the threshold
            if($plot{$window_pair} >= $conservation_threshold)
            {
                #Each window pair is split by "_" in the Alignment_Plot data structure
                my @window_split = split(/_/, $window_pair);
                
                #Now, add this conserved window pair to our trackers
                #Each window here is an array ref to all linking windows in the other species
                #So the species 1 window is linked to the species 2 window...
                push(@{$species_1_conserved_windows{$window_split[0]}}, $window_split[1]);
                #And the species 2 window is linked to the species 1 window...
                push(@{$species_2_conserved_windows{$window_split[1]}}, $window_split[0]);
            }
        }
        
        #Here we find the contiguous conserved regions in both the sequences
        my $seq_1_regions_ref = $self->find_contiguous_conserved_regions(\%species_1_conserved_windows, length($seq_1), $species_1, $species_1_gene_acc, $fiveprime_1, $window_size);
        my $seq_2_regions_ref = $self->find_contiguous_conserved_regions(\%species_2_conserved_windows, length($seq_2), $species_2, $species_2_gene_acc, $fiveprime_2, $window_size);
        
        #What is the last ID we added?
        $end_id = $self->cur_ID() - 1;
        
        #Now we've found the conserved regions, we want to check which parts of which regions match which parts of which regions in each species
        $self->match_contiguous_conserved_regions(\%species_1_conserved_windows, \%species_2_conserved_windows, $seq_1_regions_ref, $seq_2_regions_ref, $window_size);
        
        #Merge each region into the bundler
        #And add the sequences
        foreach my $region (@{$seq_1_regions_ref})
        {
            $region->seq($seq_1);
            $self->add_conserved_region_into_bundler($region);
        }
        foreach my $region (@{$seq_2_regions_ref})
        {
            $region->seq($seq_2);
            $self->add_conserved_region_into_bundler($region);
        }
        
	} # add #
    
    #<== find_contiguous_conserved_regions ==>#
    #Given a list of conserved windows, a sequence length, and a bunch of metadata,
    #This method will find contiguous conserved regions in a sequence and return them as an arrayref of WPTools::Conserved_Region
    
    #<== Start method ==>#
    method find_contiguous_conserved_regions (HashRef $conserved_windows, Int $sequence_length, Str $species, Str $gene_acc, Int $fiveprime, Int $window_size) {
        
        #If we find any conserved regions, we add them to this array
        my @regions = ();
        
        #Now, start searching
        my $start = -1;
        my $end = -1;
        
        #Go through sequence
        for(my $i=0;$i<$sequence_length+1;$i++)
        {
            #If this window is conserved above our threshold
            if(defined($conserved_windows->{$i}))
            {
                #If start is -1, then we aren't currently going through a conserved region
                if($start == -1)
                {
                    $start = $i;
                    $end = $i;
                }
                #Else we're still in a conserved region so we need to push the end of it forward
                else
                {
                    $end = $i;
                }
            }
            #If this window isn't conserved (or we're at length(seq1))
            else
            {
                #If we were in a conserved region - we are no longer in a conserved region so we must end it
                if($start != -1)
                {
                    #Make a new conserved region with the correct parameters
                    my $current_region = WPTools::Conserved_Region->new;
                    
                    #Set the meta information
                    $current_region->ID($self->cur_ID);
                    $current_region->species($species);
                    $current_region->gene_accession($gene_acc);
                    $current_region->five_prime_pos($fiveprime);

                    #Set the region information
                    $current_region->start($start);
                    $current_region->end($end + ($window_size-1));
                    
                    #Add this to the output array
                    push(@regions, $current_region);
                    
                    #Move up an ID
                    $self->cur_ID($self->cur_ID + 1);
                    
                    #Start looking for new regions again
                    $start = -1;
                }
            }
        }
        #Here we go!
        return \@regions;
    }
    
    method match_contiguous_conserved_regions (HashRef $conserved_windows_1, HashRef $conserved_windows_2, ArrayRef $conserved_regions_1, ArrayRef $conserved_regions_2, Int $window_size)    {
        
        #Translate location mappings to conserved region mapping
        my %region_map_1 = ();
        my %region_map_2 = ();
        
        #For each window map, we find which region it belongs to, and which it maps to
        #First for sequence 1 window mappings
        foreach my $seq_1_window (keys %{$conserved_windows_1})
        {            
            #We match the sequence 2 window(s) to the correct region(s)
            #So we check each region
            foreach my $region_2 (@{$conserved_regions_2})
            {
                #To see which seq 2 window it maps to
                foreach my $seq_2_window (@{$conserved_windows_1->{$seq_1_window}})
                {
                    #If we find a match
                    if($seq_2_window >= $region_2->start && $seq_2_window <= $region_2->end - ($window_size - 1))
                    {
                        #Then add to the translation map
                        push(@{$region_map_1{$seq_1_window}}, $region_2->ID);
                    }
                }
            }
        }
        
        #Do the same for sequence 2 mapping back to sequence 1
        foreach my $seq_2_window (keys %{$conserved_windows_2})
        {
            foreach my $region_1 (@{$conserved_regions_1})
            {
                foreach my $seq_1_window (@{$conserved_windows_2->{$seq_2_window}})
                {
                    if($seq_1_window >= $region_1->start && $seq_1_window <= $region_1->end - ($window_size - 1))
                    {
                        push(@{$region_map_2{$seq_2_window}}, $region_1->ID);
                    }
                }
            }
        }
        
        my %cur_regions = ();
        
        #Now we have the matches for each window pair, we simply need to assemble them into contiguous sections
        #And match them with the appropriate regions
        foreach my $seq_1_region (@{$conserved_regions_1})
        {
            #Go through each region
            for(my $i=$seq_1_region->start;$i<$seq_1_region->end - ($window_size - 3);$i++)
            {
                #If we have a map, then start building a contiguous sequence
                if(defined($region_map_1{$i}))
                {
                    foreach my $seq_2_region (@{$region_map_1{$i}})
                    {
                        #Cur regions is the current regions we're building, so fix the start and end as appropriate if we're still in it
                        if(defined($cur_regions{$seq_2_region}))
                        {
                            $cur_regions{$seq_2_region}->{"end"} = $i;
                        }
                        else
                        {
                            $cur_regions{$seq_2_region}->{"start"} = $i;
                            $cur_regions{$seq_2_region}->{"end"} = $i;
                        }
                    }
                }
                
                #Add the completed contiguous matches to the region itself
                foreach my $cur_region (keys %cur_regions)
                {
                    if($cur_regions{$cur_region}->{"end"} < $i)
                    {
                        my $current_match = WPTools::Region_Match->new;
                        
                        #Set the match region information
                        $current_match->ID_match($cur_region);
                        $current_match->start($cur_regions{$cur_region}->{"start"});
                        $current_match->end($cur_regions{$cur_region}->{"end"} + ($window_size - 1));
                        push(@{$seq_1_region->{"other_species_matches"}}, $current_match);
                        delete $cur_regions{$cur_region};
                    }
                }
            }
        }
        
        #Do the same for seq 2
        foreach my $seq_2_region (@{$conserved_regions_2})
        {
            for(my $i=$seq_2_region->start;$i<$seq_2_region->end - ($window_size - 3);$i++)
            {
                if(defined($region_map_2{$i}))
                {
                    foreach my $seq_1_region (@{$region_map_2{$i}})
                    {
                        if(defined($cur_regions{$seq_1_region}))
                        {
                            $cur_regions{$seq_1_region}->{"end"} = $i;
                        }
                        else
                        {
                            $cur_regions{$seq_1_region}->{"start"} = $i;
                            $cur_regions{$seq_1_region}->{"end"} = $i;
                        }
                    }
                }
                
                foreach my $cur_region (keys %cur_regions)
                {
                    if($cur_regions{$cur_region}->{"end"} < $i)
                    {
                        my $current_match = WPTools::Region_Match->new;
                        
                        $current_match->ID_match($cur_region);
                        $current_match->start($cur_regions{$cur_region}->{"start"});
                        $current_match->end($cur_regions{$cur_region}->{"end"} + ($window_size - 1));
                        
                        push(@{$seq_2_region->{"other_species_matches"}}, $current_match);
                        delete $cur_regions{$cur_region};
                    }
                }
            }
        }
    }
    
    #<== merge_conserved_region_into_bundler ==>#
    #Given a conserved region, this method will merge it into the bundler
    method add_conserved_region_into_bundler (WPTools::Conserved_Region $new_region)   {
        $self->{"conserved_regions"}->{$new_region->{"ID"}} = $new_region;
    }
    
    #<== do_regions_overlap ==>
    #Returns 1 if two regions overlap, 0 if they don't
    #Note: the regions must both be with reference to the same gene
    method do_regions_overlap (WPTools::Conserved_Region $one, WPTools::Conserved_Region $two) {
        
        if($one->{"gene_accession"} ne $two->{"gene_accession"})
        {
            return 0;
        }
        
        my $as = $one->{"five_prime_pos"} + $one->{"start"};
        my $ae = $one->{"five_prime_pos"} + $one->{"end"};
        
        my $bs = $two->{"five_prime_pos"} + $two->{"start"};
        my $be = $two->{"five_prime_pos"} + $two->{"end"};
        
        #If a start overlaps b
        if($as <= $be && $as >= $bs)
        {
            return 1;
        }
        
        #If a end overlaps b
        if($ae <= $be && $ae >= $bs)
        {
            return 1;
        }
        
        #If b start overlaps a
        if($bs <= $ae && $bs >= $as)
        {
            return 1;
        }
        
        #If b end overlaps a
        if($be <= $ae && $be >= $as)
        {
            return 1;
        }
        
        return 0;
    }
    
    #<== merge_regions ==>
    #Merges two overlapping regions
    #Note: they MUST overlap, so check first
    method merge_regions    (WPTools::Conserved_Region $one, WPTools::Conserved_Region $two, Int $newID) {
        
        #Create a new region
        my $new_region = WPTools::Conserved_Region->new;
        $new_region->ID($newID);
        $new_region->species($one->species);
        $new_region->gene_accession($one->gene_accession);
        
        

        #We're going to take the leftmost fiveprime region
        #Make a table to translate co-ordinates
        my $translation_table = {"1" => 0, "2" => 0};
        if($one->{"five_prime_pos"} >= $two->{"five_prime_pos"})
        {
            $new_region->five_prime_pos($two->{"five_prime_pos"});
            $translation_table->{"1"} = $one->{"five_prime_pos"} - $two->{"five_prime_pos"};
        }
        else
        {
            $new_region->five_prime_pos($one->{"five_prime_pos"});
            $translation_table->{"2"} = $two->{"five_prime_pos"} - $one->{"five_prime_pos"};
        }
        
        #Now we want to create the actual region
        #Get real starts and ends
        my $as = ($one->five_prime_pos + $one->start);
        my $ae = ($one->five_prime_pos + $one->end);
        
        my $bs = ($two->five_prime_pos + $two->start);
        my $be = ($two->five_prime_pos + $two->end);

        
        #Which starts and ends are better?
        my $new_start = min($as, $bs) - $new_region->five_prime_pos;
        my $new_end = max($ae, $be) - $new_region->five_prime_pos;
        
        $new_region->start($new_start);
        $new_region->end($new_end);
        
        
        #Concatenate the strings correctly
        my @new_seq;
        
        #Where does seq 1 start in the new seq?
        #Just start of seq 1, minus start of new seq to get the translation factor
        my $seq_1_start_in_new_seq = ($one->start + $one->five_prime_pos) - ($new_start + $new_region->five_prime_pos);
        my $actual_seq_1 = $self->get_conserved_region_sequence($one);
        
        for(my $i=0;$i<length($actual_seq_1);$i++)
        {
            #Actual position in the sequence is
            $new_seq[$seq_1_start_in_new_seq + $i] = substr($actual_seq_1, $i, 1);
        }
        
        #Do same for seq2
        my $seq_2_start_in_new_seq = ($two->start + $two->five_prime_pos) - ($new_start + $new_region->five_prime_pos);
        my $actual_seq_2 = $self->get_conserved_region_sequence($two);

        for(my $i=0;$i<length($actual_seq_2);$i++)
        {
            $new_seq[$seq_2_start_in_new_seq + $i] = substr($actual_seq_2, $i, 1);
        }
        
        #Set the new sequence
        $new_region->seq(join("", @new_seq));
        
        #Now fix the species matches
        my @other_spec = (@{$one->other_species_matches}, @{$two->other_species_matches});
        foreach my $other (@other_spec)
        {
            my $new_match = WPTools::Region_Match->new;
            $new_match->new_from_region_match($other);
            push(@{$new_region->{"other_species_matches"}}, $new_match);
        }
        
        return $new_region;
    }   # merge_regions #
    
    method get_conserved_region_sequence (WPTools::Conserved_Region $reg) {
        return substr($reg->seq, $reg->start, ($reg->end - $reg->start) + 1);
    }
	
} # WPBundler #
