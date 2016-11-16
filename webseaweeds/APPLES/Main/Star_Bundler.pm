### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Star_Bundler Class ###
# bundles results of window-pair algorithms based on pairwise comparisons of n species against one "centre species"

use MooseX::Declare;

class Star_Bundler {		
	use Cons_Elem_Raw_Data;
	use ReMo_Data;	
	use Parameters;
	use Partial_Threshold_Matrix;
	use APPLES_Datatypes qw(APPLESSpeciesName);
	use constant {FALSE => 0,
			TRUE => 1};
	use General_Utilities;
	use Genomic_Interval_Set;
	use Data::Dumper;

	has 'sr_tables' => (is => 'rw', isa => 'ArrayRef'); 

	my $GU = General_Utilities->new();
    
    method truncated_bundle (Genomic_Interval $centre_gi, Genomic_Interval_Set $gi_set_for_comparison, ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, ArrayRef $alignment_results)    {
        
        #Make GI set
        # print "Star_Bundler line 28.\n";
        push (my @combined_gi_array, $centre_gi);
        # print "Star_Bundler line 30.\n";
        push (@combined_gi_array, @{$gi_set_for_comparison->{genomic_interval_set}});
        # print "Star_Bundler line 32.\n";
        my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@combined_gi_array);
        # print "Star_Bundler line 34.\n";
        #Get remo datalist
        my @remo_datalist = $self->get_truncated_datalist($parameters, $gi_set, $alignment_results);
        # print "Star_Bundler line 37.\n";
        my @single_bundles = $self->private_do_bundling(\@remo_datalist, $parameters->star_bundler_parameters, $gi_set);
        # print "Star_Bundler line 39.\n";
        my @remo_sets = $self->private_turn_bundles_into_remo_sets(\@single_bundles, \@remo_datalist, $gi_set);
        # print "Star_Bundler line 41.\n";
        return @remo_sets;
    }
    
    # method output_remo_sets (Str $filename, ArrayRef )
    #{
    #[nd]
    #}
	
	method make_remo_sets (Genomic_Interval $centre_gi, Genomic_Interval_Set $gi_set_for_comparison, ReMo_Set_Phylogenetic_Constructor_Parameters $parameters) {
	  # 1 - fetch (or create) datalist, 2 - do bundling, 3 - make remo_sets
	  push (my @combined_gi_array, $centre_gi);
	  push (@combined_gi_array, @{$gi_set_for_comparison->{genomic_interval_set}});
	  my $gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@combined_gi_array);	    
	  my @remo_datalist = $self->private_get_datalist($parameters, $gi_set);
	  my @remo_sets;
	  if (@remo_datalist) {
	    my @single_bundles = $self->private_do_bundling(\@remo_datalist, $parameters->star_bundler_parameters, $gi_set);
	    @remo_sets = $self->private_turn_bundles_into_remo_sets(\@single_bundles, \@remo_datalist, $gi_set);
	  }
	  return @remo_sets;
	} # make_remo_sets #
	
    
    method get_truncated_datalist (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, Genomic_Interval_Set $gi_set, ArrayRef $alignment_results) {

        my @remo_datalist;

        my @gi_set = @{$gi_set->genomic_interval_set};
        my $first_gi = $gi_set[0];
        $first_gi->get_sequence();
        
        my $todo = scalar (@gi_set) -1;
        my $counter = 1;
        while ($counter <= $todo) {
            my @alignment_result;
            my $second_gi = $gi_set[$counter];
            my @windowpair_result = @{$alignment_results->[$counter-1]};
            $counter++;
            
            $second_gi->get_sequence();
            
            #my @windowpair_result = @{$alignment_result[0]}; # does this deref work?
            #my @windowpair_result = @{$alignment_results->[0]};
           
            my $rmd = ReMo_Data_Object->new();
                my $s1 = Sequence->new();
                my $s2 = Sequence->new();
                $rmd->windowlength($parameters->window_pair_algorithm_parameters->{windowlength});
                $rmd->stepwidth1($parameters->window_pair_algorithm_parameters->{stepwidth});
                $rmd->stepwidth2($parameters->window_pair_algorithm_parameters->{stepwidth});
                
                $s1->species($first_gi->genome_db->{alias});
                $s1->unmaskedsequence($first_gi->get_sequence());
                $s1->maskedsequence($first_gi->get_repeatmasked_sequence);
                $s2->species($second_gi->genome_db->{alias});
                $s2->unmaskedsequence($second_gi->get_sequence());
                $s2->maskedsequence($second_gi->get_repeatmasked_sequence());
                $rmd->firstsequence($s1);
                $rmd->secondsequence($s2);
                $rmd->windowpairs(\@windowpair_result); # this is an arrayref[windowpairs]
                push (@remo_datalist, $rmd);
        }

        return @remo_datalist;
	} # private_get_datalist #
	
	method private_do_bundling (ArrayRef $remo_data_list_ref, Star_Bundler_Parameters $parameters, Genomic_Interval_Set $gi_set) {

	  my $overlaptolerance = $parameters->overlap_tolerance;
	  my @SRtables;
	  undef @SRtables;
	  $self->sr_tables(\@SRtables);	
	  
	  # -- PASSED ARGUMENTS
	  my ($beliefvalue) = $parameters->belief_value;
#	  my ($repeatcut) = $parameters->repeat_cut;
	  my @ReMoDataList = @{$remo_data_list_ref};
	  $GU->user_info ( 3, "belief value: $beliefvalue\n");
#	  $GU->user_info ( 3, "repeatcut: $repeatcut\n");

	  # -- DECLARE PRIVATE VARIABLES
	  my ($allspeciesdefinedforautomatic, $stepwidthsandwindowlengthsok) = (TRUE, TRUE);
	  my ($windowlength);
	  my ($stepwidth1);
	  my ($stepwidth2);
	  my ($acceptpiece) = FALSE;
	  my ($seqreg);
	  my (@singlebundles);
	  my (@repeatscores);
	  my ($onemaxscore);
	  my ($bestbeliefscore);
	  my ($beliefscore);
	  my ($centersequencelength);
	  my (@sequenceregions);
	  my (@conservedelementlinkers);
	  my ($actualcomparison, $maskedsequence, $unmaskedsequence);
	  my (@maxscores);
	  # 1 ~ BUILD SR-TABLES ------------------------#
	  
	  $GU->user_info( 2, "Building SR-tables..." );
	  for my $i( 0..$#ReMoDataList ) {
	    my ($actualcomparison) = $ReMoDataList[$i];
	    
	    my ($nextcheck) = $self->private_build_one_sr_table($i, $actualcomparison->firstsequence->species, $actualcomparison->secondsequence->species, $parameters->partial_threshold_matrix);

	    if ($nextcheck eq FALSE) {
	      $allspeciesdefinedforautomatic = FALSE;
	    }
	    if ($i == 0) {
	      $windowlength = $actualcomparison->windowlength;
	      $stepwidth1 = $actualcomparison->stepwidth1;
	      $stepwidth2 = $actualcomparison->stepwidth2;
	      $centersequencelength = length($actualcomparison->firstsequence->unmaskedsequence);
	      if($stepwidth1 != $stepwidth2) {
		$stepwidthsandwindowlengthsok = FALSE;
	      }
	    }
	    else {
	      if (($windowlength != $actualcomparison->windowlength) || ($stepwidth1 != $actualcomparison->stepwidth1) || ($stepwidth2 != $actualcomparison->stepwidth2)) {
		$stepwidthsandwindowlengthsok = FALSE;
	      }
	    }
	  }
	  $GU->user_info( 2, "DONE\n" ) unless (($allspeciesdefinedforautomatic eq FALSE) || ($stepwidthsandwindowlengthsok eq FALSE));

	  if (($allspeciesdefinedforautomatic eq TRUE) && ($stepwidthsandwindowlengthsok eq TRUE)) {
	    # 2 ~ COMPUTE ALL REPEAT SCORES -------------- #
	    $GU->user_info( 2, "Computing all repeat scores..." );
    
	    for my $i( 0..$#ReMoDataList+1 ) {
	      if($i == 0) {
		$actualcomparison = $ReMoDataList[$i];
		$maskedsequence = $actualcomparison->firstsequence->maskedsequence;
		$unmaskedsequence = $actualcomparison->firstsequence->unmaskedsequence;
	      }
	      else {
		$actualcomparison = $ReMoDataList[$i-1];
		$maskedsequence = $actualcomparison->secondsequence->maskedsequence;
		$unmaskedsequence = $actualcomparison->secondsequence->unmaskedsequence;
	      }
            
           
	      
	      for(my $j=0; $j < (length($unmaskedsequence)/$stepwidth1); $j++) {
		if((($j*$stepwidth1+$windowlength)-1) < length($unmaskedsequence)) {
            
		  my ($substringunmasked) = substr($unmaskedsequence, $j*$stepwidth1, $windowlength);
		  my ($substringmasked) = substr($maskedsequence, $j*$stepwidth1, $windowlength);
		  my ($Ncountunmasked) = $self->private_count_ns($substringunmasked);
		  my ($Ncountmasked) = $self->private_count_ns($substringmasked);
		  $repeatscores[$i][$j] = (100*($Ncountmasked-$Ncountunmasked))/$windowlength;
		}
	      }
	    }
	    $GU->user_info( 2, "DONE\n" );

        

	    # 3 ~ FIND ALL SINGLE BUNDLES ---------------- #
	    $GU->user_info( 2, "Starting to bundle..." );
	    for my $i( 0..($centersequencelength/$stepwidth1) ) {
	      undef $singlebundles[$i][$#ReMoDataList+1];
	    }
	    for my $i( 0..$#ReMoDataList ) {
	      $actualcomparison = $ReMoDataList[$i];
	      
	      for my $j( 0..$#{$actualcomparison->windowpairs} ) {
		my ($actualpair) = @{$actualcomparison->windowpairs}[$j]; # have to dereference arrayref windowpairs
		my ($index1) = $actualpair->offset1 / $stepwidth1;
		my ($index2) = $actualpair->offset2 / $stepwidth2;
		my ($repeatratio) = $repeatscores[0][$index1];
		
		if($repeatscores[$i+1][$index2] > $repeatratio) {
		  $repeatratio = $repeatscores[$i+1][$index2];
		}
		
		my ($scorepercentage) = (100 * $actualpair->score)/$windowlength;
		my @SRtables = @{$self->sr_tables};
		my ($beliefscore) = $SRtables[$i][$scorepercentage][$repeatratio][0];
		
		my ($integrate) = $SRtables[$i][$scorepercentage][$repeatratio][1];
		
		if($beliefscore > 0) {
		  my ($nextconselem);
		  if(!($singlebundles[$index1][0][0])) {
		    $nextconselem = Cons_Elem_Raw_Data->new();
		    $nextconselem->startbase($actualpair->offset1);
		    $nextconselem->endbase($actualpair->offset1 + $windowlength - 1);
		    $nextconselem->conservation(-13);
		    $nextconselem->targetstartbase(-13);
		    $nextconselem->targetendbase(-13);
		    $nextconselem->beliefscoreisdefined(FALSE);
		    $nextconselem->beliefscore(-13);
		    $nextconselem->repeatratio(-1);
		    $singlebundles[$index1][0][0] = $nextconselem;
		  }
		  $acceptpiece = TRUE;
		  for(my $k=0; $k < 2; $k++) {
		    for(my $l = $#{$singlebundles[$index1][$i+1]}; $l >= 0; $l--) {
		      my ($overlapResult) = $self->private_overlap(	$singlebundles[$index1][$i+1][$l]->startbase,
									$singlebundles[$index1][$i+1][$l]->endbase,
									$actualpair->offset2,
									$actualpair->offset2+$windowlength-1,
									0, $parameters);
		      if($overlapResult eq TRUE) {
			my ($scoreisbetterResult) = $self->private_score_is_better(	$singlebundles[$index1][$i+1][$l]->beliefscore,
											$singlebundles[$index1][$i+1][$l]->integrate,
											$beliefscore,
											$integrate);
			if($scoreisbetterResult eq TRUE) {
			  $acceptpiece = FALSE;
			}
			else {
			  if(($acceptpiece eq TRUE) && ($k == 1)) {
			    splice(@{$singlebundles[$index1][$i+1]}, $l, 1);
			  }
			}
		      }
		    }
		  }
		  if($acceptpiece eq TRUE) {
		    $nextconselem = Cons_Elem_Raw_Data->new();
		    $nextconselem->conservation($actualpair->score);
		    $nextconselem->repeatratio($repeatratio);
		    $nextconselem->startbase($actualpair->offset2);
		    $nextconselem->endbase($actualpair->offset2+$windowlength-1);
		    $nextconselem->targetstartbase($actualpair->offset1);
		    $nextconselem->targetendbase($actualpair->offset1+$windowlength-1);
		    $nextconselem->beliefscoreisdefined(TRUE);
		    $nextconselem->beliefscore($beliefscore);
		    $nextconselem->integrate($integrate);
              
		    push @{$singlebundles[$index1][$i+1]}, $nextconselem;
		  }
		}
	      }
	    }
          
          
          


          
	    $GU->user_info( 2, "DONE\n" );
	    # 4 ~ SELECT SINGLE BUNDLES ABOVE THRESHOLD -- #
	    $GU->user_info( 2, "Selecting single bundles above threshold ($beliefvalue)..." );
	    for(my $i=0; $i < ($centersequencelength/$stepwidth1); $i++) {
	      if($singlebundles[$i][0][0]) {
		for my $j( 0..$#ReMoDataList ) { 
		  $maxscores[$j] = 0;
		}
		for my $j( 1..$#{$singlebundles[$i]} ) {
		  for my $k (0 .. $#{$singlebundles[$i][$j]}) {
		    if(($singlebundles[$i][$j][$k]->beliefscore) > $maxscores[$j-1]) {
		      $maxscores[$j-1] = $singlebundles[$i][$j][$k]->beliefscore;
		    }
		  }
		}
		$bestbeliefscore = $self->private_compute_belief_score(\@maxscores);
		if($bestbeliefscore < $beliefvalue) {
		  for my $j( 0..$#{$singlebundles[$i]} ) {
		    undef $singlebundles[$i][$j];
		  }
		}
		else {
		  $singlebundles[$i][0][0]->beliefscoreisdefined(TRUE);
		  $singlebundles[$i][0][0]->beliefscore($bestbeliefscore);
		  for my $j( 1..$#{$singlebundles[$i]} ) {
		    $onemaxscore = $maxscores[$j-1];
		    for( my $l=$#{$singlebundles[$i][$j]}; $l >= 0; $l-- ) {
		      $maxscores[$j-1] = $singlebundles[$i][$j][$l]->beliefscore;
		      $beliefscore = $self->private_compute_belief_score(\@maxscores);
		      if($beliefscore < $beliefvalue) {
			splice(@{$singlebundles[$i][$j]}, $l, 1);
		      }
		    }
		    $maxscores[$j-1] = $onemaxscore;
		  }
		}
	      }			
	    }
	    $GU->user_info( 2, "DONE\n" );


	    # 5 ~ ASSEMBLE BUNDLES ----------------------- #
	    $GU->user_info( 2, "Assembling bundles..." );
	    my ($i)=0;
	    my ($foundbundle) = FALSE;
	    while(($i < ($centersequencelength/$stepwidth1)) && ($foundbundle eq FALSE)) {
	      if($singlebundles[$i][0][0]) {
		$foundbundle = TRUE; # FINDS FIRST BUNDLE IN LIST
	      }
	      else {
		$i += 1;
	      }
	    }
	    if($foundbundle eq TRUE) {
	      my ($j) = $i+1;
	      while($j < ($centersequencelength/$stepwidth1)) {
		if(!$singlebundles[$j][0][0]) {
		  $j += 1; # FINDS NEXT BUNDLE IN LIST
		}
		else {
		  my ($overlapResult) = $self->private_overlap(	$singlebundles[$i][0][0]->startbase,
								$singlebundles[$i][0][0]->endbase,
								$singlebundles[$j][0][0]->startbase,
								$singlebundles[$j][0][0]->endbase,
								$overlaptolerance, $parameters);
		  if($overlapResult eq TRUE) {
		    for my $k( 0..$#{$singlebundles[$i]} ) {
		      for my $m(0 .. $#{$singlebundles[$j][$k]}) {
			for(my $l=$#{$singlebundles[$i][$k]}; $l >= 0; $l--) {
			  $overlapResult = $self->private_overlap(	$singlebundles[$i][$k][$l]->startbase,
									$singlebundles[$i][$k][$l]->endbase,
									$singlebundles[$j][$k][$m]->startbase,
									$singlebundles[$j][$k][$m]->endbase,
									$overlaptolerance, $parameters);
			  if ($overlapResult eq TRUE) {
			    if($singlebundles[$i][$k][$l]->conservation > $singlebundles[$j][$k][$m]->conservation) {
			      $singlebundles[$j][$k][$m]->conservation($singlebundles[$i][$k][$l]->conservation);
			    }
			    
			    if($singlebundles[$j][$k][$m]->startbase > $singlebundles[$i][$k][$l]->startbase) {
			      $singlebundles[$j][$k][$m]->startbase($singlebundles[$i][$k][$l]->startbase);
			    }
			    
			    if($singlebundles[$i][$k][$l]->endbase > $singlebundles[$j][$k][$m]->endbase) {
			      $singlebundles[$j][$k][$m]->endbase($singlebundles[$i][$k][$l]->endbase);
			    }
			    
			    if($singlebundles[$j][$k][$m]->targetstartbase > $singlebundles[$i][$k][$l]->targetstartbase) {
			      $singlebundles[$j][$k][$m]->targetstartbase($singlebundles[$i][$k][$l]->targetstartbase);
			    }
			    
			    if($singlebundles[$i][$k][$l]->targetendbase > $singlebundles[$j][$k][$m]->targetendbase) {
			      $singlebundles[$j][$k][$m]->targetendbase($singlebundles[$i][$k][$l]->targetendbase);
			    }
			    
			    $singlebundles[$j][$k][$m]->beliefscoreisdefined(TRUE);
			    
			    if($singlebundles[$i][$k][$l]->beliefscore > $singlebundles[$j][$k][$m]->beliefscore) {
			      $singlebundles[$j][$k][$m]->beliefscore($singlebundles[$i][$k][$l]->beliefscore);
			    }
			    splice(@{$singlebundles[$i][$k]}, $l, 1);
			  }
			}
			push(@{$singlebundles[$i][$k]}, $singlebundles[$j][$k][$m]);
		      }
		      undef $singlebundles[$j][$k];
		    }
		  }
		  else {
		    $i = $j;
		    $j = $i+1;
		  }
		}
	      }
	    }
	    $GU->user_info( 2, "DONE\n" );
	  }
        
	  # 6 ~ FIND SEQUENCE CONSERVATION ------------- #
	  $GU->user_info( 2, "Finding conservation of bundles..." );
	  for my $j( 0..$#singlebundles ) {
	    if( $singlebundles[$j][0][0] ) {
	      $singlebundles[$j][0][0]->targetstartbase(-13);
	      $singlebundles[$j][0][0]->targetendbase(-13);
	      my ($k) = 1;
	      while( !$singlebundles[$j][$k][0] ) {
		$k += 1;
	      }
	      my ($maxconservation) = $singlebundles[$j][$k][0]->conservation;
	      for my $m( 1..$#{$singlebundles[$j][$k]} ) {
		if( $singlebundles[$j][$k][$m]->conservation > $maxconservation ) {
		  $maxconservation = $singlebundles[$j][$k][$m]->conservation;
		}
	      }
	      $singlebundles[$j][0][0]->conservation($maxconservation);
	    }
	  }
	  $GU->user_info( 2, "DONE\n" );
	  

	  # 7 ~ FIND REPEAT RATIO ------------- #
	  $GU->user_info( 2, "Finding repeat ratios of bundles..." );
	  for my $i( 0..$#ReMoDataList ) {
	    my ($sequencelength);
	    my ($allofsequence);
	    if( $i == 0 ) {
	      $actualcomparison = $ReMoDataList[0];
	      $sequencelength = length($actualcomparison->firstsequence->unmaskedsequence);
	    }
	    else {
	      $actualcomparison = $ReMoDataList[$i-1];
	      $sequencelength = length($actualcomparison->secondsequence->unmaskedsequence);
	    }
	    for my $j( 0..$#singlebundles ) {
	      if( $singlebundles[$j][0][0] ) {
		for my $k ( 0..$#{$singlebundles[$j][$i]} ) {
		  if ( $i == 0 ) {
		    $allofsequence = $actualcomparison->firstsequence;
		  }
		  else {
		    $allofsequence = $actualcomparison->secondsequence;
		  }
		  my ($conselseqmasked) = substr($allofsequence->maskedsequence, $singlebundles[$j][$i][$k]->startbase, 1+($singlebundles[$j][$i][$k]->endbase-$singlebundles[$j][$i][$k]->startbase));
		  my ($conselsequnmasked) = substr($allofsequence->unmaskedsequence, $singlebundles[$j][$i][$k]->startbase, 1+($singlebundles[$j][$i][$k]->endbase-$singlebundles[$j][$i][$k]->startbase));
		  my ($conselseqlength) = length($conselsequnmasked);
		  my ($Ncountunmasked) = $self->private_count_ns($conselsequnmasked);
		  my ($Ncountmasked) = $self->private_count_ns($conselseqmasked);
		  my $rr = (100*($Ncountmasked-$Ncountunmasked))/$conselseqlength;
		  $singlebundles[$j][$i][$k]->repeatratio($rr);
		}
	      }
	    }
	  }
	  $GU->user_info( 2, "DONE\n" );
        
	  return @singlebundles;
	  
	} # private_do_bundling #
	
	method private_print_bundles_to_file (ArrayRef $single_bundles_ref, Str $filename, ArrayRef $remo_datalist_ref) {
	  print "PRINTING BUNDLES\n";
	  my @bundles = @{$single_bundles_ref};
	  my @ReMoDataList = @{$remo_datalist_ref};
	  open(OUT, '>', $filename) or die "Could not open ".$filename.", died";
	  my @header = qw/ID startbase endbase targetstartbase targetendbase conservation repeatratio beliefscore masked_sequence unmasked_sequence/;
	  print OUT join("\t", @header)."\n";
	  for my $i(0 .. $#bundles) {
	    if( $bundles[$i][0][0] ) {
	      my @line;
	      for my $j( 0 .. $#{$bundles[$i]} ) {
		for my $k( 0 .. $#{$bundles[$i][$j]} ) {
		  if( $bundles[$i][$j][$k] ) { # IF A BUNDLE IS DEFINED IN REFERENCE SPECIES
		    my ($wholesequencemasked, $wholesequence);
		    my ($taketo) = $bundles[$i][$j][$k]->endbase;
		    my ($takefrom) = $bundles[$i][$j][$k]->startbase;
		    
		    if ( $j == 0 ) {
		      $wholesequence = $ReMoDataList[0]->firstsequence->unmaskedsequence;
		      $wholesequencemasked = $ReMoDataList[0]->firstsequence->maskedsequence;
		    }
		    else {
		      $wholesequence = $ReMoDataList[$j-1]->secondsequence->unmaskedsequence;
		      $wholesequencemasked = $ReMoDataList[$j-1]->secondsequence->maskedsequence;
		    }
		    
		    print OUT $ReMoDataList[$j]->firstsequence->species if $j == 0; #was id, not needed
		    print OUT $ReMoDataList[$j-1]->secondsequence->species if $j != 0;#was id, not needed
		    print OUT "	";
		    print OUT $bundles[$i][$j][$k]->startbase."	";
		    print OUT $bundles[$i][$j][$k]->endbase."	";
		    print OUT $bundles[$i][$j][$k]->targetstartbase if $j != 0;
		    print OUT "	";
		    print OUT $bundles[$i][$j][$k]->targetendbase if $j != 0;
		    print OUT "	";
		    print OUT $bundles[$i][$j][$k]->conservation."	";
		    print OUT $bundles[$i][$j][$k]->repeatratio."	";
		    print OUT $bundles[$i][$j][$k]->beliefscore."	";
		    print OUT substr($wholesequencemasked, $takefrom, ($taketo-$takefrom)+1)."	";
		    print OUT substr($wholesequence, $takefrom, ($taketo-$takefrom)+1)."\n";
		  }
		}
	      }
	    }
	  }
	  close(OUT);
	} # private_print_bundles_to_file #
	
	method private_overlap (Int $startA, Int $endA, Int $startB, Int $endB, Int $mindistance, Star_Bundler_Parameters $parameters) {
		
	  my ($result) = FALSE;
	  
	  if(	( $endA+$mindistance >= $startB ) &&
		( $startA <= $endB+$mindistance )) {
	    $result = TRUE;
	  }
	  if($mindistance == $parameters->overlap_tolerance ) {
	  }
	  return $result;
	} # private_overlap #
	
	method private_score_is_better (Num $beliefscoreA, Num $integrateA, Num $beliefscoreB, Num $integrateB) {
	  
	  my ($result) = FALSE;
	  
	  if ( $beliefscoreA > $beliefscoreB ) {
	    $result = TRUE;
	  }
	  if (( $beliefscoreA == $beliefscoreB ) && ( $integrateA >= $integrateB )) {
	    $result = TRUE;
	  }
	  return $result;
	} # private_score_is_better #
	
	method private_compute_belief_score (ArrayRef $beliefscores_arrayref) {
	  my (@beliefscores) = @{$beliefscores_arrayref};
	  
	  my ($result) = 1;
	  
	  for my $i( 0..$#beliefscores ) {
	    $result *= 1-$beliefscores[$i];
	  }
	  
	  $result = 1-$result;
	  return $result;
	} # private_compute_belief_score #
	
	method private_count_ns (Str $sequence ) {
	  my @bases = split //, $sequence;
	  my $nCount = 0;
	  
	  foreach( @bases ) {
	    if( /N/i ) { # i flag to match N and n
	      $nCount++;
	    }
	  }
	  return $nCount;
	} # private_count_ns #
	
	method private_build_one_sr_table (Int $position, APPLESSpeciesName $firstspecies, APPLESSpeciesName $secondspecies, Partial_Threshold_Matrix $ptm ) {
	  my @SRtables = @{$self->sr_tables};# deref the sr array		
	  my ($A, $B, $nextscore, $integrate) = 0;
	  my ($S, $R);
	  
        #[nd] Here is where we get the thresholds
        #my @threshold = $ptm->get_thresholds($firstspecies, $secondspecies);
        #$A = $threshold[0]; #min threshold
        #$B = $threshold[1]; #max threshold
        
        if($secondspecies eq "apis mellifera")
        {
            $A = 80;
            $B = 94;
        }
        elsif($secondspecies eq "atta cephalotes")
        {
            $A = 80;
            $B = 94;
        }
        elsif($secondspecies eq "bombyx mori")
        {
            $A = 80;
            $B = 94;
        }
        elsif($secondspecies eq "homo sapiens")
        {
            $A = 80;
            $B = 94;
        }

	# LB added pairwise thresholds for plant species, from here.
        # !currently assuming here that first species is rice! (May) need to change these hardwired values if using different first species!, e.g.grape-banana 
	
	elsif ($secondspecies eq "musa acuminata"){
	    $A = 78;
	    $B = 100;
	}
	elsif ($secondspecies eq "vitis vinifera") {
	    $A = 78;
	    $B = 100;
	}
	elsif ($secondspecies eq "arabidopsis thaliana") {
	    $A = 78; 
	    $B = 100;
	}
	elsif ($secondspecies eq "sorghum bicolor") {
	    $A = 78;
	    $B = 100;
	}
	elsif ($secondspecies eq "phoenix dactylifera") {
	    $A = 78;
	    $B = 100;
	}
	elsif ($secondspecies eq "oryza sativa") {
            $A = 78;
            $B = 100;
        }
	elsif ($secondspecies eq "populus trichocarpa") {
            $A = 78;
            $B = 100;
        }
	elsif ($secondspecies eq "medicago truncatula") {
            $A = 78;
            $B = 100;
        }
	elsif ($secondspecies eq "ArabidopsisTAIR10") {
            $A = 78;
            $B = 100;
        }
	elsif ($secondspecies eq "Niben101") {
            $A = 78;
            $B = 100;
        }
	# to here
        else
        {
            die "\nUnknown species";
        }
        
	  for my $S( 0..100 ) {
	    for my $R( 0..100 ) {
	      if( $R < 100 ) {
		$integrate = 100 / (100 - $R);
	      }
	      else {
		$integrate = 1.5;
	      }
	      if( 1.5 < $integrate ) {
		$integrate = 1.5;
	      }
	      if( $S >= $R ) {
		$integrate *= ($S - $R);
	      }
	      else {
		$integrate = 0;
	      }
	      if( $integrate > 100 ) {
		$integrate = 100;
	      }
	      	      
	      if( $integrate <= $A ) {
		$nextscore = 0; # BELIEF SCORE IS 0 BELOW LOWER THRESHOLD
	      }
	      else {
		if( $integrate >= $B ) {
		  $nextscore = 1; # BELIEF SCORE IS 1 ABOVE UPPER THRESHOLD	
		}
		else {
		  $nextscore = 1 / (1+exp(5-(10*($integrate-$A)/($B-$A)))); # BELIEF SCORE IS LOGISTIC BETWEEN THRESHOLDS
		}
	      }
	      $SRtables[$position][$S][$R][0] = $nextscore; # y-axis
	      $SRtables[$position][$S][$R][1] = $integrate; # x-axis
	    }
	  }
	  $self->sr_tables(\@SRtables);# make the changes to sr_table attribute permanent
	  return TRUE;
	} # private_build_one_sr_table #
	
    method private_turn_bundles_into_remo_sets (ArrayRef $single_bundles_ref, ArrayRef $remo_datalist_ref, Genomic_Interval_Set $gi_set) {
	  my @bundles = @{$single_bundles_ref};
	  my @ReMoDataList = @{$remo_datalist_ref};
	  my @gis = @{$gi_set->genomic_interval_set};
      print "Star_Bundler line 679.\n";
	  my @remo_sets;
	  for my $i(0 .. $#bundles) {
	    if( $bundles[$i][0][0] ) {
	      my @remos;
	      for my $j( 0 .. $#{$bundles[$i]} ) {
		my $gi = $gis[$j];
		for my $k( 0 .. $#{$bundles[$i][$j]} ) {
		  my $remo;
		  if( $bundles[$i][$j][$k] ) { # IF A BUNDLE IS DEFINED IN REFERENCE SPECIES
		    print "Star_Bundler line 689.\n";
		    my $db = $gi->genome_db;
			print "Star_Bundler line 691.\n";
		    my $five = $gi->five_prime_pos+$bundles[$i][$j][$k]->startbase;
		    my $three = $gi->five_prime_pos+$bundles[$i][$j][$k]->endbase;

		    if ($gi->strand eq 'negative') {
		      $five = $gi->five_prime_pos-$bundles[$i][$j][$k]->startbase;
		      $three = $gi->five_prime_pos-$bundles[$i][$j][$k]->endbase;
		    }
			
			my $label = 'not defined';
			if (defined $gi->label) {
				$label = $gi->label;
			}
			print "Star_Bundler line 704.\n";
		    $remo = ReMo->new('genome_db' => $db,
				      'coord_sys_name' => $gi->coord_sys_name,
				      'region' => $gi->region,
				      'five_prime_pos' => $five,
				      'three_prime_pos' => $three,
				      'strand' => $gi->strand,
				      'belief_score' => $bundles[$i][$j][$k]->beliefscore,
				      'repeat_ratio' => $bundles[$i][$j][$k]->repeatratio,
				      'conservation' => $bundles[$i][$j][$k]->conservation, 
					  'working_sequence' => $gi->working_sequence(),
					  'label' => $label
				     );
		    print "Star_Bundler line 717.\n";
		    # print Dumper($remo);
              $remo->get_sequence();# fill out sequences
            print "Star_Bundler line 720.\n";
		    push (@remos, $remo);
		  } # makes a single remo from a bundle
		}
	      }
	      # this is one remo_set
	      my $remo_set = ReMo_Set->new(remo_set => \@remos);
	      
	      push (@remo_sets, $remo_set);
	    }
	  }
	  return @remo_sets;
	} # private_turn_bundles_into_remo_sets #
	
    
    #Implemented
	method private_belief_score() {
	  ### belief-score as individual method 
	  die 'method not implemented yet!';
	} # private_belief_score #
    
    method get_conservation_profiles (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, Genomic_Interval $first_gi, Genomic_Interval $second_gi) {
        
        $GU->user_info ( 2, "getting conservation profiles\n" );
        my $job_handler = Job_Handler->new();
        my $cache_result = TRUE;
        my $cache_duration = 365; # days to keep result in cache
        $first_gi->get_working_sequence();
        $second_gi->get_working_sequence();
        my $exception_info = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler
        my $stats_required = FALSE;
        
        my @alignment_result;
        my @profile_pair;
        
        my $alignment_algorithm;
        if ($parameters->window_pair_algorithm_parameters->isa('Ott_Algorithm_Parameters') ) {
            $alignment_algorithm = "ott";
        }
        elsif ($parameters->window_pair_algorithm_parameters->isa('Seaweed_Algorithm_Parameters') ) {
            $alignment_algorithm = "seaweed";
        }
        eval {
            @alignment_result = $job_handler->handle_alignment_job($parameters, $first_gi, $second_gi, $alignment_algorithm, $cache_result, $cache_duration); # returns windowpairs, profile1, profile2
        };
        if ($@) { # deal with any errors first
            my $exception_content = $@;
            if (!UNIVERSAL::can($exception_content, 'isa')) {
                die $exception_content; # throw error string
            }
            else {
                if ($exception_content->isa('Job_Information_Exception')) {
                    $stats_required = TRUE;
                    $exception_info->merge($exception_content); #, $CPU
                }
                else {
                    die $exception_content; # throw error object
                }
            }
        }
        else { # no errors were thrown, proceed with obtaining profiles from windowpair result
            # push both profiles into an array
            my @profile1 = @{$alignment_result[1]};
            my @profile2 = @{$alignment_result[2]};
            push (@profile_pair, \@profile1, \@profile2);
        }
        if ($stats_required eq TRUE) {
            die $exception_info;
        }
        return @profile_pair;
	} # get_conservation_profiles #
    
    method private_get_datalist (ReMo_Set_Phylogenetic_Constructor_Parameters $parameters, Genomic_Interval_Set $gi_set) {
        $GU->user_info ( 2, "getting datalist\n" );
        
        my $job_handler = Job_Handler->new();
        my @remo_datalist;
        my $cache_result = TRUE;
        my $cache_duration = 365; # days to keep result in cache
        my @gi_set = @{$gi_set->genomic_interval_set};
        my $first_gi = $gi_set[0];
        $first_gi->get_sequence();
        my $exception_info = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler
        my $stats_required = FALSE;
        my $todo = scalar (@gi_set) -1;
        my $counter = 1;
        while ($counter <= $todo) {
            my @alignment_result;
            my $second_gi = $gi_set[$counter];
            $GU->user_info( 2, "COMPARING GI ".$counter." of ".$todo."\n" );
            $counter++;
            $second_gi->get_sequence();
            my $alignment_algorithm;
            if ($parameters->window_pair_algorithm_parameters->isa('Ott_Algorithm_Parameters') ) {
                $alignment_algorithm = "ott";
            }
            elsif ($parameters->window_pair_algorithm_parameters->isa('Seaweed_Algorithm_Parameters')) {
                $alignment_algorithm = "seaweed";
            }
            else {
                die 'cannot handle that type of window_pair_aligorithm_parameters yet!';
            }
            eval {
                @alignment_result = $job_handler->handle_alignment_job($parameters, $first_gi, $second_gi, $alignment_algorithm, $cache_result, $cache_duration); # also needs to return profile1, profile2
            };
            if ($@) {
                my $exception_content = $@;
                if (!UNIVERSAL::can($exception_content, 'isa')) {
                    die $exception_content; # throw error string
                }
                else {
                    if ($exception_content->isa('Job_Information_Exception')) {
                        $stats_required = TRUE;
                        $exception_info->merge($exception_content);
                        next;
                    }
                    else {
                        die $exception_content; # throw error object
                    }
                }
            } # end if ($@)
            else {
                # dereference first element of the alignment_result array - this is the array of windowpairs
                
                my @windowpair_result = @{$alignment_result[0]}; # does this deref work?
                
                my $rmd = ReMo_Data_Object->new();
                my $s1 = Sequence->new();
                my $s2 = Sequence->new();
                $rmd->windowlength($parameters->window_pair_algorithm_parameters->{windowlength});
                if ($parameters->window_pair_algorithm_parameters->isa('Ott_Algorithm_Parameters') ) {
                    $rmd->stepwidth1($parameters->window_pair_algorithm_parameters->{stepwidth1});
                    $rmd->stepwidth2($parameters->window_pair_algorithm_parameters->{stepwidth2});
                }
                elsif ($parameters->window_pair_algorithm_parameters->isa('Seaweed_Algorithm_Parameters') ) {
                    $rmd->stepwidth1($parameters->window_pair_algorithm_parameters->{stepwidth});
                    $rmd->stepwidth2($parameters->window_pair_algorithm_parameters->{stepwidth});
                }
                $s1->species($first_gi->genome_db->{alias});
                $s1->unmaskedsequence($first_gi->get_sequence());
                $s1->maskedsequence($first_gi->get_repeatmasked_sequence);
                $s2->species($second_gi->genome_db->{alias});
                $s2->unmaskedsequence($second_gi->get_sequence());
                $s2->maskedsequence($second_gi->get_repeatmasked_sequence());
                $rmd->firstsequence($s1);
                $rmd->secondsequence($s2);
                $rmd->windowpairs(\@windowpair_result); # this is an arrayref[windowpairs]
                push (@remo_datalist, $rmd);
            }
        } # end while loop	  
        if ($stats_required eq TRUE) {
            die $exception_info;
        }
        return @remo_datalist;
	} # private_get_datalist #
	
} # Star_Bundler #
