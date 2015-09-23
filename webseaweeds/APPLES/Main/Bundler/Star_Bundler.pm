### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Star_Bundler Class ###
# bundles results of window-pair algorithms based on pairwise comparisons of n species against one "centre species"

use MooseX::Declare;

class Bundler::Star_Bundler {
	use constant {FALSE => 0,
			TRUE => 1};
    use Bundler::Dependencies::APPLES_Datatypes qw(APPLESSpeciesName);
    use Bundler::Dependencies::Partial_Threshold_Matrix;
	use Bundler::Dependencies::General_Utilities;
	use Data::Dumper;


	has 'sr_tables' => (is => 'rw', isa => 'ArrayRef');
    
    #<=== Check if score is better ===>#
    #Takes two scores, A and B, and checks if A is better than B
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

    
    #<==== Compute belief score =====>#
    #Takes in an array of numbers
    #Computes a belief score of 0 to 1 based on this array
    method private_compute_belief_score (ArrayRef $beliefscores_arrayref) {
        my (@beliefscores) = @{$beliefscores_arrayref};
        
        #Start with 1
        my ($result) = 1;

        #Multiply 1 by each score, in turn
        #If we have a lot of numbers, we will thus have a very low score
        for my $i( 0..$#beliefscores ) {
            $result *= 1-$beliefscores[$i];
        }
        
        #Now we do 1 minus the score. A low score is thus turned into a high score
        $result = 1-$result;
        
        #Return the result
        return $result;
	} # private_compute_belief_score #
    
    
    #<==== Count N's in a sequence =====>#
    method private_count_ns (Str $sequence ) {
        
        #Split sequence into characters
        my @bases = split //, $sequence;
        my $nCount = 0;
        
        #Check if each base is an N or n, and add up the count
        foreach( @bases ) {
            if( /N/i ) { # i flag to match N and n
                $nCount++;
            }
        }
        return $nCount;
	} # private_count_ns #
    
    
    #<==== Unimplemented method =====>#
    method private_belief_score() {
        ### belief-score as individual method
        die 'method not implemented yet!';
	} # private_belief_score #
    
    
    #In progress
    method private_build_one_sr_table (Int $position, APPLESSpeciesName $firstspecies, APPLESSpeciesName $secondspecies, Bundler::Dependencies::Partial_Threshold_Matrix $ptm ) {
        my @SRtables = @{$self->sr_tables};# deref the sr array
        my ($A, $B, $nextscore, $integrate) = 0;
        my ($S, $R);
        
        my @threshold = $ptm->get_thresholds($firstspecies, $secondspecies);
        $A = $threshold[0]; #min threshold
        $B = $threshold[1]; #max threshold
        
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
	
} # Star_Bundler #
