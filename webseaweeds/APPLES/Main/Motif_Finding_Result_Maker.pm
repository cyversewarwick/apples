### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Motif_Finding_Result_Maker class ###
## Modified 31/07/09: Richard Hickman

use MooseX::Declare;

class Motif_Finding_Result_Maker {
  use Parameters;
  use Motif_Finding_Result;
  use Data::Dumper;
  use Statistics;
  use General_Utilities;
  use Generic_Sequence_Pattern;

  my $GU = General_Utilities->new();
	
  # This method, given a meme.txt output file, parses all the relevent data and then places it into a
  # MEME_Motif_Finding_Result object
  method make_motif_finding_result_from_meme_text_output(Str $meme_text_output_filepath) {
		
    #$meme_text_output_filepath = '/Users/richardhickman/meme.txt';
    my $meme_data = $self->parse_MEME_text_output($meme_text_output_filepath);
    
		
    my @MEME_motif_finding_result;
	  
    # For each motif found by MEME, load info into a new motif_finding_result object
    foreach my $motif (@{$meme_data}) {
	    
      # Create wm object for this motif
      my @wm = @{$motif->{WM}};
      my $weight_matrix = Generic_Weight_Matrix->new(wm_identifier=>$motif->{MOTIF_NUM}, wm_freq_dist=>\@wm);
			
      # Calculate the closest matching WM in TRANSFAC
      
    	
      # Define a new motif_finding_result object to hold info on individual motif
     
      my $motif_finding_result = 
	  MEME_Motif_Finding_Result->new(num_sites => $motif->{SITES},
				       motif_width => $motif->{WIDTH},
				       occurrence_ratio => $motif->{RATIO},
				       e_value => $motif->{EVALUE},
				       positional_bias_pvalue => $motif->{POS_BIAS},
				       strand_bias_pvalue => $motif->{STRAND_BIAS},
				       WM => $weight_matrix,
				       which_genes => \@{$motif->{GENES}},	
				       );
					
			
      # Add this motif result object to the list of motifs found in this MEME run.
      push(@MEME_motif_finding_result, $motif_finding_result);
			
    }	
    #Retain contents of meme output
    open(MEME_FILE, "<$meme_text_output_filepath") || die("Cannot open $meme_text_output_filepath $!\n");
    my $meme_txt = do { local $/; <MEME_FILE>};
	close(MEME_FILE);	
	unshift(@MEME_motif_finding_result, $meme_txt);
    return \@MEME_motif_finding_result;		
  } # make_motif_finding_result_from_meme_text_output #
	
  method parse_MEME_text_output(Str $inputFile) {	
	
	open (FILEHANDLE, "< " . $inputFile) or die "oops: Cannot open $inputFile $!";
    my $lineNumber = 0;
    my $nmotifs = 0;
    my $numSeqs;
    my %promotor_lengths = ();
    my $longest_promoter = 0;
    my @motif_data;
    #my $cn = substr($d, -1,1);
    #$GU->user_info( 3, "\nCLUSTER $cn\n" ); # where is $d defined?
    while (<FILEHANDLE>) {
    	my $lineOfInput = $_;
	
      	$lineNumber++;
      	
      	if ($lineOfInput =~ m/^TRAINING SET/) {
      		
      		#get promotor lengths
  
      		 while (<FILEHANDLE>){
      		 	
      		 	my $lineOfInput = $_;
      		 	
      		 	if ($lineOfInput =~ m/\S+\s+[\d\.]+\s+\d+\s+/) {
      		 		
      		 		while($lineOfInput =~ m/(\S+)\s+[\d\.]+\s+(\d+)\s+/g) {
      		 		
      		 			$promotor_lengths{$1} = $2;
      		 			if($2 > $longest_promoter){
      		 				$longest_promoter = $2;
      		 			}
      		 			
      		 		}
      		 		
      		 	}
      		 	
      		 	if ($lineOfInput =~ m/^COMMAND/) {
      		 		last;
      		 	}
      		 	
      		 }
      	}
	
      	if ($lineOfInput =~ m/^model:/) {
		
			my @split = split(/\s+/,$_);
			$nmotifs = $split[4];
      	}
	
      	if ($lineOfInput =~ m/^data/) {
			
			my @split = split(/\s+/,$_);
			$numSeqs = $split[4];
			
			for my $m (1 .. $nmotifs) {
				
	 	 		my $rightCount = 0;
	  			my $leftCount = 0;
	  			my $plusStrandCount = 0;
	  			my $minusStrandCount = 0;
	  			#my $png;
				
	  			# Initialise variables for individual motif data
	  			my $motifNum;
	  			my $width;
	  			my $sites;
	  			my $evalue;
	  			my $ratio;
	  			my @matrix;
	  			my $pos_bias;
	  			my $strand_bias_pvalue;
	  			my @which_gene_ids = ();
	  
	  			while (<FILEHANDLE>) {
	    
	  				my $lineOfInput = $_;
	    
	    			if ($lineOfInput =~ m/^MOTIF/) {
	      
	      				my @split = split(/\s+/,$_);
	      
	      				$motifNum = $split[1];
	      				my $filepath = `pwd`;
	      				chomp($filepath);
	      				#$png = $filepath."/".$d."/logo".$motifNum.".png";
	      
	      				$width = $split[4];
	      				$sites = $split[7];
	      				$evalue = $split[13];
	      				$ratio = $sites / $numSeqs;
	      				$ratio = substr($ratio, 0, 5);
	      				$GU->user_info( 3, "MOTIF $motifNum\tWidth: $width\tSites: $sites\tRatio: $ratio\tE-value:$evalue\n" );
	      
	   				 }
	     
	   				 if ($lineOfInput =~ m/\s{2,}(\+|-)\s{2,}/) {  # Sascha: changed this as Richard's previous pattern ($lineOfInput =~ m/\s{10,}(\+|-)/))
                                                          #         did not work on fungal data where gene-IDs are too long to leave 10 white-space
                                                          #         characters
						$GU->user_info( 3, "****$_****" );
	      				my @split = split(/\s+/,$_);
	      				my $strand = $split[1];
	     				my $position = $split[2];
	     				my $gene_id = $split[0];
	     				$GU->user_info( 3, $position."\n" );
	     				my $scaled_position = 0;
	      				# check if numeric
	      				unless($position =~ m/[^0-9.]/) {
	      					#scale position as each promotor a different length
							$scaled_position = ($position/($promotor_lengths{$gene_id}-$width+1)) * ($longest_promoter-$width+1);
				#			if ($position >= ($promotor_lengths{$gene_id}/2)) {
		  		#				$rightCount++;
				#			} else {
		  		#				$leftCount++;
				#			}	    
						}
	      				# update strand counts
	      				if ($strand eq '+') {
		  					$plusStrandCount = $plusStrandCount + 1;
	      				} else {
		  					if ($strand eq '-') {
		      					$minusStrandCount = $minusStrandCount + 1;	  
		  					} else {
		      					$GU->user_info( 2,"Strand was neither plus nor minus! (it was: ".$strand.")\n");
		  					}
	      				}
	      				push(@which_gene_ids, {'id', $gene_id, 'strand', $strand, 'position', $position, 'scaled_position', $scaled_position, 'promoter_length', $promotor_lengths{$gene_id} });
	    			}
					
	    			if ($lineOfInput =~ m/Time/) {
	      				# Quick! Work out the binomial probabilty of getting positional distribution
	      				#my $p;
						my $statistics = Statistics->new();
						
						my @positions = ();
						for(my $g=0;$g<=$#which_gene_ids;$g++){
							unless($which_gene_ids[$g]{'position'} =~ m/[^0-9.]/){
								push(@positions, $which_gene_ids[$g]{'scaled_position'});
							}
						}
						#Use K-S test against a linear (random) distribution
						my @cdf = (1 .. ($longest_promoter-$width+1));	#theoretical linear distribution range
						$pos_bias = $statistics->one_sided_ks_test(\@positions, \@cdf, 'linear');
						
				#		$GU->user_info( 3, "Right count: ".$rightCount."\n" );
				#		$GU->user_info( 3, "Left count: ".$leftCount."\n" );
				#		if ($rightCount >= $leftCount) {
		    	#			$p = $statistics->get_binomial_probability($rightCount+$leftCount,$rightCount,0.5);
				#		} else {
		    	#			$p = $statistics->get_binomial_probability($rightCount+$leftCount,$leftCount,0.5);
				#		}
		
						# Compute strand bias p-value
						$GU->user_info( 3, "Positive strand count: ".$plusStrandCount."\n");
						$GU->user_info( 3, "Negative strand count: ".$minusStrandCount."\n");
						if ($plusStrandCount >= $minusStrandCount) {
		    				$strand_bias_pvalue = $statistics->get_binomial_probability($plusStrandCount+$minusStrandCount,$plusStrandCount,0.5);
						} else {
		    				$strand_bias_pvalue = $statistics->get_binomial_probability($plusStrandCount+$minusStrandCount,$minusStrandCount,0.5);
						}
	     	      						
						$GU->user_info( 3, "Prob. of distribution: ".$pos_bias."\n" );
						$GU->user_info( 3, "Prob. of strand bias: ".$strand_bias_pvalue."\n" );
						#$GU->user_info( 3, "$png\n" );
	     	 			last;
						
	    			}
					
	    			if ($lineOfInput =~ m/letter-probability matrix/) {
							
	      				$GU->user_info( 3, "Hit Matrix Data\n" );
	      				while (<FILEHANDLE>) {
							chomp;
							
							if ($_ =~ m/^-/) {
			  					last;	
							}
								
							my @bases = split(/\s+/, $_);
							shift(@bases);
							foreach (@bases) {
		  						$GU->user_info( 3, "...$_...\t" );	
							}
							$GU->user_info( 3, "\n" );
								
							push(@matrix, \@bases);
								
	      				}	
	    			}
				} #end while
				
				#Sort gene ids
				my @sorted_gene_ids = sort {$a->{'id'} cmp $b->{'id'}} @which_gene_ids;
				
	  			# Store all motif info in anonomous hash
	  			my $motif_record = {
			      	MOTIF_NUM   => $motifNum,
			      	WIDTH  => $width,
			      	SITES  => $sites,
			      	RATIO  => $ratio,
			      	EVALUE => $evalue,
			      	POS_BIAS => $pos_bias,
			      	WM => \@matrix,
			      	STRAND_BIAS => $strand_bias_pvalue,
			      	GENES => \@sorted_gene_ids,
				};
	  			$GU->user_info( 3, "$motifNum\n" );
	  			push(@motif_data, $motif_record);
			}	#end for 
		} # end if	
	} # end while
	return \@motif_data;
} # parse_MEME_text_output #

} # Motif_Finding_Result_Maker #

