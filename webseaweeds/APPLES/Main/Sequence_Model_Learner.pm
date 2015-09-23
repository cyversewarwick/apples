### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Sequence_Model_Learner Class ###
use MooseX::Declare;

class Sequence_Model_Learner {
    use Storable qw(nstore retrieve);
    use APPLES_Datatypes qw (Boolean);
    use Parameters;
    use General_Utilities;
    use Job_Handler;
    use Genomic_Interval_Set;
    use Reg_Loc_Maker;
    use ReMo_Set_Maker;
    use Genome_DB_Utilities;
    use File::Temp qw (tempdir);
    use File::Path qw (rmtree);
    use constant {FALSE => 0,
		  TRUE  => 1};	

    # has 'model_learning_parameters' => (is => 'ro', isa => 'Partial_Sequence_Model_Learning_Params', required => 1);
    # class not fully implemented yet
    
    my $GU = General_Utilities->new();
    
    method learn_sequence_model(Partial_Sequence_Model_Learning_Params $model_learning_params, Genome_Sequence_Database_Parameters $genome_db){
	
	# Get gi_set from which to learn model
	my $sequence_class = $model_learning_params->learning_restrictions;
	my $gi_set;
	my $model;
	my $remo_set;
		
	if ($sequence_class eq '500bp_promoter') {
	    
	    $GU->user_info(3,"Getting all reg locs for genome:\n");	  
	    my $reg_loc_maker = Reg_Loc_Maker->new();
	    my @reg_locs = $reg_loc_maker->make_all_reg_locs_for_a_genome($genome_db); # this call will involve Job_Handler
	    $GU->user_info(3, "Got all reg locs for genome\n");

	    #----------------------------------------------------------------------------------------------------------------------
	    # Retrieval of all the upstream sequences for every reg loc is performed through the job handler as this does take at least an hour to run
	    #-----------------------------------------------------------------------------------------------------------------------

	    my $GDBU = Genome_DB_Utilities->new();	 
	    my @gi_set;	
	    my $core_promoter_length = 500;
	
	    my $remo_core_promoter_constructor_parameters = ReMo_Core_Promoter_Constructor_Parameters->new(stop_at_neighbouring_gene => FALSE,
													   length => $core_promoter_length);
	    my $remo_set_maker = ReMo_Set_Maker->new();
	    $remo_set = $remo_set_maker->make_remo_set_from_core_promoter(\@reg_locs, $remo_core_promoter_constructor_parameters);
	    $gi_set = Genomic_Interval_Set->new(genomic_interval_set => $remo_set->remo_set);
		
	    my $x=0;
	
	    $GU->user_info(3, "Genomic interval set to be used for learning sequences model created\n");
		
	}
	elsif ($sequence_class eq '1000bp_promoter') {
	    die "Functionality not implemented yet\n";
	}
	elsif ($sequence_class eq 'intergenic') {
	    die "Functionality not implemented yet\n";
	}
	elsif ($sequence_class eq 'coding') {
	    die "Functionality not implemented yet\n";
	}
	elsif ($sequence_class eq 'whole_genome') {
	    my $GDBU = Genome_DB_Utilities->new();
	    $gi_set = $GDBU->get_all_chromosomes($genome_db);
	}
	else{
	    die "We should not have got to here! The genomic sequence class is not of the correct type\n";	
	}
	
	# Now generate model
	if ( $model_learning_params->model->isa('Markov_Model') ){		
	    
	    my $order = $model_learning_params->model->order;
	    $model = $self->generate_markov_model($gi_set, $order);
		
	    
	}
	else {
	    die "Not a recognised model\n";
	}
	
	#$GU->user_info(3,Dumper($model));
	
	return $model;

    } # learn_sequence_model #
    
    method learn_sequence_model_through_job_handler(Partial_Sequence_Model_Learning_Params $model_learning_params, Genome_Sequence_Database_Parameters $genome_db) {
	my $function = 'learn_sequence_model';
	my $SML = Sequence_Model_Learner->new();
	my @parameters = ($model_learning_params,$genome_db);
	my $high_memory = TRUE;
	my @result = $GU->standard_call_to_job_handler($SML,$function,\@parameters,$high_memory,TRUE);
	return $result[0];
    } # learn_sequence_model_through_job_handler #

	method generate_markov_model(Genomic_Interval_Set $gi_set, Int $order){
		
		# Generate random directory in which to run markov model c script
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $tempdir = $GU->get_temp_random_directory(FALSE);

		#$GU->user_info(3,"tempdir: $tempdir\n");
		my $sequences_fp = $self->print_sequences_to_file($tempdir, $gi_set);
		my $path_to_executable = $APPLES_C_binaries."/learn_markov_modelREV493";
		my $sequence_model_fp = $tempdir."model.txt";
		$GU->user_info(3,"$path_to_executable $sequences_fp $order $sequence_model_fp\n"); 
		system( "$path_to_executable", "$sequences_fp", "$order", "$sequence_model_fp")  == 0 || die "Cannot run the markov model generation C script\n";
		
		my $model_from_file = $self->private_get_raw_model($sequence_model_fp);
		rmtree($tempdir);
		return $model_from_file;
	} # generate_markov_model #
	
	method print_sequences_to_file(Str $tempdir, Genomic_Interval_Set $gi_set){
		
		# Produce sequence input file for empirical scoring
		my $sequences_output_file = $tempdir."sequences.in";
		#$GU->user_info(3, "Printing sequence file to: $sequences_output_file\n");
		
		open SEQUENCES_OUTFILE, ">".$sequences_output_file or die "Cannot write sequences to file\n";
		#$GU->user_info(3,"->".scalar(@{$gi_set->genomic_interval_set}) ."\n");
		
		for (my $i = 0; $i < @{$gi_set->genomic_interval_set}; $i++){
			
			my $gi = ${$gi_set->genomic_interval_set}[$i];	
            my $sequence_to_print = $gi->get_working_sequence();
			#print SEQUENCES_OUTFILE length($sequence_to_print)."\t".$sequence_to_print; # No need to print seq. length
			print SEQUENCES_OUTFILE $sequence_to_print;
			
			unless( ($i + 1) == scalar(@{$gi_set->genomic_interval_set}) ){
				print SEQUENCES_OUTFILE "\n";
			}
			
		}
		close SEQUENCES_OUTFILE;
		
		return $sequences_output_file;
	} # print_sequences_to_file #	
	
	method validate_sequence_model(Str $sequence_filepath, Int $order){
			
		# Generate random directory in which to run markov model c script
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		#$GU->user_info(3,"tempdir: $tempdir\n");

		#my $sequences_fp = $self->print_sequences_to_file($tempdir, $gi_set);
		my $path_to_executable = $APPLES_C_binaries."/learn_markov_modelREV493";
		my $sequence_model_fp = $tempdir."model.txt";
		#$GU->user_info(3,"$path_to_executable $sequence_filepath $order $sequence_model_fp\n"); 
		system( "$path_to_executable", "$sequence_filepath", "$order", "$sequence_model_fp") == 0 || die "System error!"; 
		my $model_from_file = $self->private_get_raw_model($sequence_model_fp);		
		rmtree($tempdir);
		return $model_from_file;		
	} # validate_sequence_model #
	
	method learn_markov_model(Genomic_Interval_Set $gi_set, Int $order){
		
		die "This method has been made redundant by generate_markov_model\n";
		use Data::Dumper;
		#my $order;# = $self->order;
		my $baseCount = 0;
		my $mer = '';
		my %markov;
		my $gi_count;
		my %ntfreq = ('A' => 0,
			      'T' => 0,
			      'G'=> 0,
			      'C'=> 0
		    );
		
		$GU->user_info(3,"Generating all kmers for markov model\n");
		# create all kmers and initialise the hash
		my $kmers = $self->get_all_possible_kmers($order);
		foreach my $kmer (@$kmers){
			
			$markov{$kmer} = {%ntfreq};
		}
		$GU->user_info(3,"Kmers generated\n");
		$GU->user_info(3,"Generating model\n");
		my $debug = 0;
		foreach my $gi (@{$gi_set->genomic_interval_set}){
			
			my $sequence = $gi->get_sequence();
			my @sequence_array = split('',$sequence);
			
			#if($gi_count > 5000){
			#	last;
			#}
			$gi_count++;
			#$GU->user_info(3,$gi_count.".\n");
			foreach my $base (@sequence_array){
				$baseCount++;
				if($baseCount>$order){
					
					if(length($mer) == $order){
						#$GU->user_info(3,$mer ." $base\n");
						
						
						if($markov{$mer}){
							#$GU->user_info(3,"$mer   $base DOES have entry\n");
							# Update nucleotide count for this MER
							# If base is A,T,G or C increase count, else ignore
							if($base =~ m/A|T|C|G/){
								
								$markov{$mer}{$base}++;
								# $GU->user_info(3,$markov{$mer}{$base} ."\n");
							}
							else{
								#$GU->user_info(3,"Illegal: $base\n");
							}
						}
						else{
							# Add this curently unique MER as hash key
							#$GU->user_info(3,"$mer   $base does not have entry\n");
							if($base =~ m/A|T|C|G/){
								
								$markov{$mer} = { %ntfreq };
								# If base is A,T,G or C increase count, else ignore
								$markov{$mer}{$base}++;
							}
							else{
								#$GU->user_info(3,"Illegal: $base\n");
							}
						}
						$mer = reverse($mer);
						chop($mer);
						$mer = reverse($mer);
					}
					else{
						#$GU->user_info(3,"mer is not being created properly\n");
						die "k-mer is not being created properly\n";
					}		
					$mer = $mer . $base;	
				}
				else{
					# Build initial $order-MER
					$mer = $mer . $base;
					#$GU->user_info(3,"Initial mer = $mer\n");
				}
			}
		}
		
		# Remove all kmers that contain non A, T, C or G characters
		while ( (my $key, my $val) = each %markov) {
				
			if ( $key =~ /[^ATCG]+/) {
				delete $markov{$key};
			}
			
		}
		my $c=0;
		while ( (my $key, my $val) = each %markov){
			
			if ( $key =~ /[^ATCG]+/) {
				$GU->user_info(3,$key."\n");
			}
			else {
				$c++;
			}
			
		}
		$GU->user_info(3,$c."\n");
		# Generate frequencies for markov model	
		while ( (my $key, my $val) = each %markov) {
			
			
			my $ntcount=0;
			# Get total frequency for each MER 
			while ( (my $k, my $v) = each %$val) {
				
				$ntcount += $v;
			}
			
			while ( (my $k, my $v) = each %$val) {
				
				if ($markov{$key}{$k} == 0) {
					
				}
				else {
					$markov{$key}{$k} = $v/$ntcount;
				}
			}		
		}
		
		return \%markov;
	} # learn_markov_model #
	
	method get_all_possible_kmers(Int $k){
			
		my $width = $k;
		my $max = 4;
		my @values;
		my @val;
		my @kmers;
		&_print_loop(\@values, $width, 0, $max, \@kmers);
		
		for (my $i = 0; $i <= $#kmers; $i++){
			$kmers[$i] =~ tr/0123/ACGT/;
		}
		@kmers = sort(@kmers);
		
		return \@kmers;	
		
		#------------------------------------------------
		sub _print_loop($$$$){
			
			my $v = shift;
			my @values = @{$v};
			my $width = shift;
			my $cur_col = shift;
			my $max = shift;
			my $KMERS = shift;
			my @kmers;
			my $word;
			
			if ($cur_col == $width){
				for (my $ i=0; $i < $width; $i++){
					$word .= $values[$i];	
				}
				push(@{$KMERS}, $word);
				
			}
			else{
				for (my $i=0; $i < $max; $i++){
					
					$values[$cur_col] = $i;
					&_print_loop(\@values, $width, $cur_col + 1, $max, $KMERS);	
				}
			}	
		} # _print_loop #

	} # get_all_possible_kmers #

    method private_get_raw_model(Str $resultsfile) { 
		
	my $FH;
	$GU->wait_for_file_existence ($resultsfile); # will wait up to 20 seconds for file to exist
	open FH,$resultsfile or die "Can't open $resultsfile";
	my %data;
	
	while ( <FH> ){    
	    chomp;
	    my @split = split, /\t/;
	    unless ( $split[0] =~ /\d/ ){
		$data{$split[0]} = $split[1];
	    }		
	}	
	close FH;		
	return \%data;
    } # private_get_raw_model #

} # Sequence_Model_Learner #



