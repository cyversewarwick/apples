### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### WM_Utilities Class ###

use MooseX::Declare;

class WM_Utilities {	
	use Storable qw (nstore retrieve);
	use Parameters;
	use BiFa_Server_Interface; 
	use Data::Dumper;
	use Generic_Sequence_Pattern;
	#  The following two 'uses' might theortically be required but in practise appear 
	#  not to be, because they are already 'used' by the Perl modules that 'use' WM_Utilities
	#  They have been commented out because their presence causes recursive 'use' loops which
	#  the Eclipse debugger cannot cope with.
	# use Generic_Pattern_Matching_Model;	
	# use Genomic_Data_Miner;
	use BiSi;
	use General_Utilities;
	use Genomic_Interval_Set;
	use APPLES_Datatypes qw (BiFaAlgorithm Boolean);
	use WM_Scores;
	use Sequence_Model_Learner;
	#use Error qw(:try);
	use File::Temp qw (tempdir);
	use File::Path qw (rmtree);
	use Cwd;
	use constant {FALSE => 0,
		      TRUE => 1};
	
	has 'bifa_if' => (is => 'rw', isa => 'BiFa_Server_Interface', clearer => 'clear_bifa_if');
	
	my $GU = General_Utilities->new();
	
	method create_counts_from_sequence_alignment( ArrayRef[Str] $aligned_sequences ){
		
		my @counts;
		my @sequence;
		my @counts_for_matrix; # This is returned
		
		# Check if sequences in alignment are the same length
		
		my $pssm_length = length( $$aligned_sequences[0] );
		my $all_sequences_same_length = TRUE;
		for my $index ( 1 .. $#$aligned_sequences ){
			
			if ( $pssm_length != length($$aligned_sequences[$index]) ){
				
				$all_sequences_same_length = FALSE;
			}
		}
		
		if ( $all_sequences_same_length eq FALSE ){
				
			die "This method cannot deal with sequence alignments where all sequences are not the same length\n";
		}
		else{
			
			for (my $i=0; $i<$pssm_length; $i++ ){
				my %nt_counts = ( A => 0,
									C => 0,
									G => 0,
									T => 0,
				);
				$counts[$i] = \%nt_counts;
			}
			
			my $rep=0;
			foreach ( @$aligned_sequences ){
				my @split = split('', $_);
				for (my $i=0; $i<scalar(@split); $i++ ){
					$counts[$i]->{$split[$i]}++;	
				}	
			}
			
			my $i = 0;
			foreach (@counts){
				
				$counts_for_matrix[$i][0] = $_->{A};
				$counts_for_matrix[$i][1] = $_->{C};
				$counts_for_matrix[$i][2] = $_->{G};
				$counts_for_matrix[$i][3] = $_->{T};	
				$i++;
			}
		}
		
		return @counts_for_matrix;
		
	} # create_counts_from_sequence_alignment #
	
	method get_server_version(){	
		my $sv = $self->bifa_if->ServerVersion();
		return $sv;		
	} # get_server_version #
	
	method get_custom_pssm_version(){			
		my $cpv = $self->bifa_if->CustomPssmVersion();
		return $cpv;
	} # get_custom_pssm_version #
	
	method get_transfac_version(){
		
		$self->private_establish_bifa_connection();
		my $tv = $self->bifa_if->TransfacVersion();
		return $tv;
	} # get_transfac_version #
	
	method get_transfac_version_parameter(){
		my $tvp = $self->get_transfac_version();
		$tvp .= ".".$self->get_custom_pssm_version();
		
		return $tvp;
	} # get_transfac_version_parameter #
	
	method get_all_bifa_pssms(Str $tvp, Boolean $transfac_only, Boolean $want_freq, Boolean $want_count){
		
		my $pssms = $self->private_retrieve_all_bifa_wm_objects_job_handler($tvp, $transfac_only, $want_freq, $want_count);
		
		return $pssms;
		
	} # get_all_bifa_pssms #	
		
	method get_bifa_server_wm_objects(ArrayRef[Str] $bifa_wm_id_list, Str $transfac_version , Boolean $want_freq, Boolean $want_count){
		my @wm_objects;
		my $count = 0;
		foreach my $wm_id (@$bifa_wm_id_list){
			
			my $bifa_server_wm = BiFa_Server_Weight_Matrix->new(wm_identifier => $wm_id, wm_source => 'BiFa', transfac_version => $transfac_version);
			
			if( $want_freq == 1 ){
				my $freqs = $self->private_get_bifa_server_wm_freq_dist($wm_id, $transfac_version);
				$bifa_server_wm->{wm_freq_dist} = $freqs;
				$bifa_server_wm->private_has_frequency_matrix();
					
			}
			if( $want_count == 1 ){
				my $counts = $self->private_get_bifa_server_wm_counts($wm_id, $transfac_version);
				$bifa_server_wm->{wm_counts} = $counts;				
			}
			else{
				# Do not add frequencies or counts
			}
			$count++;
			$GU->user_info(3,"Retrieved ".$count." PSSMs.\n");
			push(@wm_objects, $bifa_server_wm);
		}
		
		return \@wm_objects;		
	} # get_bifa_server_wm_objects #
		
	method get_bifa_pssm_info(Str $pssm_name, Str $tvp) {		
		#my $soapIf = $self->private_connect_to_bifa();
		$self->private_establish_bifa_connection();
		my $pssm_info = $self->bifa_if->PssmInfo($pssm_name, $tvp);
		
		return $pssm_info;
	} # get_bifa_pssm_info #
	
	method get_wm_length (Generic_Weight_Matrix $wm) {		
	    my $wm_length;
	    $wm_length = $wm->get_pssm_length();	
	    return $wm_length;
	} # get_wm_length #
	
	method list_bifa_server_pssms(Boolean $useConsensus, Str $species, Str $pssm_filter, ArrayRef $pssmSets, Str $tv) {
		
		$self->private_establish_bifa_connection();
		
		#my $tv = $self->bifa_if->TransfacVersion();
		
		my @bifa_if_pssms = $self->bifa_if->Pssms($useConsensus,$species,$pssm_filter,$tv,$pssmSets);
		
		my @pssms;
		
		foreach( @bifa_if_pssms){
				
			#push(@pssms, scalar($_));
			my $sub = substr($_, 0);
			push(@pssms, $sub);
		}
		
		return \@pssms;
	} # list_bifa_server_pssms #
	
	method print_markov_model_to_file_for_empirical_scoring(HashRef $model, Int $model_order, Str $model_file_path){
		die;	
		my @keys = keys %{$model};
		
		my @sorted_keys = sort(@keys);
		my $model_output_file = $model_file_path."sequence.model";
		open MODEL_OUTFILE, ">".$model_output_file or die "Cannot write model to file\n";
		$GU->user_info(3, "Printing to $model_output_file\n");
		#$GU->user_info(3,Dumper($wmpm->scoring_model->sequence_model_source->model->order));
		print MODEL_OUTFILE $model_order."\n";
		
		for (my $i = 0; $i < @sorted_keys; $i++){
			
			my $key = $sorted_keys[$i];
			print MODEL_OUTFILE $key."A"."\t".${$$model{$key}}{A}."\n";
			print MODEL_OUTFILE $key."C"."\t".${$$model{$key}}{C}."\n";
			print MODEL_OUTFILE $key."G"."\t".${$$model{$key}}{G}."\n";
			print MODEL_OUTFILE $key."T"."\t".${$$model{$key}}{T};
			
			unless( $i + 1 == scalar(@sorted_keys) ){
				print MODEL_OUTFILE "\n";
			}
		}
		close (MODEL_OUTFILE);
		
		return $model_output_file;
	} # print_markov_model_to_file_for_empirical_scoring #
	
	method new_print_markov_model_to_file_for_empirical_scoring(HashRef $model, Int $model_order, Str $model_file_path){
		
		my @keys = keys %{$model};
		
		my @sorted_keys = sort(@keys);
		
		foreach (@sorted_keys){
			$GU->user_info(3,$_."\n");
		}
		
		
		my $model_output_file = $model_file_path."sequence.model";
		open MODEL_OUTFILE, ">".$model_output_file or die "Cannot write model to file\n";
		$GU->user_info(3, "Printing to $model_output_file\n");
		#$GU->user_info(3,Dumper($wmpm->scoring_model->sequence_model_source->model->order));
		print MODEL_OUTFILE $model_order."\n";
		
		for (my $i = 0; $i < @sorted_keys; $i++){
			
			my $key = $sorted_keys[$i];
			
			print MODEL_OUTFILE $key."\t".$$model{$key}."\n";
			
			unless( $i + 1 == scalar(@sorted_keys) ){
				print MODEL_OUTFILE "\n";
			}
		}
		close (MODEL_OUTFILE);

		return $model_output_file;
	} # new_print_markov_model_to_file_for_empirical_scoring #
	
	method print_multiple_pssms_to_single_file( ArrayRef[Generic_Weight_Matrix] $GWMs, Str $outputdir){
		
		# Produce pssm file that is input to empirical scoring
		my $pssm_output_file = $outputdir."multiple_pssms";
		
		open PSSM_OUTFILE, ">".$pssm_output_file or die "Cannot write pssm to file\n";
		#$GU->user_info(3, "Printing to $pssm_output_file\n");
		
		foreach my $gwm ( @{$GWMs} ){

			print PSSM_OUTFILE $gwm->wm_identifier."\t".$self->get_wm_length($gwm)."\n";
			
			foreach my $row (@{$gwm->wm_freq_dist}){
				
				foreach my $col (@$row){
					print PSSM_OUTFILE $col."\t";
				}
				
				print PSSM_OUTFILE "\n";
			}
			#print PSSM_OUTFILE "X\n";
		}
		
		close PSSM_OUTFILE;
		
		return $pssm_output_file;
		
	} # print_multiple_pssms_to_single_file # 
	
	method print_pssm_to_file_for_empirical_scoring(Generic_Weight_Matrix $gwm, Str $outputdir){
			
		# Produce pssm file that is input to empirical scoring
		my $pssm_output_file = $outputdir.$gwm->wm_identifier.".pssm";
		
		open PSSM_OUTFILE, ">".$pssm_output_file or die "Cannot write pssm to file\n";
		$GU->user_info(3, "Printing to $pssm_output_file\n");
		print PSSM_OUTFILE $gwm->wm_identifier."\t".$self->get_wm_length($gwm)."\n";
		
		foreach my $row (@{$gwm->wm_freq_dist}){
			
			foreach my $col (@$row){
				print PSSM_OUTFILE $col."\t";
			}
			
			print PSSM_OUTFILE "\n";
		}
		
		close PSSM_OUTFILE;
		
		return $pssm_output_file;
	} # print_pssm_to_file_for_empirical_scoring #
	
	method print_pssm_profile_to_file_for_empirical_scoring(Weight_Matrix_Profile $WMP, Str $out_dir){
		
		my $pssm_profile_fp = $out_dir.$WMP->weight_matrix->wm_identifier.".profile";
		
		open PSSM_PROFILE, ">".$pssm_profile_fp or die "Cannot write pssm profile to file\n";
		
		print PSSM_PROFILE "ID\t".$WMP->weight_matrix->wm_identifier."\n";
		print PSSM_PROFILE "NI\t".$WMP->number_of_intervals."\n";
		
		my $i = 0;
		foreach(@{$WMP->hit_profile}){
				
			print PSSM_PROFILE "$i\t".$_."\n";
			$i++;
		}
		
		close PSSM_PROFILE;
		return $pssm_profile_fp;
	} # print_pssm_profile_to_file_for_empirical_scoring #
	
	method call_empirical_scoring_create_profile(Str $model_file_path, Str $pssm_file_path, Int $random_sequence_length, Int $burn_in, Int $seed, Str $running_dir, Empirical_Scoring_Method $ESM){
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_executable = $APPLES_C_binaries."/empirical_scoringREV830";
		
		$GU->user_info(3, "trying to create a weight matrix profile ... \n");
		
		
		my $prev_dir = getcwd();
		chdir($running_dir);
		
		if ( $ESM->isa('Biobase_Additive') ){
			my $entropy_param = $ESM->entropy_param;
			$GU->user_info(3, "$path_to_executable p $model_file_path $random_sequence_length $burn_in $pssm_file_path $seed $running_dir A $entropy_param\n");
			system("$path_to_executable", "p", "$model_file_path", "$random_sequence_length", "$burn_in", "$pssm_file_path", "$seed", "$running_dir", "A", "$entropy_param") == 0 || die "System error!";
			
		}
		elsif( $ESM->isa('Multiplicative') ){

		    die 'functionality not fully implemented yet';

			$GU->user_info(3, "$path_to_executable p $model_file_path $random_sequence_length $burn_in $pssm_file_path $seed $running_dir M\n");
			system("$path_to_executable", "p", "$model_file_path", "$random_sequence_length", "$burn_in", "$pssm_file_path", "$seed", "$running_dir", "M") == 0 || die "System error!";
		}
		else{
			die "There is currently no other type of empirical scoring method available\n";
		}
		
		chdir($prev_dir);	
	} # call_empirical_scoring_create_profile #
	
	method run_empirical_scoring(Weight_Matrix_Profile $WMP, Genomic_Interval_Set $gi_set){
		
		$GU->user_info(3, "Running empirical scoring\n");
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		
		my $pssm_profile_fp = $self->print_pssm_profile_to_file_for_empirical_scoring($WMP, $tempdir);
		my $pssm_fp = $self->print_pssm_to_file_for_empirical_scoring($WMP->weight_matrix, $tempdir);

		#$GU->user_info(3,"Testing ".$gi->label."\n");
		#push(my @single_gi_in_array, $gi);
		#my $single_gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@single_gi_in_array);
			
		my $sequences_fp = $self->print_sequences_to_file_for_empirical_scoring($gi_set, $tempdir);
			
		my $scores_fp = $self->call_empirical_scoring_score_sequences("s", $pssm_fp, $sequences_fp, $pssm_profile_fp, $tempdir);
			
		my $parsed_scores = $self->parse_empirical_scores($scores_fp);
		#$GU->user_info(3,Dumper($parsed_scores));
		unlink( $sequences_fp );
		unlink( $scores_fp);
		rmtree($tempdir);
			
		return $parsed_scores;
	} # run_empirical_scoring #
	
	method score_set_of_sequences(WM_Pattern_Model $WMPM, Weight_Matrix_Profile $WMP, Genomic_Interval_Set $gi_set){
		
		die "This method is not yet correctly implemented\n";

		# Sascha: did not amend call to overrepresentation_test below as there is this die-statement anyway

		my $data_miner = Genomic_Data_Miner->new();
		
		my @bp_array;
		
		foreach my $gi (@{$gi_set->genomic_interval_set}){
			
			#$GU->user_info(3,"Testing ".$gi->label."\n");
			push(my @single_gi_in_array, $gi);
			my $single_gi_set = Genomic_Interval_Set->new(genomic_interval_set => \@single_gi_in_array);
			my $raw_scores = $self->run_empirical_scoring($WMP, $single_gi_set);
			# Process scores
			my $wm_scores = WM_Scores->new(raw_result => $raw_scores, gi_set => $single_gi_set);
			#$GU->user_info(3,Dumper($wm_scores));
			my $l = $self->get_wm_length($WMPM->pattern);
			#$GU->user_info(3,"len: ".$l."\n");
			my $bp = $data_miner->overrepresentation_test($wm_scores, $WMPM, 3, $l);
			
			#nstore($bp, "/cluster/richardhickman/$gi->label.bpdump");
			push(@bp_array, $bp);
			#$GU->user_info(3,Dumper($bp));
			
		}
		
		return \@bp_array;
	} # score_set_of_sequences #
	
	method call_empirical_scoring_score_sequences(Str $scoring_switch, Str $pssm_filepath, Str $sequences_filepath, Str $profile_filepath, Str $running_dir, WM_Pattern_Model $WMPM){
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_executable = $APPLES_C_binaries."/empirical_scoringREV830";
		$GU->user_info(3,"SCORING\n");
		$GU->user_info(3, "$path_to_executable $scoring_switch $pssm_filepath $sequences_filepath $profile_filepath $running_dir");
		
		my $prev_dir = getcwd();
		chdir($running_dir);
		
		# Call empirical scoring Based on the choice of scoring metric (defined in WMPM)
		
		if ( $WMPM->scoring_model->scoring_method->isa('Biobase_Additive') ){


		my $entropy_param = $WMPM->scoring_model->scoring_method->entropy_param();
		$GU->user_info(3, " A $entropy_param\n");
		system("$path_to_executable", "$scoring_switch", "$pssm_filepath", "$sequences_filepath", "$profile_filepath", "$running_dir", "A", "$entropy_param") == 0 || die "System error!";

		}
		elsif ( $WMPM->scoring_model->scoring_method->isa('Multiplicative') ){
			
			system("$path_to_executable", "$scoring_switch", "$pssm_filepath", "$sequences_filepath", "$profile_filepath", "$running_dir", "M") == 0 || die "System error!";
        		$GU->user_info(3, " M\n");
		}
		else{
			
		}
		
		chdir($prev_dir);

		my $output = $sequences_filepath.".scores";
		
		return $output;
	} # call_empirical_scoring_score_sequences #
	
	method generate_PSSM_profile(HashRef $model, WM_Pattern_Model $wmpm){
		
		$GU->user_info(3, "GENERATING PSSM PROFILE\n");
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		$GU->user_info(3, "temporary directory for PSSM-profiling: $tempdir\n");
		
		my $profile = $self->call_PSSM_profile($model, $wmpm, $tempdir, $tempdir);
		$GU->user_info(3,ref($profile));
		
		$self->clear_bifa_if();
		
		unlink($tempdir);
		return $profile;
	} # generate_PSSM_profile #
	
	method call_PSSM_profile(HashRef $model, WM_Pattern_Model $wmpm, Str $model_output_filepath, Str $pssm_output_filepath){
		
		#my $model_output_file_name = "sequence_model.model";
		#my $pssm_output_file_name = "$outputdir.$wmpm->pattern->wm_identifier.".pssm";
		$GU->user_info(3, "Printing model to file\n");
		my $model_file_path = $self->new_print_markov_model_to_file_for_empirical_scoring($model, $wmpm->scoring_model->sequence_model_source->model->order, $model_output_filepath);
		$GU->user_info(3, "Printing PSSM to file\n");
		my $pssm_file_path = $self->print_pssm_to_file_for_empirical_scoring($wmpm->pattern, $pssm_output_filepath);
		
		# Now call empirical scoring and profile the PSSM:
		#$GU->user_info(3,Dumper($wmpm));
		my $random_sequence_length = $wmpm->scoring_model->training_sequence_length;
		my $burn_in = $wmpm->scoring_model->burn_in;
		my $seed = $wmpm->scoring_model->seed;
		#my $seed = 1080472374;
		
		#system("./emp_scoring_ds_mod", "p", "$model_output_file", "$random_sequence_length", "$burn_in", "$pssm_output_file", "$seed") == 0 || die "System error!"; 
		
		$self->call_empirical_scoring_create_profile($model_file_path, $pssm_file_path, $random_sequence_length, $burn_in, $seed, $model_output_filepath, $wmpm->scoring_model->scoring_method);
		my $profile_file_path = $pssm_output_filepath.$wmpm->pattern->wm_identifier.".profile";
		
		my $profile = $self->get_pssm_profile($profile_file_path, $wmpm->pattern);
		
		return $profile;
	} # call_PSSM_profile #
	
	method get_pssm_profile(Str $profile_filepath, Generic_Weight_Matrix $gwm){
			
		$GU->user_info(3, "Reading in PSSM profile\n");
		
		open FH, $profile_filepath or die "Cannot open the pssm profile\n";
		
		my $id;
		my $ni;
		my @hit_profile;
		while(<FH>){
				
			chomp;
			
			if ( m/^ID/){
					
				my @split = split('\t');
				$id = $split[1];
				
			}
			
			if ( m/^NI/){
					
				my @split = split('\t');
				$ni = $split[1];
				
				while(<FH>){
						
					chomp;
					my @split = split('\t');
					push(@hit_profile, $split[1]);
				}
			}
			
			
		}
		
		my $WMP = Weight_Matrix_Profile->new(weight_matrix => $gwm, number_of_intervals => $ni, hit_profile => \@hit_profile);
		close FH;
		return $WMP;
	} # get_pssm_profile #
	
	method print_sequences_to_file_for_empirical_scoring(Genomic_Interval_Set $gi_set, Str $sequences_output_file){
			
		# Produce sequence input file for empirical scoring
		$sequences_output_file = $sequences_output_file."sequences.in";
		$GU->user_info(3, "Printing sequence file to: $sequences_output_file\n");
		
		open SEQUENCES_OUTFILE, ">".$sequences_output_file or die "Cannot write sequences to file\n";
		#$GU->user_info(3,"->".scalar(@{$gi_set->genomic_interval_set}) ."\n");
		
		for (my $i = 0; $i < @{$gi_set->genomic_interval_set}; $i++){
			
			my $gi = ${$gi_set->genomic_interval_set}[$i];	
			print SEQUENCES_OUTFILE length($gi->get_working_sequence())."\t".$gi->get_working_sequence();
			
			unless( ($i + 1) == scalar(@{$gi_set->genomic_interval_set}) ){
				print SEQUENCES_OUTFILE "\n";
			}
			
		}
		close SEQUENCES_OUTFILE;
		
		return $sequences_output_file;
	} # print_sequences_to_file_for_empirical_scoring #
	
	method parse_empirical_processed_results(Str $results_input_file){
		
		open PROCESSED_RESULTS_FILE, $results_input_file or die "Cannot open file containing processed results\n";
		
		my @array_recs;
		
		while( <PROCESSED_RESULTS_FILE> ){
			
			chomp;
			
			my $ID;
			my $BP = 5;
			my @hit_details;
			
			my @split = split("\t", $_);
			
			if ($split[0] eq "ID"){
				$ID = $split[1];
				
				while( <PROCESSED_RESULTS_FILE> ){
						
					chomp;
					
					my @split = split("\t", $_);
					
					if ($split[0] eq "BP"){
						
						$BP = $split[1];
						
					}
					
					if ($split[0] eq "MI"){
						
						while( <PROCESSED_RESULTS_FILE> ){
							
							chomp;
							
							if($split[0] eq "XX"){
									
								last;
							}
							else{
								my @signif_hit_details = split ("\t", $_);
								
								push(@hit_details, \@signif_hit_details);
							}
							
						}
						last;
					}
						
				}
			}
			
			my $rec = {
				ID => $ID,
				BP => $BP,
				HITS => \@hit_details,
			};
			
			push(@array_recs, $rec);
		}
		my $i=0;
		foreach(@array_recs){
			
			$GU->user_info(3,Dumper($_));
			$i++;
			if( $i == 1){
				last;
			}
			
		}
		close PROCESSED_RESULTS_FILE;
	} # parse_empirical_processed_results #
	
	method parse_empirical_scores(Str $scores_input_file){
		
	    $GU->user_info(3,$scores_input_file."\n");
		open SCORES_INFILE, $scores_input_file or die "Cannot open file that contains empirical scores\n";
		
		my @all_empirical_scores;
		my $scores_count=0;
		while(<SCORES_INFILE>){
			
			chomp;
			#$GU->user_info(3,$scores_count."\n");
			my @scores = split(',', $_);
			#$GU->user_info(3,@scores."\n");
			push(@all_empirical_scores, \@scores);
			$scores_count++;
		}
		#$GU->user_info(3,Dumper($all_empirical_scores[18]));
		
		my @empirical_bifa_style_scores;
		
		
		# Turn the empirical scores into the same formate as bifa (this should change and is a quick hack for now)
		for (my $i = 0; $i < scalar(@all_empirical_scores); $i +=2 ){
			
			my @emp_scores;
			
			#$GU->user_info(3,scalar($all_empirical_scores[$i]));
			
			for (my $j = 0; $j < scalar(@{$all_empirical_scores[$i]}); $j++){
				
				push(@emp_scores, $all_empirical_scores[$i][$j]);
				push(@emp_scores, $all_empirical_scores[$i+1][$j]);
			}
			push(@empirical_bifa_style_scores, \@emp_scores);
		}
		
		#return \@all_empirical_scores;
		#$GU->user_info(3,Dumper($empirical_bifa_style_scores[18]));
		close SCORES_INFILE;

		return \@empirical_bifa_style_scores;
		
	} # parse_empirical_scores #
	
	method score_wm_on_sequences(WM_Pattern_Model $wmpm, Genomic_Interval_Set $gi_set){
		
		if( $wmpm->scoring_model->isa('BiFa_WM_Scoring_Parameters') ){
			
			if( $self->private_is_bifa_scoring_compatible($wmpm) ){
				$GU->user_info(3, "Can be scored by BiFa\n");
				$self->private_establish_bifa_connection();
				
				my $pssm_id = $wmpm->pattern->wm_identifier;
				
				my $raw_scores = $self->private_score_using_bifa_server($pssm_id, $gi_set, $wmpm->scoring_model->algorithm);
				
				return $raw_scores;
			}
			else{
				
				die "This pattern cannot be scored using the BiFa scoring method\n";
			}
		}
		elsif( $wmpm->scoring_model->isa('Empirical_WM_Scoring_Parameters') ){
			
			#die "Empirical scoring method is not yet implemented\n";
			
			if( $self->private_is_empirical_scoring_compatible($wmpm) ){
				$GU->user_info(3, "Can be scored by empirical scoring method\n");
				
				if ( defined($wmpm->pattern->wm_freq_dist) ){
					$GU->user_info(3,"Frequency matrix defined\n");
					
				}
				else {
					die "Empirical scoring requires that the pattern has a frequency matrix associated with it.\n";	
				}
				# Learn sequence model
				
				my $genome_db;
				my $model;
				if( $wmpm->scoring_model->sequence_model_source->isa('Actual_Sequence_Model') ){
					die "This functionality not yet implemented\n";
				}
				elsif( $wmpm->scoring_model->sequence_model_source->isa('Partial_Sequence_Model_Learning_Params') ){
					$GU->user_info(3, "Using Partial_Sequence_Model_Learning_Params\n");
					# Retrieve genome database depending on whether we explicitly or not. 
					my $genome_db;
					if( $wmpm->scoring_model->sequence_model_source->isa('Full_Sequence_Model_Learning_Params') ){
						$genome_db = $wmpm->scoring_model->sequence_model_source->genome_db();
					}
					else{
						$genome_db = ${$gi_set->genomic_interval_set}[0]->genome_db;
					}
					
					$GU->user_info(3,"Thinking about making the model now\n");
					my $SML = Sequence_Model_Learner->new();
					#$model = retrieve("/Volumes/cluster/richardhickman/april_apples_working_copy/generated_sequence_model.dump") or die;
					
					$model = $SML->learn_sequence_model_through_job_handler($wmpm->scoring_model->sequence_model_source,$genome_db);
				}
				else{
					die "The sequence_model_source type checking is not working correctly!\n";
				}
				
				# Generate random directory in which to run empirical scoring
				my $tempdir = $GU->get_temp_random_directory(FALSE);
				$GU->user_info(3, "$tempdir\n");
				#------
				my $cache = TRUE; # we want to cache the result
				my $cache_duration = 90; # number of days to store result
				my $job_handler = Job_Handler->new();
				$job_handler->get_config_settings();
				my $function = 'generate_PSSM_profile'; # name of the function you wish to call
				#push (my @parameters, $seed, $int, $genome_sequence_database_parameters); # make array of parameters required for the function you 	are calling, in the correct order!
				
				push(my @profile_parameters, $model, $wmpm);
				
				#foreach(@profile_parameters){ $GU->user_info(3,$_."\n");}
				
				my $profile_job_parameters = Job_Parameters->new(memory_requirement_high => FALSE,
										 wall_time_estimate => 172800,
										 cache_in_memory => TRUE
				); # make a Job_Parameters object 
				
				# CALL THE FUNCTION A THE JOB_HANDLER
				$self->clear_bifa_if();
				my $aggregate_exception = Aggregate_Exception->new(); # collects Job_Information_Exception objects thrown by Job_Handler			
				my @pssm_profile_cache_result = eval {
				    $job_handler->handle_APPLES_function($function, $self, \@profile_parameters, $cache, $cache_duration, $profile_job_parameters);
				};
				#deal with any Job_Information_Exception objects
				if ($@) {
				    my $exception_content = $@;
				    my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
				    # no parallelism here, though
				    $aggregate_exception = $outcome[1];
				    if ($outcome[0]) {
					rmtree($tempdir);
					$aggregate_exception->print_statistics_and_die;
				    } else {
					die 'did not expect this could happen at design time of this code.';
				    }
				}
				
				# THEN TO HANDLE THE RETURNED RESULT
				my $pssm_profile_object = shift(@pssm_profile_cache_result); # first element is the calling object, shift it from the array
				my $pssm_profile = $pssm_profile_cache_result[0];
				#------
				#my $pssm_profile = $self->call_PSSM_profile($model, $wmpm, $tempdir, $tempdir);
				my $profile_fp = $self->print_pssm_profile_to_file_for_empirical_scoring($pssm_profile, $tempdir);
				my $pssm_fp = $self->print_pssm_to_file_for_empirical_scoring($wmpm->pattern, $tempdir);
				
				my $sequences_fp = $self->print_sequences_to_file_for_empirical_scoring($gi_set, $tempdir);			

				my $scores_input_file;
				
				my $scoring_option = "s"; # Default for now
				
				my $raw_scores;
				
				if( $scoring_option eq "s"){
					
					$scores_input_file = $self->call_empirical_scoring_score_sequences("s", $pssm_fp, $sequences_fp, $profile_fp, $tempdir, $wmpm);
					
					$GU->user_info(3, "$scores_input_file\n");
					
					$raw_scores = $self->parse_empirical_scores($scores_input_file);
					
				}
				elsif( $scoring_option eq "f"){
					die "This functionality not implemented yet\n";
					$scores_input_file = $self->call_empirical_scoring_score_sequences("f", $pssm_fp, $sequences_fp, $profile_fp, $tempdir);	
				}
				
				rmtree($tempdir);
				#$GU->user_info(3,"*******************************".$tempdir."\n");

				return $raw_scores;
				
			    #my $raw_scores = $self->private_score_using_bifa_server($pssm_id, $gi_set, $wmpm->scoring_model->algorithm);
				
				#return $raw_scores;
			}
			else{
				
				die "This pattern cannot be scored using the Emprical scoring method\n";
			}
			
		}
		
		else{
			
			die "This weight matrix pattern model cannot be scored\n";
		}
	} # score_wm_on_sequences #
	
	method score_multiple_wm_on_sequences(ArrayRef[WM_Pattern_Model] $array_wm_objects, Genomic_Interval_Set $gi_set){
		
		# Check all WM are BiFa scoring compatible
		my @pssm_set;
		
		foreach my $wmpm ( @{$array_wm_objects} ){
			
			push(@pssm_set, $wmpm->pattern->wm_identifier);
			
		}
		
		$self->private_establish_bifa_connection();
		my $tv = $self->bifa_if->TransfacVersion();
		my $algorithm = BiFa_Server_Interface::ALG_OTT;
		
		my @seqs = $gi_set->return_unmasked_sequences();
		
		my @all_raw_scores = $self->bifa_if->scorePssmsOnSequences(\@pssm_set,\@seqs,$algorithm,$tv);
		
		return \@all_raw_scores;
		#die "Method not implemented yet\n";
		
	} # score_multiple_wm_on_sequences #
	
	method compute_wm_profiles(ArrayRef[Generic_Weight_Matrix] $pssms,HashRef $sequence_model,Empirical_WM_Scoring_Parameters $empirical_scoring_params) {
	    # computes profiles for all given weight matrices through Job_Handler, can run in parallel
	    # returns the profiles
	    
	    $GU->user_info(2,"Beginning of compute_wm_profiles.\n");
	    my @all_profiles;
	    my $aggregate_exception = Aggregate_Exception->new();
	    my $stats_required = FALSE;
	    foreach my $PSSM (@{$pssms}) {
		my $object = WM_Utilities->new();
		$object->clear_bifa_if(); # may not be necessary as $object was only just instantiated
		                          # (background to this: an existing bifa_if cannot be passed to
                              		  # to the compute nodes, the compute nodes need to contact the
		                          # BiFa-server themselves)
		my $function = 'generate_PSSM_profile';
		my $wmpm = WM_Pattern_Model->new(pattern => $PSSM, scoring_model => $empirical_scoring_params);
		my @parameters = ($sequence_model, $wmpm);
		my $high_memory = FALSE;
		my @result;
		eval {
		    @result = $GU->standard_call_to_job_handler_without_exception_handling($object,$function,\@parameters,$high_memory,FALSE);
		};
		if ($@) {
		    my $exception_content = $@;
		    my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
		    if ($outcome[0]) {
			$stats_required = TRUE;
		    }
		    $aggregate_exception = $outcome[1];
		} else {
		    push(@all_profiles,$result[0]);
		}
	    }
	    if ($stats_required) {
		$aggregate_exception->print_statistics_and_die;
	    }
	    $GU->user_info(2,"End of compute_wm_profiles.\n");
	    return @all_profiles;
	} # compute_wm_profiles #

	method get_all_bifa_pssm_set_names() {
	    $self->private_establish_bifa_connection();
	    my @result = $self->bifa_if->PssmSetNames();
	    return @result;
	} # get_all_bifa_pssm_set_names #

	method private_score_using_bifa_server(Str $pssm_id, Genomic_Interval_Set $gi_set, BiFaAlgorithm $algorithm) {
		# returns result (array) of pattern matching
		
		$self->private_establish_bifa_connection();
		my $tvp = $self->get_transfac_version_parameter();
		
		push (my @pssms, $pssm_id);
		
		if ($algorithm eq "ALG_OTT"){
			$algorithm = 0;
		}
		elsif ($algorithm eq "ALG_BAYESIAN"){
			$algorithm = 1;
		}
		else{
			die 'unknown BiFa scoring method';
		}
		
		my @seqs = $gi_set->return_unmasked_sequences();
		
		my @all_raw_scores = $self->bifa_if->scorePssmsOnSequences(\@pssms,\@seqs,$algorithm,$tvp);
		#$GU->user_info( 3, "raw scores from scorePssmsOnSequences:\n" );
		#$GU->user_info( 3,  Dumper (@all_raw_scores) );
		
		return \@all_raw_scores;
		
	} # private_score_using_bifa_server #
		
	method private_is_bifa_scoring_compatible(WM_Pattern_Model $wmpm){
		
		if ( $wmpm->pattern->isa('BiFa_Server_Weight_Matrix') ){
			
			return TRUE;
		}
		else{
			return FALSE;
		}
		
		
	} # private_is_bifa_scoring_compatible #
	
	method private_is_empirical_scoring_compatible(WM_Pattern_Model $wmpm){
		
		if ( $wmpm->pattern->isa('Generic_Weight_Matrix') ){
			
			return TRUE;
		}
		else{
			return FALSE;
		}
	} # private_is_empirical_scoring_compatible #
	
	method private_get_bifa_access_details() {
		use Config::General;
		# read APPLES config file + make hash of settings
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'}); 
		my %APPLES_config = $APPLES_conf->getall();
		# determine username and password from APPLES config file
		my $username = $APPLES_config{username};
		my $password = $APPLES_config{password};
		
		return ($username, $password);
		# then return the username and password!
	} # private_get_bifa_access_details #
	
	method private_connect_to_bifa() {		
		my ($username, $password) = $self->private_get_bifa_access_details();		
		package main;
		
		my $soapIf = new BiFa_Server_Interface();
				
		if ($soapIf -> init($username,$password))
		{
			$GU->user_info(3,"Connected to BiFa server\n");
			#$GU->user_info(3,"Used password: ".$password."\n");
			return $soapIf;
		}
		else
		{ 
			$GU->user_info( 1, "User validation failed" );
			die "Something wrong - is your BiFa-password correct?\n";
		}
	} # private_connect_to_bifa #
	
	method private_init_bifa_connection() {		
		my $bifa = $self->private_connect_to_bifa();
		$self->{bifa_if} = $bifa;
		#$GU->user_info(3,$self->bifa_if."\n");
	} # private_init_bifa_connection #
	
	method private_establish_bifa_connection(){		
		if (!defined($self->bifa_if)){
			
			$self->private_init_bifa_connection();
			#$GU->user_info(3,"bifa_if NOT defined -> initialising now\n");
		}		
		else{
			# Check connection is active
			#$GU->user_info(3,"Checking connection active\n");		
     		
			eval{
				$self->get_server_version();
			};
			if ($@){
				# Need to establish fresh bifa server interface
				$self->private_init_bifa_connection();
			}
			else{
				#$GU->user_info(3,"nothing wrong appparently!\n");
			}
			
		}
		
		return 1;		
	} # private_establish_bifa_connection #

	method private_get_bifa_server_wm_freq_dist(Str $pssmName, Str $transfacVersion) {		
		$self->private_establish_bifa_connection();
		my @freqs = $self->bifa_if->PssmFreqs($pssmName, $transfacVersion);
		
		return \@freqs;
	}# private_get_bifa_server_wm_freq_dist #
	
	method private_get_bifa_server_wm_counts(Str $pssmName, Str $transfacVersion){		
		$self->private_establish_bifa_connection();
		my @counts = $self->bifa_if->PssmCounts($pssmName, $transfacVersion);
		
		return \@counts;
	} # private_get_bifa_server_wm_counts #      

	method private_retrieve_all_bifa_wm_objects(Str $tvp, Boolean $transfac_only, Boolean $want_freq, Boolean $want_count){
		
		my @all_wms;
		
		if ( $transfac_only ){
			@all_wms = @{$self->private_get_all_transfac_pssms($tvp, $want_freq, $want_count)};
		}
		else{
			my $transfac_wms = $self->private_get_all_transfac_pssms($tvp, $want_freq, $want_count);
			my $custom_wms = $self->private_get_all_custom_pssms($tvp, $want_freq, $want_count);
			@all_wms = (@{$transfac_wms}, @{$custom_wms});
		}
		
		return \@all_wms;
		
	} # private_retrieve_all_bifa_wm_objects #
	
	method private_get_all_transfac_pssms(Str $tvp, Boolean $want_freq, Boolean $want_count){
		
		my @pssmSets = ();
		my $pssm_filter = ".";
		my $species = "PVIFNB"; 
		my $useConsensus = TRUE;
		
		@pssmSets = ('transfac'); #  ('transfac','all-custom')
		
		if ( !defined( $self->bifa_if() ) ){
			$self->private_establish_bifa_connection();
		}	
		
		my @pssms = $self->bifa_if->Pssms($useConsensus,$species,$pssm_filter,$tvp,\@pssmSets,);
		my $wms = $self->get_bifa_server_wm_objects(\@pssms, $tvp, $want_freq, $want_count);
		$self->clear_bifa_if();
		return $wms;

	} # private_get_all_transfac_pssms #
		
	method private_get_all_custom_pssms(Str $tvp, Boolean $want_freq, Boolean $want_count){

		my @pssmSets = ();
		my $pssm_filter = ".";
		my $species = ""; 
		my $useConsensus = TRUE;
		
		@pssmSets = ('all-custom'); #  ('transfac','all-custom')
		
		if ( !defined( $self->bifa_if() ) ){
			$self->private_establish_bifa_connection();
		}
		
		my @pssms = $self->bifa_if->Pssms($useConsensus,$species,$pssm_filter,$tvp,\@pssmSets,);
		my $wms = $self->get_bifa_server_wm_objects(\@pssms, $tvp, $want_freq, $want_count);
		$self->clear_bifa_if();
		return $wms;
	} # private_get_all_custom_pssms #
	
	method private_retrieve_all_bifa_wm_objects_job_handler(Str $tvp, Boolean $transfac_only, Boolean $want_freq, Boolean $want_count){
		
	    my $cache = TRUE;
	    my $cache_duration = 30;
	    my $job_handler = Job_Handler->new();
	    $job_handler->get_config_settings();
	    my $function = 'private_retrieve_all_bifa_wm_objects';
	    push (my @parameters, $tvp, $transfac_only, $want_freq, $want_count);
	    
	    my $job_parameters = Job_Parameters->new(memory_requirement_high => TRUE,
						     wall_time_estimate => 172800,
		                                     recache => FALSE);
	    $self->clear_bifa_if();
	    my $aggregate_exception = Aggregate_Exception->new();
	    my @cache_result = eval {
		$job_handler->handle_APPLES_function($function, $self, \@parameters, $cache, $cache_duration, $job_parameters);
	    };
	    if ($@) {
		my $exception_content = $@;
		my @outcome = $GU->standard_exception_handling_for_parallelisation($exception_content,$aggregate_exception);
		# no parallelism here, though
		$aggregate_exception = $outcome[1];
		if ($outcome[0]) {
		    $aggregate_exception->print_statistics_and_die;
		} else {
		    die 'did not expect this could happen at design time of this code.';
		}
	    }
	    shift @cache_result;

	    return $cache_result[0];
	} # private_retrieve_all_bifa_wm_objects_job_handler #	

}; # WM_Utilities #
