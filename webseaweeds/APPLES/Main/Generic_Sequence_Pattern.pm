### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Generic_Sequence_Pattern Class ###
# class for any type of sequence pattern
# This may be, for example, consensus sequences, weight matrices,
# graph-based motifs (MotifScan),  instrinsic DNA bending/conformational or
# physiochemical properties of DNA, 3D-structures of TF-DNA binding,
# variable length motifs, and motifs including optional bases and spacing rules.

use MooseX::Declare;

class Generic_Sequence_Pattern {

    method list_binding_molecules() {
	# return type is ArrayRef[Molecule]

	my @empty_array = ();
	return \@empty_array;
    } # list_binding_molecules #

} # Generic_Sequence_Pattern #

class K_mer extends Generic_Sequence_Pattern {

	# class not fully implemented yet
	method get_kmer_length {
		
		die "Method not implemented yet\n";	
	} # get_kmer_length #

} # K_mer #

class Generic_Weight_Matrix extends Generic_Sequence_Pattern{
	use APPLES_Datatypes qw (WeightMatrixSource JobHandlerIdentifier Boolean);
	use Statistics;
	use General_Utilities;
	use Exception;
	use Data::Dumper;
	use File::Path qw (rmtree);
	use constant {FALSE => 0,
		      TRUE	=> 1};
	
	has 'wm_identifier' => (is => 'rw', isa => 'Str', required => 1); # Trigger (trigger => \&private_clear_matrices) removed, as it prevents the adding of counts to a custom matrix. 
	has 'wm_source' => (is => 'ro', isa => WeightMatrixSource); # e.g. BiFa,
	has 'wm_freq_dist' => (is => 'rw', isa => 'ArrayRef[ArrayRef]', clearer => 'private_clear_wm_freq_dist');
	has 'wm_counts' => (is => 'rw', isa => 'ArrayRef[ArrayRef]', clearer => 'private_clear_wm_counts');
	has 'wm_name' => (is => 'rw', isa => 'Str');
	
	my $GU = General_Utilities->new();
	
	method print_to_custom_pssm_format( Str $outputdir ){

		my $custom_pssm_id = $self->wm_identifier();
		open FH, ">$outputdir/$custom_pssm_id.pssm" or die "Cannot open custom pssm output file\n";
		
		my $ID = $self->wm_identifier();
		my $NA = $self->wm_name();
		my $WI = $self->get_pssm_length();
		
		print FH "ID  $ID\n";
		print FH "NA  $NA\n";
		print FH "WI  $WI\n";
		print FH "PO  ";
		for ( 1 .. $WI ){
			
			if ( $_ < 10 ){
				
				print FH "0".$_." ";	
			}
			else{
				print FH $_." ";	
			}
		}  
		print FH "\n";
		
		my $aline = "CA  ";
		my $cline = "CC  ";
		my $gline = "CG  ";
		my $tline = "CT  ";
		
		my $position_counter = 0;
		
		foreach my $position ( @{$self->wm_counts()} ){
			
			$position_counter++;
			
			if ( scalar(@$position) !=4 ){
					
				die "This position in the array of counts does not have 4 elements\n";
			}
			else {
				
				$aline .="$$position[0] ";
				$cline .="$$position[1] ";
				$gline .="$$position[2] ";
				$tline .="$$position[3] ";
	
			}
		}
		
		print FH $aline."\n";
		print FH $cline."\n";
		print FH $gline."\n";
		print FH $tline."\n";
		
		close FH;
		
	} # print_to_custom_pssm_format #
	
	method generate_logo(Str $output_dir){

	    use Graphics;
	    
	    if (!$self->private_has_frequency_matrix) {
		die 'cannot produce a sequence logo if the frequency matrix is not defined.';
	    }
	    #$GU->user_info(3,Dumper($self->wm_freq_dist));
	    my @icm_matrix = $self->private_compute_icm();
	    my $pssm_id = $self->wm_identifier();
	    my $graphics = Graphics->new();
	    my $logo_filepath = $graphics->generate_PSSM_logo(\@icm_matrix, $pssm_id, $output_dir);
	    
	    return $logo_filepath;
	    
	} # generate_logo #
	
	method get_information_content() {
	    # returns information content for weight matrices that have a frequency matrix
	    
	    if (!$self->private_has_frequency_matrix) {
		die 'require frequency matrix to compute information content';
	    }
	    my $matrix = $self->wm_freq_dist;
	    my $result = 0;
	    foreach my $position (@{$matrix}) {
		my $tic = $self->private_compute_total_information_content_for_position($position);
		$result = $result + $tic;
	    }
	    return $result;
	} # get_information_content #

	method private_compute_icm(){
		
		my $matrix = $self->wm_freq_dist;
		my @icm;
		
		foreach my $position (@{$matrix}){
			
			my $tic = $self->private_compute_total_information_content_for_position($position);
			my @ic4p;
			
			foreach my $base (@{$position}){
				
				
				my $ic = $base * $tic;
				#$GU->user_info(3,$base." * ".$tic." = ". $ic." ic= ".$ic."\n");
				#$GU->user_info(3,$ic."\t");
				push(@ic4p, $ic);
			}
			push(@icm, \@ic4p);
			#$GU->user_info(3,"\n");
			
		}
		
		return @icm;
		
	} # compute_icm #
	
	method private_compute_total_information_content_for_position(ArrayRef $position) {
	    
	    my $addme = 0;
	    
	    foreach my $base (@{$position}){
		if ($base != 0) {
		    $addme += ( $base * ( log($base) / log(2) )  );
		}
	    }
	    
	    $addme = $addme + 2;
	    
	    return $addme;
	    
	} # private_compute_total_information_content_for_position #
	
	method generate_freq_dist_from_counts(WM_Pseudocount_Parameters $PCP){
		
		$GU->user_info(3,"generating FD for $self->wm_identifier\n");
		
		if ( $self->private_has_counts() ){
    		
			my @frequency_distribution;
			
			foreach (my $i = 0; $i < scalar(@{$self->wm_counts}); $i++ ) {
				my $addme = 0;
				
				my $position = $self->wm_counts()->[$i];
				foreach (@$position){
					
					$addme += $_;
				}
				# Add pseudocount to each count for each base
				
				for (my $j = 0; $j < scalar(@$position); $j++) {
					
					my $N = $addme;
					my $pseudocount = $self->private_get_pseudocount($PCP, $N);
					
					#$GU->user_info(3,$self->wm_counts()->[$i][$j] ."\t");
					$frequency_distribution[$i][$j] = ( $self->wm_counts()->[$i][$j] + $pseudocount ) / ($N + ( 4 * $pseudocount ) );
					#$GU->user_info(3,$frequency_distribution[$i][$j] ."\t");	
				}
				#$GU->user_info(3,"\n");
					
				$addme = 0;
				#push(@frequency_distribution, \@freqs);
				
			}
			#$GU->user_info(3,"Setting dist\n");
			$self->wm_freq_dist( \@frequency_distribution );
		}
		else{
			die "Weight matrix does not have any counts assigned to it\n";
		}
	} # generate_freq_dist_from_counts #     
	
	# Returns a hash containing the top n matches to $self in $wm_list. Key => value pair is 'pssm_id' => 'distance' 
	method get_top_n_closest_matches_in_list (ArrayRef[Generic_Weight_Matrix] $wm_to_test, Int $top_n){
		
		my @wm_list;
		push( @wm_list, $self );
		foreach ( @$wm_to_test ){
			push( @wm_list, $_);
		}
		
		unless( $self->private_has_frequency_matrix() ) {
		    die "Frequency matrix is NOT defined\n";
		}
		$GU->user_info(3,"Beginning of function get_top_n_closest_matches_in_list.\n");
		
		my $total_number_wms = @wm_list;
		
		
		my $WMU = WM_Utilities->new();
		my $GU = General_Utilities->new();
		# Sort out the location where intermediate files are going to be dumped
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		my $wms_fp = $WMU->print_multiple_pssms_to_single_file(\@wm_list, $tempdir);

		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_matrix_distance_c_exe = $APPLES_C_binaries."/pssm_distance_functionsREV789";
		#my $path_to_matrix_distance_c_exe = "/cluster/richardhickman/april_apples_working_copy/C_Programs/pssm_distance_functions";
		#print ("$path_to_matrix_distance_c_exe D $wms_fp $total_number_wms $tempdir");
		system("$path_to_matrix_distance_c_exe", "-f", "D", "-i", "$wms_fp", "-n", "$total_number_wms", "-o", "$tempdir", "-d", "k") == 0 || die "System error!";
	
		#`$path_to_matrix_distance_c_exe D $wms_fp $total_number_wms $tempdir`;
		
		open FH, $tempdir."/distances" or die "Cannot open pairwise distance matrix file\n";
		my @distances;
		
		while(<FH>){
				
			chomp;
			
			if ( /\d/ ){
				push(@distances, $_);
			}
			else{
			    $GU->user_info(1,"directory: ".$tempdir."\n");
			    $GU->user_info(1,"value: ".$_."\n");
			    die "Distance is not a number\n";		
			}
			
		}
		my %pssm_dist;	
		for ( my $i=0; $i<scalar(@distances); $i++){
			
			$pssm_dist{$$wm_to_test[$i]->wm_identifier} = $distances[$i];
		}
		
		if ($total_number_wms < $top_n){
			$top_n = $total_number_wms;
		}
		
		my @sorted_distances = sort { $pssm_dist{$a} <=> $pssm_dist{$b} } keys %pssm_dist; 
		my %distances_to_return;
		for (my $i=0; $i<$top_n; $i++){
			$distances_to_return{$sorted_distances[$i]} = $pssm_dist{$sorted_distances[$i]};
		}

		rmtree($tempdir);
		
		return \%distances_to_return;
		
	} # get_top_n_closest_matches_in_list #
	
	method identify_closest_n_known_motifs_in_list (ArrayRef[Generic_Weight_Matrix] $wm_list, Int $top_n){
		
		unless( $self->private_has_frequency_matrix() ) {
		    die "Frequency matrix is NOT defined\n";
		}
		
		my $total_number_wms = @{$wm_list};
		if ($total_number_wms < $top_n){
			$top_n = $total_number_wms;
		}
		my %closest_matrices_ranked;
		my %closest_matrices;
		my %pssm_distance;
		my @closest_motifs_ranked;	
		
		my @shortest_distances;
		foreach my $wm (@$wm_list) {

			my $dist = $self->get_shortest_distance_between_matrix($wm,1);
		
			$GU->user_info( 3, "$wm->wm_identifier -> $dist\n" );
			
			my $record = { WM => $wm, DIST => $dist };
			push(@shortest_distances, $record);
				
		}
		
		my $count = 0;
		
		
		my @sorted_scores;
		@sorted_scores = sort { $a->{DIST} <=> $b->{DIST} } @shortest_distances;
	
		
		for (my $i = 0; $i < $top_n; $i++) {
			
			push (@closest_motifs_ranked, $sorted_scores[$i]);
		}
		
		return \@closest_motifs_ranked;
	} # identify_closest_n_known_motifs_in_list # 
	
	method get_shortest_distance_between_matrix (Generic_Weight_Matrix $matrixB, Bool $include_revcomp){
		
		my @sense_score = $self->align_matrix($matrixB);
		
		# Find shortest distance using sense score
		my $shortest_distance = $sense_score[0];
		
		foreach (@sense_score) {
			
			if ($_ < $shortest_distance){
				$shortest_distance = $_;		
			}	
		}		
		
		# Use reverse complement by default
		
		if ( $include_revcomp == 1 ){
			
			# Get reverse complement of matrix B
			my $matrixB_revcomp = get_revcomp($matrixB);
			
			my @antisense_score = $self->align_matrix($matrixB_revcomp);
			
			# Find shortest distance
			foreach (@antisense_score) {
				if ($_ < $shortest_distance){
					$shortest_distance = $_;		
				}	
			}
		}	
		
		return $shortest_distance;
		
	} # get_shortest_distance_between_matrix #
	
	method align_matrix (Generic_Weight_Matrix $matrixB) {
		
		my @matrixA = @{$self->wm_freq_dist}; # Get frequency distribution this matrix (matrix A);
		my @matrixB = @{$matrixB->wm_freq_dist}; # Get frequency distribution for matrix B;
		my $n = @matrixA; # Get length of matrixA
		my $m = @matrixB; # Get length of matrixB
		
		#initialise matrix A
		
		@matrixA = @{append_Ns_to_start_of_matrix(\@matrixA, $m-1)};
		
		@matrixA = @{append_Ns_to_end_of_matrix(\@matrixA, $m-1)};
		
		my $Ax = $m-1;
		my $Ay = $Ax+$n;
		my $Bx;
		my $By = $Ax+1;# Initialise the position of the last base in B which is aligned with A; this will always be $Ax+1
		#$GU->user_info( 3, "Within the extended matrix, the original matrix A starts at position $Ax -> $Ay\n" );
		my @distances;
		
		for ( my $i = 0; $i <= $n+$m-2; $i++) {
			
			$Bx = $i;
			my $alignment_A;
			my $alignment_B;
			# Create 'modified matrix B' variable by appending appropriate N's to start/end of matrix
			my @mod_B = @matrixB;
			
			if ( ($By > $Ax) && ($By < $Ay) ) {
				my $ns = $Ay-$By;
				@mod_B = @{append_Ns_to_end_of_matrix(\@mod_B, ($Ay-$By))};
				
			}
			elsif ($Bx > $Ax) {
				my $ns = $Bx-$Ax;
				@mod_B = @{append_Ns_to_start_of_matrix(\@mod_B, ($Bx-$Ax))};
			}
			else {
				my $ns = 0;
				
				@mod_B = @mod_B;
				
			}
			my @sub_A;		
			
			#Take submatrix of A to align with B
			if ($Bx > $Ax) {
				@sub_A = @matrixA[$Ax..($Ax + $#mod_B)];
				
				my $s = $#mod_B - $#matrixB;
				my $j = ($Bx-$Ax);
				
			}
			else {
				@sub_A = @matrixA[$i..($i + $#mod_B)];
			}
			
			
			$By++; # update position of the last element in matrix B aligned with respect to matrix A
			
			# Calculate local matrix score for matrixB against matrixA at position i
			
			my $local_score_p_q = $self->matrix_distance(\@sub_A, \@mod_B);
			my $local_score_q_p = $self->matrix_distance(\@mod_B, \@sub_A);
			my $average_local_score = ( $local_score_p_q + $local_score_q_p ) / 2 ;
			#$GU->user_info(3,"(P,Q) = $local_score_1 | (Q,P) = $local_score_2\n");
			push(@distances, $average_local_score);
			
		}# for each i
		
		return @distances;
		
		
		# ------------------------------------------------
		# -------- SUBROUTINES FOR align_matrix() --------
		# ------------------------------------------------
		sub append_Ns_to_start_of_matrix{
			
			my $a = shift;
			my $num = shift;
			my @mod_a = @$a;
			
			for (my $i = 0; $i < $num; $i++) {
				
			    my @N = (0.32,0.18,0.18,0.32);
			    unshift(@mod_a, \@N);
				
			}
			
			return \@mod_a;	
		}# append_Ns_to_start_of_matrix()
		
		sub append_Ns_to_end_of_matrix{
			
			my $a = shift;
			my $num = shift;
			my @mod_a = @$a;
			
			for (my $i = 0; $i < $num; $i++) {
				
			    my @N = (0.32,0.18,0.18,0.32);
			    push(@mod_a, \@N);
				
			}
			
			return \@mod_a;	
		}# append_Ns_to_end_of_matrix()
		
	} # align_matrix #
	
	method matrix_distance (ArrayRef $fd1, ArrayRef $fd2){
		#warn "Local Matrix Distance method not fully implemented yet\n";
		
		my $distance = 0;
		my $N = @$fd1;
		my $M = @$fd2;
		my $stats = Statistics->new();
		
		if($N != $M){	
		    die "Distance cannot be computed for two matrices of different lengths;\nmatrix A = $N, matrix B = $M";
		}
			
		for (my $i = 0; $i < $N; $i++){
		    $distance += $stats->kullback_leibler_distance($$fd1[$i], $$fd2[$i]);	
		} 
		
		return $distance;
		
	} # matrix_distance # 
	
	method get_revcomp (){
		
		#warn "Getting reverse complement of matrix\n";
		my @comp_freq_dist;
		
		foreach my $position (@{$self->wm_freq_dist}) {	
			my @rc = ($$position[3], $$position[2], $$position[1], $$position[0]);
			push(@comp_freq_dist, \@rc);
		}
		my @revcomp_freq_dist = reverse(@comp_freq_dist);
		my $new_id = $self->wm_identifier."_REVCOMP";
		
		my $revcomp_wm = Generic_Weight_Matrix->new(wm_identifier=>$new_id, wm_freq_dist=>\@revcomp_freq_dist);
		
		return $revcomp_wm;
	} # get_revcomp #
	
	method get_consensus (ArrayRef[ArrayRef[Num]] $freq){
		
		my $string;
		
		foreach my $position (@$freq) {
			
			if( ($$position[1] && $$position[2] && $$position[3]) < $$position[0]){ $string .='A';}# If A is highest
			elsif( ($$position[0] && $$position[2] && $$position[3]) < $$position[1]){ $string .='C';}# If C is highest
			elsif( ($$position[0] && $$position[1] && $$position[3]) < $$position[2]){ $string .='G';}# If G is highest
			elsif( ($$position[0] && $$position[1] && $$position[2]) < $$position[3]){ $string .='T';}# If T is highest
			elsif( ($$position[0] && $$position[1] && $$position[2] && $$position[3]) == 0.25){ $string .='N';}# If A or C or G or T are equal
			elsif( ($$position[2] == $$position[3]) && ($$position[2] > ($$position[0] && $$position[1])) ){ $string .='K';}# If G or T are joint highest
			elsif( ($$position[0] == $$position[1]) && ($$position[0] > ($$position[2] && $$position[3])) ){ $string .='M';}# If A or C are joint highest
			elsif( ($$position[0] == $$position[2]) && ($$position[0] > ($$position[1] && $$position[3])) ){ $string .='R';}# If A or G are joint highest
			elsif( ($$position[1] == $$position[2]) && ($$position[1] > ($$position[0] && $$position[3])) ){ $string .='S';}# If C or G are joint highest
			elsif( ($$position[0] == $$position[3]) && ($$position[0] > ($$position[1] && $$position[2])) ){ $string .='W';}# If A or T are joint highest
			elsif( ($$position[1] == $$position[3]) && ($$position[1] > ($$position[0] && $$position[2])) ){ $string .='Y';}# If C or T are joint highest
			
			else{ $string .= "N";}	
		}
		
		return $string;
		
	} # get_consensus #
	
	method get_pssm_length () {
	    my $wm_length;
	    
	    if (!$self->private_has_frequency_matrix()) {
		if ($self->private_has_counts()) {
		    $wm_length = @{$self->wm_counts};
		} else {
		    die "cannot establish weight matrix length.";
		}
	    } else {
		$wm_length = @{$self->wm_freq_dist};
	    }
	    return $wm_length;
	} # get_pssm_length #

	method private_get_bifa_access_details () {
		use Config::General;
		# read APPLES config file + make hash of settings
		my $APPLES_conf = new Config::General("APPLES.dat");
		my %APPLES_config = $APPLES_conf->getall();
		# determine username and password from APPLES config file
		my $username = $APPLES_config{username};
		my $password = $APPLES_config{password};
		
		return ($username, $password);
		# then return the username and password!
	} # private_get_bifa_access_details #
	
	method private_connect_to_bifa () {
		
		my ($username, $password) = $self->private_get_bifa_access_details();
		
		package main;
		
		my $soapIf = new BiFa_Server_Interface();
		
		
		if ($soapIf -> init($username,$password))
		{
			#$GU->user_info(3,"Connected to BiFa server\n");
			return $soapIf;
		}
		else
		{ 
			$GU->user_info( 1, "User validation failed" );
			die "Something wrong\n";
		}
	} # private_connect_to_bifa #

	method private_has_frequency_matrix (){
			
	    if (defined $self->wm_freq_dist) {
		return TRUE;
	    }
	    else {
		return FALSE;
	    }
	} # private_has_frequency_matrix #
	
	method private_has_counts () {		
	    if (defined $self->wm_counts) {
		return TRUE;
	    }
	    else {
		return FALSE;
	    }
	} # private_has_counts #

	method private_get_pseudocount(WM_Pseudocount_Parameters $PCP, Num $N){
	    
	    my $pseudocount;
	    
	    if( $PCP->pseudocount_type() eq "p_over_n" ){
		$pseudocount = $PCP->p() / $N;
	    }
	    elsif( $PCP->pseudocount_type() eq "p_over_sqrt_n" ){
		$pseudocount = $PCP->p() / sqrt($N);
	    }
	    elsif( $PCP->pseudocount_type() eq "p" ){
		$pseudocount = $PCP->p();
	    }
	    else{
		die "DEBUG me please -> we should never have reached this point\n";
	    }
	    
	    return $pseudocount;
	} # private_get_pseudocount #

	sub private_clear_matrices() {
	    my ($self, $arg) = @_; 

	    $self->private_clear_wm_freq_dist(); # clearer methods are provided by Moose if a name is provided (as done above)
	    $self->private_clear_wm_counts();
	} # private_clear_matrices #

} # Generic_Weight_Matrix #

class BiFa_Server_Weight_Matrix extends Generic_Weight_Matrix {
    use Molecule;
    use constant {FALSE => 0,
		  TRUE	=> 1};
    use Data::Dumper; 
    
    has 'url' => (is => 'ro', isa => 'Str');
    has 'transfac_version' => (is => 'ro', isa => 'Str', required => TRUE);
    
    method set_frequency_matrix (){			
	my $wm_util = WM_Utilities->new();		
	my $tv = $wm_util->get_transfac_version();
	my $pssm_freq = $wm_util->private_get_bifa_server_wm_freq_dist($self->wm_identifier,$tv);
	$self->{wm_freq_dist} = $pssm_freq;
    } # set_frequency_matrix #
    
    override private_has_frequency_matrix () {		
	# Retrieve frequency distribution for BiFa WM if it does not exist.
	
	if (!defined $self->wm_freq_dist){
	    $self->set_frequency_matrix();	    
	}
	if (!defined $self->wm_freq_dist) {
	    die 'failed to obtain frequency matrix from BiFa-server (check-point 1).';
	}
	if (!defined ${$self->wm_freq_dist}[0]) {  # Sascha: this can happen
	    print Dumper($self->wm_freq_dist);
	    print Dumper($self->wm_identifier);
	    die 'failed to obtain frequency matrix from BiFa-server (check-point 2).';
	}
	return TRUE;
    } # private_has_frequency_matrix #
    
    override list_binding_molecules () {
	my $wm_utils = WM_Utilities->new();
	my $pssm_info = $wm_utils->get_bifa_pssm_info($self->wm_identifier, $self->transfac_version);
	my $factors = $pssm_info->{'factors'}->{'string'};
	my @molecules;
	my %map = $self->private_get_transfac_to_uniprot_map();
	foreach my $factor (@$factors) {
	    $factor  = substr($factor,0,6);
	    my @uniprot_factors = $self->private_get_uniprot_ac_for_transfac_ac($factor, \%map);
	    foreach my $uniprot_factor (@uniprot_factors) {
		my $molecule = Protein->new(uniprot_id => $uniprot_factor);
		push(@molecules, $molecule);
	    }
	}	
	return @molecules;
    } # list_binding_molecules #
    
    method private_get_transfac_to_uniprot_map () {
	my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
	my %APPLES_config = $APPLES_conf->getall();
	my $mapfile = $APPLES_config{transfac_to_uniprot_map};
	my %transfac_uniprot_map;
	open(INFILE, $mapfile)
	    or die "Can't open file $mapfile\n";
	# skip the first line
	my $line = <INFILE>;
	while ($line = <INFILE>) {
	    my ($tid, $upid) = split("\t", $line);
	    if ($upid){
		$transfac_uniprot_map{$tid} = $upid;
	    }
	    else {
		$transfac_uniprot_map{$tid} = 'NONE';
	    }
	}
        close(INFILE);
        return %transfac_uniprot_map;
    } # private_get_transfac_to_uniprot_map #
    
    method private_get_uniprot_ac_for_transfac_ac (Str $transfac_factor_ac, $transfac_uniprot_map) {       
	my @result;
	my $uniprot_ac = $transfac_uniprot_map->{$transfac_factor_ac};
	if (defined $uniprot_ac) {
	    if (($uniprot_ac ne 'NONE')&&
		($uniprot_ac ne '')) {
		push(@result,$uniprot_ac);
	    }
	}
	return @result;
    } # private_get_uniprot_ac_for_transfac_ac #
    
} # BiFa_Server_Weight_Matrix #

class Generic_Pair_Weight_Matrices extends Generic_Sequence_Pattern{
	
	has 'wm_a' => (is => 'ro', isa => 'Generic_Weight_Matrix', required => 1);
	has 'wm_b' => (is => 'ro', isa => 'Generic_Weight_Matrix', required => 1);
	has 'min_spacing' => (is => 'ro', isa => 'Int', required => 1);
	has 'max_spacing' => (is => 'ro', isa => 'Int', required => 1);
	
	method get_pattern_length (){
			
		die "Method not implemented yet\n";
	} # get_pattern_length #

} # Generic_Pair_Weight_Matrices #

class Weight_Matrix_Profile extends Generic_Sequence_Pattern{
	
	has 'weight_matrix' => ( is => 'ro', isa => 'Generic_Weight_Matrix', required => 1);
	has 'number_of_intervals' => (is => 'ro', isa => 'Int', required => 1);
	has 'hit_profile' => (is => 'ro', isa => 'ArrayRef[Int]', required => 1);
		
} # Weight_Matrix_Profile #
