=pod

=head1 Class Statistics

A collection of statistical utilities.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<entropy>

Returns an entropy/information content associated with provided distribution.

=item C<kullback_leibler_distance>

Returns Kullback-Leibler distance between 2 distributions.

=back

=head1 LICENSE

This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package and distributed 
under Academic Non-Commercial Use Licence.

=head1 COPYRIGHT

(c) Copyright University of Warwick 2009-2010

=head1 AUTHOR

=cut

use MooseX::Declare;

class Statistics {
	use Parameters;
	use General_Utilities;
	use File::Temp qw (tempdir);
	use File::Path qw (rmtree);
	use Data::Dumper;
	use APPLES_Datatypes qw (Boolean);
	use constant {FALSE => 0,
		      TRUE  => 1};	

	my $GU = General_Utilities->new();
	
=pod
    
=head2 entropy( ArrayRef $p )

Calculate entropy/information content of a distribution. See http://mathworld.wolfram.com/Entropy.html
$p is assumed to be NxM where M <= 4 and N can be anything. The distribution is summed over rows, 
then M colomns are added together for overall entropy value.

Return: Returns entropy of distribution p.

=cut	
	method entropy(ArrayRef $p) {
		
		print Dumper($p);
		my $N = @$p;
		my $M = @{$$p[0]};
		
		my @entropy = (0,0,0,0);		
		
		my $log2 = log(2);
		
		my $entropy = 0;
		
		for(my $i = 0; $i < 1; $i++) {
			for(my $j = 0; $j < $M; $j++) {
				$entropy[$j] -= $$p[$i][$j] * log($$p[$i][$j]) / $log2;
			}
		}
		
		my $result = 0;
		foreach my $ent(@entropy) {
			$result += $ent;
		}
		
		return $result;		
	}
	
	
	# Kullback_Leibler_Distance calculates the relative entropy between two discrete 
	# probability distributions p and q.
	method kullback_leibler_distance( ArrayRef $p, ArrayRef $q){
		
		my $d = 0;
		my $N = @$p;
		#$GU->user_info( 3, @$p );
		#$GU->user_info( 3, "\n+++++\n" );
		#$GU->user_info( 3, @$q );
		#$GU->user_info( 3, "\n-----\n" );
		
		for (my $i = 0; $i < $N; $i++){
			if ($$p[$i] == 0){ $$p[$i] = 0.000000001;}# Problem with log(0)!
			if ($$q[$i] == 0){ $$q[$i] = 0.000000001;}
			
			$d += $$p[$i] * ( log($$p[$i] / $$q[$i]) / log(2) );
			
		} 
		#die "-----------\n";
		#$GU->user_info( 3, $d."\n" );
		return $d;		
	} # kullback_leibler_distance #
	
	method get_binomial_probability(Int $n, Int $k, Num $p) {# Changed $p type form Int to Num
	    if ( $p == 0 ) {	# If the p-value is 0...
		if ( $k == 0 ) {	# AND the next position is 0...
		    return 1;		# Return 1.
		}
		else {				# But if the next site exists
		    return 0;		# Return 0
		}
	    }
	    if ( $p == 1 ) {	# If the p-value is 1
		return 1;			# Return 1
	    }		
	    my @firstsum;		# If the p-value NE 0 or 1, do stuff.
	    my @secondsum;
	    $#firstsum = $n;
	    $#secondsum = $n;
	    for ( my $i = 0 ; $i <= $n ; $i++ ) {
		if ( $i == 0 ) {
		    $firstsum[$i] = 0;
		    $secondsum[$i] = 0;
		}
		else {
		    $firstsum[$i] = $firstsum[ $i - 1 ] + log( ( $n - $i ) + 1 );
		    $secondsum[$i] = $secondsum[ $i - 1] + log($i);
		}
	    }
	    my $result = 0;
	    for ( my $i = $n ; $i >= $k ; $i-- ) {
		my $addme = $i * log($p);
		#$GU->user_info(3,"1. $addme\n");
		$addme += ( $n - $i ) * log( 1 - $p );
		#$GU->user_info(3,"2. $addme\n");
		$addme -= $secondsum[$i];
		#$GU->user_info(3,"3. $addme\n");
		$addme += $firstsum[$i];
		#$GU->user_info(3,"4. $addme\n");
		$addme = exp($addme);
		#$GU->user_info(3,"5. $addme\n");
		$result += $addme;
		#$GU->user_info(3,"R. $result\n\n");
		
	    }
	    #$GU->user_info(3,"RESULT = ".$result."\n");
	    return $result;
	} # get_binomial_probability #
	
	method best_binomial_pvalue(ArrayRef $n_best_pvalues, Int $max_sites) {
		
		my $best_p_value = 100;
		my $best_p;
		my $p;
		
		for ( my $i = 0; $i < scalar(@$n_best_pvalues); $i++) {
			
			$p = $self->get_binomial_probability($max_sites, $i+1, ${$n_best_pvalues}[$i]{PVALUE});
			
			$GU->user_info( 3, "p value: ".$p."\n" );
			#$GU->user_info(3,"p value $p\n");
			if ($p < $best_p_value) {
				$best_p_value = $p;
				$best_p = {
					PVALUE => $best_p_value,
					INDEX => $i,
					N_BEST_PVALUES => $n_best_pvalues
				};				
			}
		}	
		return $best_p;		
	} # best_binomial_pvalue #
	
	method phyper(Int $q, Int $m, Int $n, Int $k, Boolean $lower_tail){
	    # calls function of same name in R, see R-documentation for description

	    $GU->user_info(3,"phyper input q: ".$q."\n");
	    $GU->user_info(3,"phyper input m: ".$m."\n");
	    $GU->user_info(3,"phyper input n: ".$n."\n");
	    $GU->user_info(3,"phyper input k: ".$k."\n");
	    $GU->user_info(3,"phyper input lower_tail: ".$lower_tail."\n");

		# Find location of R binary
		my $APPLES_DAT = $ENV{'APPLES_DAT'};
		my $APPLES_conf = new Config::General($APPLES_DAT);
		# read APPLES config file
		my %APPLES_config = $APPLES_conf->getall();
		# determine job handler config file from APPLES config file
		my $job_handler_config_file = $APPLES_config{job_handler_config};
		# read job handler config file
		my $job_handler_conf = new Config::General($job_handler_config_file);
		
		#my $job_handler_conf = new Config::General($job_handler_config_file);
		# set config parameters, e.g. location of binaries
		my %job_handler_parameters = $job_handler_conf->getall();
		my $R_binary = $job_handler_parameters{R_binary};
		
		my $non_system_temp_dir = $APPLES_config{non_system_temp_dir};
		
	    	my $tempdir = $GU->get_temp_random_directory(FALSE);
		
		my $r_script = $tempdir."/foo.txt";
		$GU->user_info(3,$r_script."\n");
		
		my $lower_tail_string;
		if( $lower_tail ){
		    $lower_tail_string = "lower.tail = TRUE";		    
		}
		else{
		    $lower_tail_string = "lower.tail = FALSE";
		    $q = $q - 1;
		}
		
		my $Rcommand = "phyper($q,$m,$n,$k,$lower_tail_string,log.p = FALSE);";
		
		open RSCRIPT, ">".$r_script or die "Cannot create R input script at $r_script\n";
		
		print RSCRIPT $Rcommand;
		close RSCRIPT;
		my $r_command = "$R_binary --vanilla < $r_script > $tempdir/foo.results";
		
		$GU->user_info(3,$r_command."\n");
		
		#`$R_binary --vanilla < $r_script > $tempdir/foo.results`; #or die "Cannot run R command\n";
		my $result = `$R_binary --vanilla < $r_script > $tempdir/foo.results`;
		
		open RESULTS, "$tempdir/foo.results" or die "Cannot open the R results file\n";
		
		my $pvalue;
		
		while (my $input_line = <RESULTS>){
			
			chomp $input_line;
			
			if($input_line =~ m/\[/){
				
				my @split = split('\s+', $input_line);
				$pvalue = $split[1];
			}
			
		}
		#$GU->user_info(3,$pvalue."\n");
		close RESULTS;
	        rmtree($tempdir);
		return $pvalue;	
	} # phyper #
	
	
	method one_sided_ks_test(ArrayRef $data, ArrayRef $range, Str $type){
		
		#Added by PEB, June 2012. Tests whether a list of numbers is distributed between two specified values
		#according to the specified theoretical distribution, using one sided Kolmogorov-Smirnov test
		#Used to test if the distribution of motif start positions in a promotor region is random, ie linear
		
		#data - measured distribution to be tested
		#range - ascending list of possible values of theoretical distribution
		#type - expected shape of  theoretical distribution. 'linear' means data randomly distributed over range
		my @range = @$range;
		my $lower = $range[0];
		my $upper = $range[$#range];
		my $distribution_type = lc($type);
		
		if($lower >= $upper){
			die("ERROR: statistics::one_sided_ks_test: Invalid cdf range specified\n");	
		}
		
		my @sorted_data = sort {$a <=> $b} @$data;
		my $num_values = scalar (@sorted_data);
		
		if(($sorted_data[0] < $lower) || ($sorted_data[$num_values-1] > $upper)){
			
			if($sorted_data[0] < $lower){
				die("ERROR: statistics::one_sided_ks_test: '".$sorted_data[0]."' less than '".$lower."'");
			}
			if($sorted_data[$num_values-1] > $upper){
				die("ERROR: statistics::one_sided_ks_test: '".$sorted_data[$num_values-1]."' greater than '". $upper."'");
			}
			
			
			
			die("ERROR: statistics::one_sided_ks_test: Data falls outside cdf range specified\nData: ".$sorted_data[0]." to ". $sorted_data[$num_values-1].", Range ".$lower." to ".$upper."\n");	
		}
		
		$GU->user_info(3, 'K-S test sorted input data\n');
		$GU->user_info(3, Dumper(\@sorted_data));
		$GU->user_info(3, 'Distribution over range '.$lower.' to '.$upper.' expected to be '.$distribution_type.'\n');

		#Compare distributions to find largest difference
		my $D_value = 0;
		my $idx = 0;
		
		for(my $i = 0; $i <= $#range; $i++){
			my $cumulative_expected;
			if($distribution_type eq "linear"){
				#The expected cumulative distributon if linear
				$cumulative_expected = ($i+1)/($#range+1);
			}else{
				die("ERROR: statistics::one_sided_ks_test: Invalid distribution type ".$distribution_type." specified\n");		
			}
			
			#This is the proportion of values in data expected to be <= $range[$i]
			#find actual proportion
			my $val = $range[$i];
			while($idx < $num_values){
				if($sorted_data[$idx] <= $val){
					$idx++;
				}else{
					last;
				}
			}
			my $cumulative_observed = $idx/$num_values;
			
			if(abs($cumulative_observed-$cumulative_expected) > $D_value){
				#record largest difference between the two
				$D_value = 	abs($cumulative_observed-$cumulative_expected);
			}	
		}
		$GU->user_info(3, 'Test statistic: '.$D_value.'\n');
		
		#Now calculate the corresponding probability of getting this deviation from linear
		#according to algorithm from Marsaglia et al, Journal of Statisitcal Software, Vol 8:18
		
		my $probability_value;
		my $s = $D_value**2 * $num_values;

		if(($s>7.24) || (($s>3.76 && ($num_values>99)))){
			#shortcut for big D and high values of n, ie highly siginificant p value
			$probability_value = 2 * exp(-(2.000071+0.331/sqrt($num_values)+1.409/$num_values) * $s)
		}else{
			my $k = int($num_values * $D_value) + 1;
			my $m = int(2 * $k - 1);
			my $h = $k - $num_values * $D_value;
		
			my @H = ();
			my @Q = ();
			my $eQ;
	
			for(my $i = 0; $i < $m; $i++){
				for(my $j = 0; $j < $m; $j++){
					if($i-$j+1 < 0){
						$H[$i*$m+$j] = 0;
					}else{
						$H[$i*$m+$j] = 1;
					}
				}
			}
	
			for(my $i = 0; $i < $m; $i++){
				$H[$i*$m] -= $h**($i+1);
				$H[($m-1)*$m+$i] -= $h**($m-$i); 	
			}
			$H[($m-1)*$m] += ((2*$h-1 > 0) ? (2*$h-1)**$m : 0);
	
			for(my $i = 0; $i < $m; $i++){
				for(my $j = 0; $j < $m; $j++){
					if($i-$j+1 > 0){
						for(my $g=1; $g<=($i-$j+1); $g++ ){
							$H[$i*$m+$j] /= $g;
						}
					}
				}
			}
	
			my $eH = 0;
			$self->mPower(\@H, $eH, \@Q, \$eQ, $m, $num_values);
			$s = $Q[($k-1)*$m+$k-1];
			for(my $i=1; $i<=$num_values;$i++){
				$s = $s*$i/$num_values;
				if($s<1e-140){
					$s*=1e140;
					$eQ-=140;
				}
			}
	
			$s*=10**$eQ;
			$probability_value = 1-$s;
			
			$GU->user_info(3, 'P value: '.$probability_value.'\n');
			
			return $probability_value;
		}
			
	} #one_sided_ks_test


	method mMultiply(ArrayRef $A, ArrayRef $B, ArrayRef $Result, Int $m) {
		
		#Helper function for K-S test
	
		for(my $i=0;$i<$m;$i++){
			for(my $j=0;$j<$m;$j++){
				my $s = 0.0;
				for(my $k=0;$k<$m;$k++){
					$s+=$$A[$i*$m+$k]*$$B[$k*$m+$j];
				}
				$$Result[$i*$m+$j] = $s;	
			}
		}	
	} #mMultiply

	method mPower(ArrayRef $A, Int $eA, ArrayRef $Result, ScalarRef $eV, Int $m, Int $n) {
	
		#Helper function for K-S test
	
		my $eB;

		if($n==1){
			for(my $i=0;$i<($m*$m);$i++){
				$$Result[$i] = $$A[$i];	
			}	
			$$eV = int($eA);
			return;
		}
	
		#calls itself recursively
		$self->mPower($A, $eA, $Result, $eV, $m, int($n/2));
	
		my @B = ();
		$self->mMultiply($Result, $Result, \@B, $m);
		$eB = int(2 * $$eV);
	
		if($n%2 == 0){
			for(my $i=0; $i<($m*$m);$i++){
				$$Result[$i] = $B[$i];	
			}
			$$eV = int($eB);
		}else{
			$self->mMultiply($A, \@B, $Result, $m);
			$$eV = int($eB+$eA);
		}
		if($$Result[($m/2)*$m + ($m/2)] > 1e140){
			for(my $i=0;$i<($m*$m);$i++){
				$$Result[$i] = $$Result[$i] * 1e-140;	
			}
			$$eV += 140;	
		}			
	} # mPower
	
} # Statistics #
