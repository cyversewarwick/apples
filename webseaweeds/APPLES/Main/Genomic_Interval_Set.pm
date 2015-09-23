### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genomic_Interval_Set class ###

use MooseX::Declare;

class Genomic_Interval_Set {
	use Bio::SeqIO;
	use Bio::Seq;
	use General_Utilities;
	use Genomic_Interval;
	use Genomic_Interval_Set_Maker;
	use Data::Dumper;
	use APPLES_Datatypes qw(NonNegativeInt WorkingSequenceType Boolean);
	use constant {FALSE => 0,
		      TRUE  => 1};

        use Moose;
	has 'genomic_interval_set' => (is => 'rw', isa => 'ArrayRef[Genomic_Interval]', required => 1);

	my $GU = General_Utilities->new();

	method gi_set_overlap (Genomic_Interval_Set $setb){
	    
	    my $overlap = FALSE;
	    
	    foreach my $gia(@{$self->genomic_interval_set}) {
		foreach my $gib(@{$setb->genomic_interval_set}) {
		    $overlap = $gia->gi_overlap($gib);
		    if ($overlap) {
			return $overlap;
		    }
		}
	    }
	    return $overlap;
	} # gi_set_overlap #
	
	method return_number_of_intervals {
	    my $result = @{$self->genomic_interval_set};

	    return $result;
	} # return_number_of_intervals #

	method return_unmasked_sequences {
	    # returns an array of the unmasked sequences

	    $GU->user_info(2,"Obtaining sequences for a GI set.\n");
	    my @sequences;
	    my $progress_count = 0;
	    my $total = @{$self->genomic_interval_set};
	    foreach my $gi ( @{$self->genomic_interval_set} ) {
		push (@sequences, $gi->get_sequence());
		$progress_count++;
		$GU->user_info(2,"Done ".$progress_count." out of ".$total.".\n");
	    }		
	    return @sequences;
	} # return_unmasked_sequences #
	
	method return_masked_sequences {
		die 'method not implemented yet!\n';
		# returns an array of the masked sequences
	} # return_masked_sequences #
	
	method total_length_of_gi_set_sequences {
	    my $total_length = 0;
	    foreach my $gi ( @{$self->genomic_interval_set} ) {
		$total_length = $total_length + $gi->get_length;
	    }
	    return $total_length;
	} # total_length_of_gi_set_sequences #
	
	method return_all_gi_lengths {
	    my @result;
	    foreach my $gi (@{$self->genomic_interval_set}) {
		my $length = $gi->get_length();
		push(@result,$length);
	    }
	    return @result;
	} # return_all_gi_lengths #

	method create_FASTA_from_genomic_interval_set(Str $output_fasta_filename) {
	    # uses the "working sequence" of its genomic intervals so it could be reference sequences,
	    # repeatmasked reference sequences or modified sequences
		
		#$output_fasta_filename = 'meme_input.fasta';
		my $output_stream = Bio::SeqIO->new(-file => '>'.$output_fasta_filename, -format => 'Fasta');
		
		foreach my $gi (@{$self->genomic_interval_set}){
			
			# Create new bioperl sequence object
			my $seq = Bio::Seq->new( -seq => $gi->get_working_sequence(),
			-id => $gi->label,
			);
			
			# Write bioperl sequence object to a .fasta file
			unless( $output_stream->write_seq($seq) ){
				die "There was a problem writing the fasta file\n";	
			}
			
		}		
		
		return $output_fasta_filename;
		
	} # create_FASTA_from_genomic_interval_set #

	method add_gi_sets(ArrayRef[Genomic_Interval_Set] $gi_sets) {
	  # adds given gi sets into $self by concatenating the arrays, overlaps of genomic intervals are not removed
	  my @all_gis = @{$self->genomic_interval_set};
	  foreach my $gi_set (@{$gi_sets}) {
	    push(@all_gis,@{$gi_set->genomic_interval_set});
	  }
	  $self->genomic_interval_set(\@all_gis);
	} # add_gi_sets #

	method merge(Boolean $turn_to_positive_strand) {
	    # produces a new Genomic_Interval_Set that covers the same sequence and in which
            # genomic intervals do not overlap
	    # all genomic intervals will be set to positive strand if $turn_to_positive_strand
		# modified sequences are dropped
		
	    my @intervals = @{$self->genomic_interval_set};
	    my $int = @intervals;
            # for each gi, see if it overlaps with any others; if it does, merge the intervals (replacing one element
	    # with the merged interval, and undefining the other);
            # once there are no more overlaps, collect remaining genomic intervals
	    $GU->user_info (2, $int." intervals to merge.\n");
	    my $mergers;
	    my $counter = 0;
	    do { 	    # do loop repeats until in the last run no mergers have occurred (so the remaining set must be overlap-free)
		$mergers = FALSE;
		for ( my $i = 0; $i < scalar (@intervals) ; $i++ ) {
		    # ignore undefined elements in the array
		    if (defined $intervals[$i]){
			for (my $j = $i+1; $j < scalar (@intervals) ; $j++ ) {
			    if (defined $intervals[$j]) {
				if ($intervals[$i]->gi_overlap($intervals[$j])) {
				    # merge the intervals together
				    $GU->user_info(3,"intervals overlap!\n");				    
				    my $merged_gi = $intervals[$i]->merge_interval($intervals[$j]);
				    $intervals[$i] = $merged_gi;
				    undef $intervals[$j];
				    $mergers = TRUE; # set merger flag to true - so do loop will repeat
				}
			    }
			} # end of inner for loop
		    } 
		} # end of outer for loop
	    } until (!$mergers); # end of do loop - if no mergers occurred in last iteration of do loop we will exit the loop
	    my @genomic_intervals;
	    $int = @intervals;
	    $GU->user_info(2, $int." intervals remain after merge.\n");
	    foreach my $gi (@intervals) {
		if (defined $gi) {
			if ($turn_to_positive_strand) {
		    my $five_prime_pos = $gi->five_prime_pos;
		    my $three_prime_pos = $gi->three_prime_pos;
		    if ($gi->three_prime_pos < $gi->five_prime_pos) {
			$five_prime_pos = $gi->three_prime_pos;
			$three_prime_pos = $gi->five_prime_pos;
		    }
			my $working_sequence = $gi->working_sequence;
		    my $genomic_interval = Genomic_Interval->new(
			'genome_db' => $gi->genome_db,
			'region' => $gi->region,
			'five_prime_pos' => $five_prime_pos,
			'three_prime_pos' => $three_prime_pos,
			'strand' => 'positive',
			'working_sequence' => $working_sequence
			);
		    if (defined $gi->coord_sys_name) {
			$genomic_interval->coord_sys_name($gi->coord_sys_name);
		    }
		    if (defined $gi->type) {
			$genomic_interval->type($gi->type);
		    }
		    
		    if(defined $gi->label) {
		    	$genomic_interval->label($gi->label);
		    }
		    
		    push (@genomic_intervals, $genomic_interval);
			} else {
			push (@genomic_intervals, $gi);
			}
		}
	    }
	    my $genomic_interval_set_maker = Genomic_Interval_Set_Maker->new();
	    my $gi_set = $genomic_interval_set_maker->make_gi_set_from_gi_array(\@genomic_intervals);
	    return $gi_set;
	} # merge #

	method first_n (Int $n) {
	    # returns a new genomic interval set containing the first n intervals
	    # (or less if there are not n intervals)

	    my @first_n_gis;
	    if ($n>0) {
		my $last = $n-1;
		@first_n_gis = @{$self->genomic_interval_set}[0..$last];
	    }
	    my $result = Genomic_Interval_Set->new(genomic_interval_set => \@first_n_gis);
	    return $result;
	} # first_n #

	method last_n (Int $n) {
	    # returns a new genomic interval set containing only the last n intervals
	    # (or less if there are not n intervals)

	    my @last_n_gis;
	    if ($n>0) {
		my $number_of_intervals = @{$self->genomic_interval_set};
		my $first = $number_of_intervals-$n;
		my $last = $number_of_intervals-1;
		if ($first<0) {
		    $first = 0;
		}
		@last_n_gis = @{$self->genomic_interval_set}[$first..$last];
	    }
	    my $result = Genomic_Interval_Set->new(genomic_interval_set => \@last_n_gis);
	    return $result;
	} # last_n #
	
	method all_intervals_are_on_same_genome_database() {
	    my $GDBU = Genome_DB_Utilities->new();
	    my $result = TRUE;
	    my $number_of_gis = $self->return_number_of_intervals();
	    if ($number_of_gis > 0) {
		my $first_gi = ${$self->genomic_interval_set}[0];
		my $first_genome_db_id = $GDBU->get_ID_for_database($first_gi->genome_db);
		foreach (my $index = 0; $index < $number_of_gis; $index++) {
		    my $current_gi = ${$self->genomic_interval_set}[$index];
		    my $next_db_ID = $GDBU->get_ID_for_database($current_gi->genome_db);
		    if ($next_db_ID ne $first_genome_db_id) {
			$result = FALSE;
		    }
		}
	    }
	    return $result;
	} # all_intervals_are_on_same_genome_database #

	method set_working_sequence(WorkingSequenceType $working_sequence_type) {
	    # sets the working sequence type for all its genomic intervals

	    my $number_of_gis = @{$self->genomic_interval_set};
	    for (my $i=0;$i<$number_of_gis;$i++) {
		${$self->genomic_interval_set}[$i]->working_sequence($working_sequence_type);
	    }
	} # set_working_sequence #

} # Genomic_Interval_Set #
