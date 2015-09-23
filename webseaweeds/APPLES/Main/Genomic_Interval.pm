### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genomic_Interval class ###
use MooseX::Declare;

class Genomic_Interval {
	use APPLES_Datatypes qw (StrandType SequenceChemicalType Boolean WorkingSequenceType PositiveInt);
	use constant {FALSE => 0,
		      TRUE  => 1};
	use Parameters;
	use Data::Dumper;
	use Genome_DB_Utilities;
	use General_Utilities;
	
	my $GU = General_Utilities->new();
	my $genome_db_utility = Genome_DB_Utilities->new();
	
	use Moose;
	has 'genome_db'	=>	(is => 'ro', isa => 'Genome_Sequence_Database_Parameters', required => 1, trigger => \&private_clear_sequences);
	has 'coord_sys_name' => (is => 'rw', isa => 'Str', trigger => \&private_clear_sequences); #
	has 'region'	=>	(is => 'ro', isa => 'Str', required => 1, trigger => \&private_clear_sequences); # must have a region identifier
	has 'five_prime_pos'	=>	(is => 'rw', isa => 'Int', required => 1, trigger => \&private_clear_sequences); # first character belonging to interval, positions to start from 1
	has 'three_prime_pos'	=>	(is => 'rw', isa => 'Int', required => 1, trigger => \&private_clear_sequences); # last character belonging to interval, positions to start from 1
	has 'strand'	=> (is => 'rw', isa => StrandType, required => 1, trigger => \&private_clear_sequences);
	has 'gi_sequence'  =>	(is => 'rw', isa => 'Str', clearer => 'private_clear_gi_sequence', reader => 'get_sequence'); # product of get_sequence method
	has 'gi_sequence_repeatmasked'  =>	(is => 'rw', isa => 'Str', clearer => 'private_clear_gi_sequence_repeatmasked', reader => 'get_repeatmasked_sequence'); # product of get_repeatmasked_sequence method
	has 'modified_sequence'  =>	(is => 'rw', isa => 'Str', clearer => 'private_clear_modified_sequence', reader => 'get_modified_sequence'); # product of get_modified_sequence method
	has 'working_sequence' => (is => 'rw', isa => WorkingSequenceType, required => 1, default => 'ref_sequence'); 
	has 'type' => (is => 'rw', isa => SequenceChemicalType); 
	has 'label' => (is=> 'rw', isa => 'Str'); # agi id e.g. at5g35609
	has 'sub_interval_annotation' => (is => 'rw', isa => 'ArrayRef', required => 0, reader => 'private_read_sub_interval_annotation', writer => 'private_write_sub_interval_annotation'); # can be used to record intervals of ORFs, protein domains, etc. and associated information; content of ArrayRef only to be exposed via public methods
	
	method get_working_sequence(){

		# Returns the working sequence as a string
		if ( $self->working_sequence() eq 'ref_sequence' ){
			return $self->get_sequence();
		}
		elsif( $self->working_sequence() eq 'ref_sequence_repeat_masked' ){
			return $self->get_repeatmasked_sequence();
		}
		elsif( $self->working_sequence() eq 'modified_sequence' ){
			return $self->get_modified_sequence();
		}
		else{
			die "Working sequence type is not recognised by this method. Debug me please!";	
		}
		
		
		
	} # get_working_sequence #
	
	method mask_genomic_interval( Genomic_Interval $GI_to_mask){
		# Method updates the modified_sequence within this genomic interval ($self), by masking out the subsequence $GI_to_mask.
		# also sets working sequence flag to modified sequence
		
		# Check genomic intervals overlap, otherwise $GI_to_mask cannot be used to mask this genomic interval
		# Assuming the GI we want to mask is actually a subsequence of $self:
		
		# Work out the position in the sequence we want to mask relative to the first character in the sequence and not
		# genomic coordinates. i.e the first position in the sequence is 0
		
		my $start_mask_position;
		my $end_mask_position;
		my $normalised_start;
		my $normalised_end;
		my $mod_sequence;
		
		# If modified_sequence() is already defined for this object then simply modify this, otherwise
		# we need to create a new modifed sequence.
		if( defined( $self->modified_sequence() ) ){
			$mod_sequence = $self->modified_sequence();
		}
		else{
			$mod_sequence = $self->gi_sequence();
		}
		
		# Check for overlap
		if ( $self->gi_overlap($GI_to_mask) ){
		
			if ( $self->strand eq 'positive' ){
				
				if ( $GI_to_mask->strand() eq 'positive' ){
					$start_mask_position =  $GI_to_mask->five_prime_pos;
					$end_mask_position = $GI_to_mask->three_prime_pos;
		
				}
				else{
					$start_mask_position =  $GI_to_mask->three_prime_pos;
					$end_mask_position = $GI_to_mask->five_prime_pos;
					
				}
				
				$normalised_start = $start_mask_position - $self->five_prime_pos;
				$normalised_end = $normalised_start + ( $end_mask_position - $start_mask_position );
				
			}
			else{
				
				if ( $GI_to_mask->strand() eq 'negative' ){
					$start_mask_position =  $GI_to_mask->three_prime_pos;
					$end_mask_position = $GI_to_mask->five_prime_pos;
				}
				else{
					$start_mask_position =  $GI_to_mask->five_prime_pos;
					$end_mask_position = $GI_to_mask->three_prime_pos;
				}
				
				$normalised_start = $start_mask_position - $self->three_prime_pos;
				$normalised_end = $normalised_start + ( $start_mask_position - $end_mask_position );
			}
			
			# Deal with gi's that we want to mask, but are only partial subsequences of $self; i.e. they 'hang over the edge' of $self->sequence.
			if( $normalised_start < 0 ){
				$normalised_start = 0;
			}
			elsif( $normalised_end > length($mod_sequence) ){
				$normalised_end = length($mod_sequence) - 1;
			}
			
			# Create the string with which we will 'mask' the subsequence
			my $mask = "N";
			for ( 1 .. ( $normalised_end - $normalised_start ) ){
				$mask .= 'N';
			}
			
			# Substitute sequence with string of N's
			substr( $mod_sequence, $normalised_start, ($normalised_end - $normalised_start + 1), $mask);
	
			$self->modified_sequence( $mod_sequence );
			$self->working_sequence('modified_sequence');
		}
		
	} # mask_genomic_interval #
		
	method get_working_sequence_type(){
			
		if (!(defined $self->working_sequence)){
				
			die "No working sequence has been defined within this object\n";
		}
		
		return $self->working_sequence();
	} # get_working_sequence_type #
	
	method get_modified_sequence(){
		
		# returns the modified sequence attribute
		
		if( defined( $self->modified_sequence() ) ){
			return $self->modified_sequence();
		}
		else{
				die "Modified sequence is not defined\n";
		}
	} # get_modified_sequence #
	
	method get_sequence() {
		# retrieves sequence from specified database
		
	    if (!(defined $self->gi_sequence)) {
		my $sequence = $self->private_get_genomic_interval_sequence();
		if (length($sequence)<500000) {
		    $self->gi_sequence( $sequence );
		}
		else {
		  $GU->user_info(3,"Did not store a large unmasked sequence within Genome_Interval on region ".$self->region.".\n");
		  return $sequence;
		}
	    }
	    return $self->gi_sequence;
	} # get_sequence #

	method get_length() {
	    my $result = $self->five_prime_pos-$self->three_prime_pos;
	    if ($result<0) {
		$result = -$result;
	    }
	    $result++;
	    return $result;
	} # get_length #

	method get_repeatmasked_sequence() {
	    if (!(defined $self->gi_sequence_repeatmasked)) {
		my $repeatmasked_sequence = $self->private_get_genomic_interval_repeatmasked_sequence();
		if (length($repeatmasked_sequence)<500000) {
		    $self->gi_sequence_repeatmasked( $repeatmasked_sequence ); 
		}
		else {
            $GU->user_info(3,"Did not store a large repeatmasked sequence within Genome_Interval on region ".$self->region.".\n");
		    return $repeatmasked_sequence;
		}
	    }
	    return $self->gi_sequence_repeatmasked;
	} # get_repeatmasked_sequence #

	method genomic_interval_details {
	
		die 'method not implemented yet!';
	
		return $self->species . ' ' . $self->assembly . ' ' . $self->region . ':' . 
		$self->start_pos . '-' . $self->end_pos;
	} # genomic_interval_details #
	
	method gi_overlap (Genomic_Interval $target_gi) {
		#$GU->user_info( 3, $self->genome_db->scientific_species_name."\n" );
		#$GU->user_info( 3, $target_gi->genome_db->scientific_species_name."\n" );
		# input another genomic interval
		# determine if the two intervals overlap
		# returns Boolean

		my $startA; #smaller of $self->five_prime_pos and $self->three_prime_pos
		my $endA; # greater of $self->five_prime_pos and $self->three_prime_pos
		my $startB; # smaller of $target_gi->five_prime_pos and $target_gi->three_prime_pos
		my $endB; # greater of $target_gi->five_prime_pos and $target_gi->three_prime_pos
		
		my $overlap = FALSE;
		my $result;
		
		my $one_ID = $genome_db_utility->get_ID_for_database($self->genome_db);
		my $other_ID = $genome_db_utility->get_ID_for_database($target_gi->genome_db);

		if ($one_ID ne $other_ID) {
		  return FALSE;
		}
		
		if ((!defined $self->coord_sys_name)||
		    (!defined $target_gi->coord_sys_name)) {
		  if ((defined $self->coord_sys_name)||
		      (defined $target_gi->coord_sys_name)) {
		    die 'Either both or none of the database must have the coordinate system defined.';
		  }
		}
		if (defined $self->coord_sys_name) {
		  if ($self->coord_sys_name ne $target_gi->coord_sys_name) {
		    return FALSE;
		  }
		}

		if ($self->region eq $target_gi->region) {
			
			if ($self->five_prime_pos > $self->three_prime_pos) {
				$startA = $self->three_prime_pos;
				$endA = $self->five_prime_pos;
			}
			else {
				$startA = $self->five_prime_pos;
				$endA = $self->three_prime_pos;
			}
		
			if ($target_gi->five_prime_pos > $target_gi->three_prime_pos) {
				$startB = $target_gi->three_prime_pos;
				$endB = $target_gi->five_prime_pos;
			}
			else {
				$startB = $target_gi->five_prime_pos;
				$endB = $target_gi->three_prime_pos;
			}
			
			if ( $startA < $startB ) {
				$result = $startB - $endA;
			}
			else {
				$result = $startA - $endB;
			}
			if ( $result < 0 ) {
				$overlap = TRUE;
			}
			
			if ( $startA == $endB ){
				$overlap = TRUE;
			}
			elsif ( $startB == $endA ){
				$overlap = TRUE;
			}
			
		}
		return $overlap;
	} # gi_overlap #

	method merge_interval (Genomic_Interval $interval_b) {
	  # takes another (overlapping) GI as input (you should only give an overlapping interval, else it will die, but we check again anyway..), and merges them together, returning the merged interval
      # the modified sequences are dropped
		
	    $GU->user_info(3,"merging two intervals\n");
	    #first a quick test to check they DO overlap: die if not
	    my $overlap = $self->gi_overlap($interval_b);
	    if (!$overlap) {
		die 'intervals do not overlap, they cannot be merged!';
	    }
	    
	    my @positions = ($self->five_prime_pos, $self->three_prime_pos, $interval_b->five_prime_pos, $interval_b->three_prime_pos);
	    @positions = sort { $a <=> $b } @positions;
	    my $min = $positions[0];
	    my $max = $positions[-1];
	    
		my $working_sequence = 'ref_sequence';
		if ($self->working_sequence eq 'ref_sequence_repeat_masked') {
			if ($interval_b->working_sequence eq 'ref_sequence_repeat_masked') {
				$working_sequence = 'ref_sequence_repeat_masked';
			}
		}
		
		my $five_prime_pos = $min;
		my $three_prime_pos = $max;
		if ($self->strand eq 'negative') {
			$five_prime_pos = $max;
			$three_prime_pos = $min;
		}
		
	    my $merged_interval = Genomic_Interval->new(
			'genome_db' => $self->genome_db,
			'region' => $self->region,
			'five_prime_pos' => $five_prime_pos,
			'three_prime_pos' => $three_prime_pos,
			'strand' => $self->strand,
			'working_sequence' => $working_sequence,
		);
		
		if (defined $self->coord_sys_name) {
			$merged_interval->coord_sys_name($self->coord_sys_name);
		}
		if (defined $self->type) {
			$merged_interval->type($self->type);
		}
		
		if(defined $self->label) {
			$merged_interval->label($self->label);
		}
	    
	    return $merged_interval;
	} # merge_interval #
	
	method add_primary_ORF() {
	    # will add likely primary ORF to sub_intervals if ORF can be found

# NEXT STEP: isolate code from CTIP that calls get_multiple_start_intervals() and this integrate here, together with the code below

	    die 'method not implemented yet!';

	    $self->private_require_DNA();
	    my $dna_sequence = UC($self->get_sequence());

 # these code reading frames containing one or more ATGs
#	my $readingframe = shift;
	
#	my @result;
#	my $startcodon = 'ATG';
#	my @stopcodons = ('TAG','TGA','TAA');
#	my $withinframe = 'false';
#	my @currentstarts;
	
	
	
	# loop through every third base in the sequence starting from each reading frame
#	for (my $i=$readingframe;$i<length($sequencestring)-2;$i=$i+3) {
#		my $codon = uc(substr($sequencestring,$i,3));
#		if ($codon eq $startcodon) {
#			$withinframe = 'true';
#			push(@currentstarts,$i); #push the position of the in-frame starts
#		}
#		if ($withinframe eq 'true') {
#			my $isstopcodon = 'false';
#			foreach my $stopcodon (@stopcodons) {
#				if ($codon eq $stopcodon) {
#					$isstopcodon = 'true';
#				}
#			}
	#		if ($isstopcodon eq 'true') {
	#			$withinframe = 'false';
		#		my $rec = {
	#				STOPPOSITION => $i,
	#			};
#				@{$rec->{STARTPOSITIONS}} = @currentstarts;
#				@currentstarts = ();
#				push(@result,$rec);
#			}
#		}
#	}
#	return @result;
	# result contains an array of references to hash table records $rec. Each element in te array is a reference
	# to a different record for a start...stop interval in the same frame

	   
	} # add_primary_ORF #

	method add_secondary_ORF() {
	    # will add one or more secondary ORFs to sub_intervals if these can be found

	    die 'method not implemented yet!';
	    $self->private_require_DNA();
	} # add_secondary_ORF #

	method get_translated_ORFs() {

	    die 'method not implemented yet!';
	    $self->private_require_DNA();
	} # get_translated_ORFs #

	method expand_interval(PositiveInt $five_prime_flank_length, PositiveInt $three_prime_flank_length) {
		# lengthens interval on both sides by given offsets
		
		my $new_five_prime_pos;
		my $new_three_prime_pos;
		
		if ( $self->strand eq 'positive' ){
			$new_five_prime_pos = $self->five_prime_pos - $five_prime_flank_length;
			$new_three_prime_pos = $self->three_prime_pos + $three_prime_flank_length;
		}
		elsif ( $self->strand eq 'negative' ){
			$new_five_prime_pos = $self->five_prime_pos + $five_prime_flank_length;
			$new_three_prime_pos = $self->three_prime_pos - $three_prime_flank_length;
		}
		else {
			die 'strand not well-defined.';
		}
		$self->five_prime_pos($new_five_prime_pos);
		$self->three_prime_pos($new_three_prime_pos);
	} # expand_interval #
	
	method private_read_sub_interval_annotation() {
	    return $self->sub_interval_annotation;
	} # private_read_sub_interval_annotation #

	method private_write_sub_interval_annotation(ArrayRef $sub_interval_annotation) {
	    $self->sub_interval_annotation = $sub_interval_annotation;
	} # private_write_sub_interval_annotation #

	method private_require_DNA() {
	    if ($self->type ne 'dna') {
		die 'genomic interval is not DNA, but DNA is required.';
	    }
	} # private_require_DNA #

	sub private_clear_sequences() {
	    my ($self, $arg) = @_; 
	    $self->private_clear_gi_sequence(); # clearer methods are provided by Moose if a name is provided (as done above)
	    $self->private_clear_gi_sequence_repeatmasked();
		$self->private_clear_modified_sequence();
	} # private_clear_sequences #

	method private_get_genomic_interval_sequence {
	    my $result;
	    if (defined $self->coord_sys_name) {		
		$result = $genome_db_utility->get_genomic_unmasked_sequence_of_gi($self->genome_db, $self->coord_sys_name, $self->region, $self->five_prime_pos, $self->three_prime_pos, $self->strand);
	    } else {
		$result = $genome_db_utility->get_genomic_unmasked_sequence_of_gi($self->genome_db, 'DUMMY', $self->region, $self->five_prime_pos, $self->three_prime_pos, $self->strand);
	    }
	    return $result;
	} # private_get_genomic_interval_sequence #

	method private_get_genomic_interval_repeatmasked_sequence {
	    
	    my $result = $genome_db_utility->get_genomic_repeatmasked_sequence_of_gi($self->genome_db, $self->coord_sys_name, $self->region, $self->five_prime_pos, $self->three_prime_pos, $self->strand);
	    
	    return $result;
	} # private_get_genomic_interval_repeatmasked_sequence #
	
} # Genomic_Interval #
