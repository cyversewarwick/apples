#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Subtasks::Mark_Conservation extends Jobs::Job {
	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	has 'pval_threshold' => (
		is => 'rw',
		isa => 'Num',
		default => 1e-3,
		documentation => "1. Pvalue threshold"
	);
	
	has 'masks' => (
		is => 'rw',
		isa => 'HashRef[Int]',
		default => "",
		documentation => "2. Masking values"
	);
	
	has 'plot' => (
		is => 'rw',
		isa => 'Datatypes::Results::Alignment_Plot',
	);
	

=head2 Overloaded validation method

=cut

	method validate () {
		if ( scalar @{ $self->selected_dataitems } < 1 ) {
			die "You need to select at least one data item.";
		}
		
		my $hr = from_json ($self->masks);
		if (!Serialization::Serializable::is_hash ($hr)) {
			die "Invalid masks string.";
		}
		$self->SUPER::validate();
	}

=head2 Make masking array for a sequence

 Parameters:
 $seq : a genomic sequence
 $window : a window length
 
 Returns:
 an array of size sequence length - window length containing 
 [
 	{
 		ann_type => fraction (0-1),
 		...
 	}, 
 	
 	...
 ]

=cut
	method get_window_counts ( Sequences::Genomic_Sequence $seq, Int $window ) {
		my $x = 0;
		my @arr = ();
		my $sequence = $seq->seq; 
		for ($x = 0; $x < length ( $sequence ) - $window + 1 ; ++$x ) {
			push @arr, { };
		}

		my $sequence_start = 
		foreach my $ann (@{ $seq->{annotation} }) {
			my $type = $ann->{type};
			## dust is a bit different from proper repeats
			if ($type eq "repeat" && $ann->{sourceid} eq "dust") {
				$type = "dust";
			}
			## so are trfs
			if ($type eq "repeat" && $ann->{sourceid} eq "trf") {
				$type = "trf";
			}
			
			my 
					
		}

		my $window_pos = 1 - $window
		for ($x = 0; $x < length ( $sequence ) ; ++$x ) {
			
			++$window_pos;
		}
		return \@arr;
	}

=head2 Run analysis job

 Parameters:
 None

 Returns:
 A serializable array of Genomic_Sequence objects

=cut
	method _run () {
		my $plot = $self->plot;
		my $s1 = $plot->sequence_1;
		my $s2 = $plot->sequence_2;
		my $pdata = $plot->plot;
		
		## create masked profiles
		foreach my $xy ( @$pdata ) {
			my @coords = split /\_/, $xy;
			if (scalar @coords != 2) {
				die "Invalid alignment plot result.";
			}
			
		}
	}
	
};