#!/usr/bin/perl 

use MooseX::Declare;

=head1 Generic Overrepresentation test template

An overrepresentation test takes as its input a list of p-values, 
and a threshold for which to accept/reject.

=cut

class Statistics::Overrepresentation_Test {
	
	## A set of p-values
	has 'pvalues' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
	);
	
	## A threshold value
	has 'threshold' => (
		is => 'rw',
		isa => 'Num',
	);
	
=head2 Run the test

 Return a list of array indices of pvalues to keep like this: 

 {
	pvalue => threshold,
	indexes => \@a,
 }

 By default, we return everything that is smaller than the threshold.

=cut
	
	method test () {
		my @a = ();
		my $i = 0;
		foreach my $val (@{$self->pvalues}) {
			push @a, $i
				if $val < $self->threshold;
			$i++;
		}
		
		return {
			pvalue => $self->threshold,
			indexes => \@a,
		};
	}
}
