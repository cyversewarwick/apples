#!/usr/bin/perl 

use MooseX::Declare;

=head1 Poisson-Binomial overrepresentation test

We use R to compute poisson-binomial distribution values and test
for overrepresentation within in a given sample size.

This is more accurate, but slower than the standard binomial 
test.

We require the R package "poibin", which must be installed in R.

=cut

class Statistics::Overrepresentation_Test::PoissonBinomial extends Statistics::Overrepresentation_Test {

	use Runtime;
	use Data::Dumper;
	use File::Temp qw(tempfile);

=head2 Run the test

 Parameters:
 $overall_samples : number of samples taken overall (N for the binomial test)

 Return a list of array indices of pvalues to keep.

 By default, we return everything that is smaller than the threshold.

=cut
	method test (Int $overall_samples) {		
		my ($fh, $filename) = tempfile();
		
		print $fh "$_\n" foreach @{$self->pvalues};
		close $fh;

        
        #Compute poisson probability for each number of hits, choose best one
		eval_R(<<END
library(poibin)
pvals  = sort(scan('$filename'))
bpvals = array(1, length(pvals))
for (i in 1:length(pvals)) bpvals[i] = 1-ppoibin(i, c(pvals[1:i], replicate($overall_samples-i, pvals[i])), method="RF")
m = which.min(bpvals)
xv = c(m, bpvals[m])
write(xv, '$filename', sep="\\n")
END
);
		my @out = ();
		{
			open $fh, "<", $filename;
			local $/ = undef;
			@out = split /\n/, <$fh>;
			close $fh;
		}
        

		if ($out[1] < $self->threshold) {		
			my @op;
			my $i = 0;
			push @op, $i++ foreach @{$self->pvalues};
		
			@op = sort {$self->pvalues->[$a] <=> $self->pvalues->[$b]} @op;
		
			my @ret = ();
			for ($i = 0; $i < $out[0]; ++$i) {
				push @ret, $op[$i];
			}
			return {
				pvalue => $out[1],
				indexes => \@ret,
			};
		} else {
			return {
				pvalue => $out[1],
				indexes => [],
			};			
		}
	}
	
}
