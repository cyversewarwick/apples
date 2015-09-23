#!/usr/bin/perl 

use MooseX::Declare;

=head1 Binomial overrepresentation test

We use R to compute binomial distribution values and test
for overrepresentation within in a given sample size.

=cut

class Statistics::Overrepresentation_Test::Binomial extends Statistics::Overrepresentation_Test {

	use Runtime;
	
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

		eval_R(<<END
pvals  = sort(scan('$filename'))
bpvals = array(1, length(pvals))
for (i in 1:length(pvals)) bpvals[i] = pbinom(i, $overall_samples, pvals[i], lower.tail=FALSE)
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
		debug ("min-count: $out[0], min-pval: $out[1]");
		
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