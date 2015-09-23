#!/usr/bin/perl 

use MooseX::Declare;

=head1 Q-Value multiple hypothesis test

We use meme's qvalue utility to compute q-values for each 
p-value, limiting the false-discovery rate.

=cut

class Statistics::Overrepresentation_Test::QValue {
	use File::Temp qw(tempfile);
	use IPC::Open3;
	use Symbol qw(gensym);
	use IO::File;
	use FindBin qw($Bin);
	
	use Runtime;
	use Configuration::AppleSeeds;
	
	## A set of p-values
	has 'pvalues' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
	);
	
	## q-value threshold
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
		my $qvalue = find_executable('qvalue');
		
		my ($fh, $filename) = tempfile();
		
		my @sp = ();
		my $i = 0;
		foreach (@{$self->pvalues}) {
			push @sp, $i++;
		}
		
		@sp = sort { $self->pvalues->[$a] <=> $self->pvalues->[$b] } @sp;
		
		my $nl_print = 0;
		for ($i = 0; $i < scalar @{$self->pvalues}; ++$i) {
			if ($nl_print == 1) {
				print $fh "\n";
			}
			print $fh (sprintf ("%.20f", $self->pvalues->[$sp[$i]]));
			$nl_print = 1;
		}
		close $fh;

		my @qvalues = ();
		local *CATCHERR = IO::File->new_tmpfile;
		my $pid = open3(gensym, \*CATCHOUT, ">&CATCHERR", "$qvalue --pi-zero $filename");
		while( <CATCHOUT> ) {
			push @qvalues, (split /\t/)[1];
		}
		waitpid($pid, 0);
		my $err = $?;
		if ($err) {
			seek CATCHERR, 0, 0;
			while( <CATCHERR> ) {
				error ($_);
			}
			die "Child process exited with $err ($qvalue --pi-zero $filename)"
		}
		
		if (scalar @qvalues != scalar @{$self->pvalues}) {
			die "Incorrect number of q-values returned.";
		}
		
		my @rx = ();

		my @sorted_pvalues = sort {$a <=> $b} @{$self->pvalues};
		# debug ("Have q-values @qvalues for @sorted_pvalues");

		## put qvalues into right order.
		my @aqvalues;
		$#aqvalues = $#qvalues;
		for ($i = 0; $i < scalar @{$self->pvalues}; ++$i) {
			my $qq = $qvalues[$i];
			$aqvalues[$sp[$i]] = $qq;
			
			if ( $qq < $self->threshold ) {
				push @rx, $i;
			}
		}

		return {
			pvalue => $self->threshold,
			indexes => \@rx,
			all_qvalues => \@aqvalues,
			all_pvalues => $self->pvalues,
		};
	}
}
