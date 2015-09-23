#!/usr/bin/perl

=head1 Set-related Utility functions

=cut

package Utilities::Set;

use strict;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(random_subset);


use Math::Random::MT::Auto 'rand';

use constant {
			   FALSE => 0,
			   TRUE  => 1
};

=head1 Pick a random subset from an array

 Parameters:
 $set : an arrayref
 $n   : the number of elements to pick out

 Returns:
 an ARRAYREF with a random selection from $set

=cut

sub random_subset {
    my $set = shift;
    my $n = shift || 1;

    my %picked = ();
    my $size = scalar @$set;

    while ($n > 0) {
        my $idx = int ( rand() * $size );
        my $v = $set->[$idx];

        if (!exists $picked{$v}) {
            $picked{$v} = 1;
            --$n;
        }
    }
    return keys %picked
        if wantarray;
    my @arr = keys %picked;
    return \@arr;
}


1;