#!/usr/bin/perl

=head1 List-related Utility functions

=cut

package Utilities::List;

use strict;

use constant {
			   FALSE => 0,
			   TRUE  => 1
};

=head2 Clear duplicates from a list 

Takes an array, removes any duplicate elements, returns array of
only the unique elements

Parameters:
 $list_ref : an ARRAYREF
 
Returns:
 $list_ref less duplicate values

=cut

sub remove_duplicates_from_list {
	my $list_ref = shift;

	my @unique_list;
	my %seen;
	foreach my $element (@$list_ref) {
		next if $seen{$element}++;
		push( @unique_list, $element );
	}
	return @unique_list;
}    # remove_duplicates_from_list #

=head2 Check if list has redundancy 

Takes an array, removes any duplicate elements, returns array of
only the unique elements

TODO: this implementation is really slow. fixme.

Parameters:
 $list_ref : an ARRAYREF
 
Returns:
true if there are redundant values

=cut

sub list_has_redundancy {
	my $list_ref             = shift;
	my $length               = @{$list_ref};
	my @non_redundant        = remove_duplicates_from_list($list_ref);
	my $result               = TRUE;
	if ( $length == scalar @non_redundant ) {
		$result = FALSE;
	}
	return $result;
}    # list_has_redundancy #

=head2 Check for list overlap 

Check if two arrays contain some of the same elements

Parameters:
 $list_ref1 : an ARRAYREF
 $list_ref2 : an ARRAYREF
 
Returns:
 TRUE if there are common values

=cut

sub lists_overlap {
	my $listref1 = shift;
	my $listref2 = shift;

	my $result = FALSE;
	foreach my $element1 ( @{$listref1} ) {
		foreach my $element2 ( @{$listref2} ) {
			if ( $element1 eq $element2 ) {
				return TRUE;
			}
		}
	}
	return FALSE;
}    # lists_overlap #

=head2 Get elements contained in both of two lists 

Check if two arrays contain some of the same elements and return these

Parameters:
 $reference_list : an ARRAYREF
 $subject_list : an ARRAYREF
 
Returns:
 an ARRAYREF containing common values

=cut

sub list_intersection {
	my $reference_list = shift;
	my $subject_list   = shift;

	my %hash;
	foreach my $element ( @{$reference_list} ) {
		$hash{$element} = 13;
	}
	my @result;
	foreach my $element ( @{$subject_list} ) {
		if ( defined $hash{$element} ) {
			push( @result, $element );
		}
	}
	return @result;
}    # list_intersection #

=head2 Get elements contained in both of two lists 

Return selected elements, given a boolean for each.

Parameters:
 $array_ref : an ARRAYREF
 $booleans : an ARRAYREF
 
Returns:
 an ARRAYREF containing common values

=cut

sub subset_array_by_booleans {
	my $array_ref = shift;
	my $booleans  = shift;

	# both arrays must be same length

	my @result;
	my $length       = @{$array_ref};
	my $other_length = @{$booleans};
	if ( $length != $other_length ) {
		die 'arrays must be of same length! (' 
		  . $length . ' vs '
		  . $other_length . ')';
	}
	for ( my $i = 0 ; $i < $length ; $i++ ) {
		if ( ${$booleans}[$i] ) {
			push( @result, ${$array_ref}[$i] );
		}
	}
	return @result;
}    # subset_array_by_booleans #

1;
