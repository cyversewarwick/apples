#!/usr/bin/perl

use MooseX::Declare;

=head1 PSSM Set class

Stores many PSSMs

=cut

class Datatypes::Motifs::PSSM_Set extends Serialization::Serializable {
    
	use JSON;
	
    has 'pssms' => (
        is => 'rw',
        isa => 'ArrayRef[Datatypes::Motifs::PSSM]',
        required => 1,
    );
    
    require Datatypes::Motifs::PSSM;
    
=head2 return the result of merging this set with another one

 Parameters:
 $other : another sequence set
 
 Returns:
 A PSSM set containing the union of all the sequences in 
 this set.

=cut
	method merged_sets ( Datatypes::Motifs::PSSM_Set $other) {
		my %vals = ();
		
		foreach my $s (@{$self->{pssms}}) {
			my $js = to_json (Serialization::Serializable::to_hash($s), {canonical => 1});
			$vals{$js} = $s;
		}
		foreach my $s (@{$other->{pssms}}) {
			my $js = to_json (Serialization::Serializable::to_hash($s), {canonical => 1});
			$vals{$js} = $s;
		}
		
		my @vals = values %vals;
		
		return Datatypes::Motifs::PSSM_Set->new ( pssms => \@vals );
	}

    
}
