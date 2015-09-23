#!/usr/bin/perl

use MooseX::Declare;

=head1 A set of Genomic Sequences (and locations)

In this class, we manage creation and storage of multiple sequences,
potentially from different databases.

=cut

class Datatypes::Sequence::Set extends Datatypes::Datatype {

	use Runtime;
	use JSON;
	use Serialization::Serializable;

## the sequences
	has 'sequences' => (
		is  => 'rw',
		isa => "ArrayRef[Sequences::Genomic_Sequence]",
		default => sub { return [] },
	);

=head2 return the result of merging this set with another one into $self

 Parameters:
 $other : another sequence set
 
 Returns:
 A sequence set containing the union of all the sequences in 
 this set.

=cut
	method merged_sets ( Datatypes::Sequence::Set $other) {
		my %vals = ();
		
		foreach my $s (@{$self->{sequences}}) {
			my $js = to_json (Serialization::Serializable::to_hash($s), {canonical => 1});
			$vals{$js} = $s;
		}
		foreach my $s (@{$other->{sequences}}) {
			my $js = to_json (Serialization::Serializable::to_hash($s), {canonical => 1});
			$vals{$js} = $s;
		}
		
		my @vals = values %vals;
		
		my $s = Datatypes::Sequence::Set->new ( sequences => \@vals );
		$s->compact_sequences;
		return $s;
	}
	
=head2 Serialization override which compacts sequences first

=cut
	
	sub to_hash {
		my $self = shift;
		my @params = @_;
		$self->compact_sequences();
		$self->SUPER::to_hash(@params);
	}
	
=head2 Try to merge sequence objects if possible

Merge sequence objects when possible

=cut

	method compact_sequences () {
		my %seqs = ();
		
		foreach my $seq (@{$self->{sequences}}) {
			if (!defined ($seqs{$seq->ident})) {
				$seqs{$seq->ident} = $seq;
			} else {
				$seqs{$seq->ident}->merge($seq); 			
			}
		}
		
		my @seqns = values %seqs;
		
		$self->sequences(\@seqns);
	}
}
