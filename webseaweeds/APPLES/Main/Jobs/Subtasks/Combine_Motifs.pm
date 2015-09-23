#!/usr/bin/perl

=head1 Combine the motifs present in multiple sequences according to 
a variety of rules.

=cut

use MooseX::Declare;

class Jobs::Subtasks::Combine_Motifs extends Jobs::Job {
	use Runtime;

	use Datatypes::Moose_Types;
	
	## The master sequence to start with
	has 'master_sequence' => (
		is => 'rw',
		isa => 'Sequences::Genomic_Sequence',
		required => 1,
		
	);
	
	## Related sequences
	has 'related_sequences' => (
		is => 'rw',
		isa => 'Datatypes::Sequences::Set',
		required => 1,
		documentation => "The sequences to work on."
	);
	
	## 	in a fraction of how many sequences does a given motif have to be present?
	##  ranges from 0 (set union) to 1 (set intersection)
	has 'combine_fraction' => (
		is => 'rw',
		isa => 'Num',
		required => 1,
		default => sub { return 1; },
		documentation => "Between 0 and 1, specify the fraction of the input sets in which the motif must be present to be retained."
	);

	## what is the type of data to work on 
	has 'ann_type' => (
		is => 'rw',
		isa => 'Num',
		required => 1,
		default => sub { return "pssmhit"; },
		documentation => "Which data to work on (use pssmhit or bifahit)."
	);

=head2 Run analysis job

 Parameters:
 None

 Returns:
 A Genomic_Sequence object based on master_sequence, but with filtered annotations

=cut

	method _run () {
        my @sequences = @{$self->related_sequences->sequences};

		my $sq = $self->master_sequence->cloneme();		
		my %anns = ();
		
		## delete all corresponding anns from copy
		$sq->delete_annotation(sub {
			my $ann = shift;
			if ($ann->{type} eq $self->ann_type) {
				return 0;
			} else {
				return 1;
			}
		});

		## make a list of the anns we deleted
		my $xanns = $self->master_sequence->get_annotations_by_type( $self->ann_type );
		foreach my $ann (@{$xanns}) {
			if (defined ($anns{$ann->{sourceid}})) {
				push @{$anns{$ann->{sourceid}}->{anns}}, $ann;
			} else {
				$anns{$ann->{sourceid}} = {
					anns => [ $ann ],
					inseqs => { $self->master_sequence->id() => 1 },
				};
			}
		}
		
		## count occurrences in other sequences
		my %allids = ();
		$allids{$self->master_sequence->id()} = 1;
		foreach my $seq (@sequences) {
			$allids{$seq->id()} = 1;
			my $xanns = $seq->get_annotations_by_type( $self->ann_type );
			foreach my $ann (@{$xanns}) {
				if (defined ($anns{$ann->{sourceid}})) {
					$anns{$ann->{sourceid}}->{inseqs}->{$seq->id()} = 1;
				}
			}
		}
		
		## add the ones which are above the threshold
		my $allcount = 1 + scalar (keys %allids);
		my $threshold_count = $allcount * $self->combine_fraction;
		
		while ( my ($aid, $ad) = each (%anns) ) {
			my $count = scalar (keys %{ $ad->{inseqs} });
			if ($count > $threshold_count) {
				foreach my $ann (@{$ad->{anns}}) {
					$sq->add_annotation($ann);
				}
			}
		}

		return $sq;
	}
}
