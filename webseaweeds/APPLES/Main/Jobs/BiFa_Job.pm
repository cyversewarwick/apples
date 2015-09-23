#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::BiFa_Job extends Jobs::Job {

	with "Jobs::Roles::Dataitem_Selection_Job";

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	use Jobs::Subtasks::BiFa;
	use Serialization::Serializable_Array;
	use Datatypes::Sequence::Set;
	use Datatypes::Motifs::Motif_Hit;
    use Data::Dumper;

	## the BiFa job has two operating modes
	## 1. we can assume all sequences are related, in which case we analyse each sequence
	##    giving the other sequences as related.
	## 2. we can analyse each sequence individually
	has 'assume_phylo' => (
		is            => "rw",
		isa           => "Bool",
		default       => sub { return 1; },
		documentation => "1. Assume all sequences are phylogenetically related."
		  . "|When this is switched off, all sequences will be analysed individually.",
	);

	## BiFa consensus parameter
	has 'species' => (
		is            => "rw",
		isa           => "Str",
		default       => sub { return "vertebrates"; },
		documentation => "2. Select the species to use.|"
		  . "Possible values: plants, vertebrates, insects, nematodes, fungi, bacteria, or any combination of these, separated by spaces).",
	);

	## BiFa consensus parameter
	has 'use_consensus' => (
		is      => "rw",
		isa     => "Bool",
		default => sub { return 1; },
		documentation =>
		  "3. Use consensus sequences as well as score matrices.",
	);

	## BiFa consensus parameter
	has 'pssm_filter' => (
		is            => "rw",
		isa           => "Str",
		default       => sub { return "."; },
		documentation => "4. PSSM Regex filter",
	);

	## BiFa consensus parameter
	has 'non_transfac_pssms' => (
		is            => "rw",
		isa           => "Str",
		default       => sub { return ""; },
		documentation => "5. Non-transfac PSSM sets, separated by spaces.|"
		  . "Examples are: PLACE, jaspar, ...",
	);

	## threshold
	has 'threshold' => (
		is            => "rw",
		isa           => "Num",
		default       => sub { return 0.005; },
		documentation => "S1. BiFa's p-value threshold|Must be between 0.001 and 0.05",
	);

	## phylo threshold
	has 'phylo_threshold' => (
		is            => "rw",
		isa           => "Num",
		default       => sub { return 0.05; },
		documentation => "S2. BiFa's p-value threshold when using phylogeny|Must be between 0.001 and 0.1",
	);

	## this is where we write the output to, not cached, and only used
	## in postprocess
	has 'nc__output_dataitem' => (
		is            => "rw",
		isa           => "Str",
		default       => sub { return "BiFa analysis result"; },
		documentation => "O. The output dataitem to write to.",
	);

=head2 Overloaded validation method

=cut

	method validate () {
		if ( scalar @{ $self->selected_dataitems } < 1 ) {
			die "You need to select at least one data item.";
		}
		$self->SUPER::validate();
	}

=head2 Run BiFa analysis job

 Parameters:
 None

 Returns:
 A serializable array of Genomic_Sequence objects

=cut

	method _run () {
		my $result = [];

		my $species_string = $self->species();

		my @pssmsets = split /[\s]+/, $self->non_transfac_pssms;

		## transform string
		$species_string =~ s/plants/p/ig;
		$species_string =~ s/vertebrates/v/ig;
		$species_string =~ s/fungi/f/ig;
		$species_string =~ s/nematodes/n/ig;
		$species_string =~ s/insects/i/ig;
		$species_string =~ s/bacteria/b/ig;
		$species_string =~ s/\s//g;
		$species_string =~ s/[^pvinbf]//ig;

		$species_string = uc($species_string);

		my @sequences = $self->selected_objects_of_kind("Sequences::Genomic_Sequence");

		foreach my $seq (@sequences) {
			my $b;
			if ( !$self->assume_phylo() ) {
				## individual analysis.
				$b = Jobs::Subtasks::BiFa->new(
					sequence        => $seq->seq(),
					threshold       => $self->threshold(),
					phylo_threshold => $self->phylo_threshold(),
					species         => $species_string,
					pssm_filter     => $self->pssm_filter(),
					use_consensus   => $self->use_consensus(),
					pssmsets        => \@pssmsets,
				);
			} else {
				## use all other sequences as related sequences
				my @related  = ();
				my $sequence = $seq->seq();

				foreach my $seq2 (@sequences) {
					my $s2 = $seq2->seq();
					if ( $s2 ne $sequence ) {
						push @related, $s2;
					}
				}

				$b = Jobs::Subtasks::BiFa->new(
					sequence        => $sequence,
					threshold       => $self->threshold(),
					phylo_threshold => $self->phylo_threshold(),
					species         => $species_string,
					pssm_filter     => $self->pssm_filter(),
					use_consensus   => $self->use_consensus(),
					pssmsets        => \@pssmsets,
					phylogeneticly_related_sequences => \@related,
				);
			}

			push @{$result},
			  {
				sequence => $seq,
				bifa     => $b->run(),
			  };
		}

		return Serialization::Serializable_Array->new ( @{$result} );
	}

=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		my @arr = ();

		foreach my $res (@{$result->data}) {
			my $seq = $res->{sequence};
			my $bifahits = $res->{bifa}->{hits};

			my $startpos = $seq->{five_prime_pos};
			my $factor   = 1;
			if ($seq->{strand} eq 'negative') {
				$factor   = -1;
			}

			foreach my $hit (@{$bifahits}) {
				my $ug = Data::UUID->new();

				my $hit_start = 0;
				my $hit_end = 0;
				my $hit_strand = "positive";

				## hit strands will be relative to the sequence strand
				if ($hit->{'Positive_strand'}) {
					$hit_start = $startpos + $factor * $hit->{'Position'};
					$hit_end = $startpos + $factor * ( $hit->{'Position'} + $hit->{'Length'} );
					$hit_strand = $seq->{strand};
				} else {
					$hit_start = $startpos + $factor * ( $hit->{'Position'} + $hit->{'Length'} );
					$hit_end = $startpos + $factor * $hit->{'Position'};
					$hit_strand = $seq->{strand} eq 'positive' ? 'negative' : 'positive';
				}

				my $gann = Sequences::Annotation->new(
					$ug->create_str(),
					"bifahit",
					"BiFa: " . $hit->{'Name'} . " pval:" . $hit->{'Value'},
					$hit->{'Identifier'},
					$hit_start,
					$hit_end,
					$hit_strand,
				);

				#$gann->{hit} = Datatypes::Motifs::Motif_Hit->new(
				#	five_prime_pos => $hit_start,
				#	three_prime_pos => $hit_end,
				#	name => $hit->{'Name'},
				#	pvalue => $hit->{'Value'},
				#	strand => $hit_strand,
				#	sequence => "",
				#);

				$seq->add_annotation($gann);
			}

			push @arr, $seq;
		}

		$self->put_dataitem(
			Datatypes::Sequence::Set->new( sequences => \@arr ),
			$self->nc__output_dataitem() );

		return $result;
	}

=head2 Expiry time (disable cache)

(all caching is handled by the created BiFa subtasks)

=cut

	method expiry_time() {
		return 0;
	}

};
