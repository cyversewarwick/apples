#!/usr/bin/perl

use MooseX::Declare;

=head1 Sequence Location Class

Stores sequence locations and databases.

=cut

class Datatypes::Concepts::Gene extends Datatypes::Sequence::Location {

	use Runtime;
	use Hash::Merge;
	use Serialization::Serializable;
	use Data::Dumper;

	require Sequences::Database::Relative_Location;
	require Datatypes::Sequence::Set;

	## force metadata retrieval
	has 'metadata' => (
		is => "rw",
		isa => "Any",
		builder => "get_meta_data",
		lazy => 1,
	);
	
	## we can associate Genomic_Sequence objects with a gene
	has 'related_sequences' => (
		is => 'rw',
		isa => 'Datatypes::Sequence::Set',
		default => sub { 
			return Datatypes::Sequence::Set->new ();
		},
	);

=head2 Merge two genes right to left

All data in $rhs is taken over into $self.

 Parameters:
 $rhs : a Datatypes::Concepts::Gene

 Returns: nothing

=cut	

	method merge (Datatypes::Concepts::Gene $rhs) {
		if ($rhs->accession ne $self->accession) {
			die "Cannot merge data from gene " . $rhs->unique_id . " with different gene " . $self->unique_id;
		}
		
		if ( Serialization::Serializable::_weakened_memory_cycle_exists($self->metadata ) ) {
			die "Cannot merge cyclic metadata object : " . Dumper ($self->metadata);
		}

		if ( Serialization::Serializable::_weakened_memory_cycle_exists($rhs->metadata ) ) {
			die "Cannot merge cyclic metadata object : " . Dumper ($rhs->metadata);
		}
		
		my $merge = Hash::Merge->new('LEFT_PRECEDENT');
		$self->metadata ($merge->merge ($self->metadata, $rhs->metadata));
		
		my @sequences = @{$self->related_sequences->sequences()};
		push @sequences, $_ foreach @{$rhs->related_sequences->sequences()};
		$self->related_sequences->sequences(\@sequences);
		$self->related_sequences->compact_sequences;
	}

=head2 Make a new Gene (shortcut)

 Parameters:
 $db : the database identifier
 $acc : the gene accession

 Returns:
 a new Gene object.

=cut

	sub new_gene {
		my $class = shift;
		my $db = shift;
		my $acc = shift;
		my $v = Datatypes::Concepts::Gene->new(
			db => $db,
			location =>
			  Sequences::Database::Relative_Location->new( identifier => $acc )
		);
		return $v;
	}

	## override this to change the color in derived classes
	method color () {
		return "#ccf";
	}
	
	## overridden to_hash to retrieve data first
	sub to_hash {
		my $self = shift;
		$self->metadata;
		$self->related_sequences;
		return $self->SUPER::to_hash(@_);
	}

=head2 Create a dot label

=cut

	method dot_style () {
		return
		    "label=\"Gene: "
		  . $self->accession . " ("
		  . $self->db . ")"
		  . "\" shape=\"ellipse\"";
	}

}
