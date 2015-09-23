#!/usr/bin/perl

=head1 Genomic sequence class

An extension of Bio::Seq to handle APPLES sequence database sources and
annotations

=cut

package Sequences::Genomic_Sequence;

use strict;

use Scalar::Util qw(blessed);
use List::Util qw(min max);
use Data::Dumper;
use POSIX;
use JSON;

use Bio::Seq;
use Bio::Species;

use Runtime;
use Carp;
use Serialization::Serializable;

our @ISA = qw(Serialization::Serializable Bio::Seq);

use Moose::Role;
with "Datatypes::Roles::Locatable";

=head2 Genomic Sequence Constructur

This constructor takes all args acceptable for Bio::Seq, and additionally:

 five_prime_pos : an absolute five prime position to associate with
 				  this sequence
 three_prime_pos : an absolute three prime position to associate with
 				   this sequence

 strand : the strand this sequence is on

 source : the source sequence database (a Database)

=cut

sub new {
	my ( $class, @args ) = @_;

	my $self = Bio::Seq->new(@args);
	return $class->new_from_seq( $self, @args );
}

=head2 Construct a Genomic_Sequence from a Bio::SeqI or such things.

 Parameters:
 $class : Genomic_Sequence
 $seq : a Bio::SeqI (or RichSeq, or ...)
 @args : BioPerl-friendly construction args

 Returns:
 a new Genomic_Sequence

=cut

sub new_from_seq {
	my ( $class, $self, @args ) = @_;

	bless $self, $class;

	$self->{annotation}      = undef;
	$self->{annotationindex} = {};

	my ( $fpp, $tpp, $species, $strand, $source ) =
	  &Bio::Root::RootI::_rearrange( $self,
		[qw(FIVE_PRIME_POS THREE_PRIME_POS SPECIES STRAND SOURCE SEQ)], @args );

	if ( defined $strand && $strand =~ m/^\-?\d+$/ ) {
		if ( $strand < 0 ) {
			$strand = "negative";
		} else {
			$strand = "positive";
		}
	}

	$self->{five_prime_pos}  = $fpp    || 0;
	$self->{three_prime_pos} = $tpp    || 0;
    $self->{species} = $species || "unknown";
	$self->{strand}          = $strand || "positive";
	$self->{source}          = $source;
	
	my $nstrand = ($self->{strand} eq "positive") ? 1 : -1;

	my $icounter = 1;
	for my $feat_object ($self->get_SeqFeatures) {          
		my $ptag = $feat_object->primary_tag;
				
		next unless ($ptag eq "CDS" || $ptag eq "mRNA" || $ptag eq "gene");

		my $sourceid = $self->accession_number();
		my $dispid = $self->display_id() || $sourceid;
		my $description = "";

		for my $tag ($feat_object->get_all_tags) { 
			if ($tag eq "gene") {
				for my $value ($feat_object->get_tag_values($tag)) {                
					if ($dispid ne "") {
						$dispid.= " ";
					}
					$dispid.= $value;
				}
			} else {
				$description.= "$tag :";
				for my $value ($feat_object->get_tag_values($tag)) {
					$description.= " " . $value;
				}
				$description.= ";";
		    }
		}
		
		my $type = "unknown";
		if ($ptag eq "mRNA" || $ptag eq "CDS") {
			$type = "exon"
		} elsif ($ptag eq "gene") {
			$type = "gene";
		}
		
		my @all_locations = $feat_object->location->each_Location;

		my $fstart = $feat_object->location->min_start;       
		my $fend = $feat_object->location->max_end;
		my $fstrand = $feat_object->location->strand;
		my $transcriptid = $sourceid;

		if ($type eq "exon") {
			$transcriptid = $sourceid . "_transcript_$icounter";
			$self->add_annotation(
				Sequences::Annotation->new(
					$dispid . "_transcript_$icounter",
					"transcript",
					$description,
					$transcriptid,
					$fstart + $fpp - 1,
					$fend + $fpp - 1,
					$fstrand < 0 ? "negative" : "positive",
					$feat_object
				)
			);
		}

		
		foreach my $location (@all_locations) {
			$fstart = $location->start;     
			$fend = $location->end;
			$fstrand = $location->strand;
			$fstrand *= $nstrand;

			$self->add_annotation(
				Sequences::Annotation->new(
					"${transcriptid}_${type}_$icounter",
					$type,
					"$dispid : $description",
					"${transcriptid}_${type}_$icounter",
					$fstart + $fpp - 1,
					$fend + $fpp - 1,
					$fstrand < 0 ? "negative" : "positive",
					$feat_object
				)
			);

			++$icounter;			
		}
		
	}

	return $self;
}

=head2 Add an annotation

 Parameters:
 $self : a Genomic_Sequence
 $entry : an Annotation

 Returns:

 Nothing

=cut

sub add_annotation {
	my ( $self, $entry, ) = @_;

	my $newnum = push @{ $self->{annotation} }, $entry;
	$self->{annotationindex}->{ $entry->{id} } = $newnum - 1;
}

## Helper: Rebuild annotation index
sub _build_ann_index {
	my $self = shift;
	my $i = 0;
	$self->{annotationindex} = {};
	foreach my $ann (@{ $self->{annotation} }) {
		$self->{annotationindex}->{ $ann->{id} } = $i++;
	}	
}

=head2 Delete one or more annotations

 Parameters:
 $self : a Genomic_Sequence
 $id   : an Annotation id or a coderef that takes an annotation as its parameter

 Returns:

 Nothing

=cut

sub delete_annotation {
	my ( $self, $id ) = @_;
	
	if (ref ($id) eq "CODE") {
		@{ $self->{annotation} } = grep { $id->($_) } @{ $self->{annotation} };			
	} else {
		@{ $self->{annotation} } = grep {$_->{id} eq $id} @{ $self->{annotation} };	
	}

	$self->_build_ann_index();
}

=head2 Get all annotation names

 Parameters:
 $self : a Genomic_Sequence

 Returns:
 A list of annotation ids.

=cut

sub get_annotation_names {
	my $self      = shift;
	my @ann_names = ();

	foreach my $ann ( @{ $self->{annotation} } ) {
		push @ann_names, $ann->{id};
	}
	return \@ann_names;
}

=head2 Get annotations by type

 Paramters:
 $self : a self object
 $type : a type string

 Returns:
 an ArrayRef to all annotations with the given type.

=cut

sub get_annotations_by_type {
	my $self      = shift;
	my $type      = shift;
	my @anns      = ();

	foreach my $ann ( @{ $self->{annotation} } ) {
		push @anns, $ann
			if (!defined ($type) || $ann->{type} eq $type);
	}
	return \@anns;

}

=head2 Reverse the strand (including locations of all annotations)

 Parameters:
 $self : a Genomic_Sequence

 Returns:
 A list of annotation ids.

=cut

sub reverse_strand {
	my $self = shift;

	my $new_five_prime_pos  = $self->{three_prime_pos};
	my $new_three_prime_pos = $self->{five_prime_pos};
	my $new_strand = $self->{strand} eq "positive" ? "negative" : "positive";

	$self->{five_prime_pos}  = $new_five_prime_pos;
	$self->{three_prime_pos} = $new_three_prime_pos;
	$self->{strand}          = $new_strand;

	$self->seq( reverse_dna_sequence_strand( $self->seq ) );

	if ( defined( $self->{sequence} ) ) {
		$self->{sequence} = $self->seq();
	}
	if ( defined( $self->{masked_sequence} ) ) {
		$self->{masked_sequence} =
		  reverse_dna_sequence_strand( $self->{masked_sequence} );
	}
}

=head2 Reverse the strand of a DNA Sequence (static method)

 Parameters:
 $DNA : the sequence to reverse the strand of

 Returns:
 Reverse complement of $DNA

=cut

sub reverse_dna_sequence_strand($) {
	my $DNA = shift;

	my $revcom = reverse $DNA;
	$revcom =~ tr/ACGTacgtNn/TGCAtgcaNn/;

	return $revcom;
}

=head2 This function returns a string which is compared to decide whether two genomic sequence objects can be merged

 Parameters:
 $self : a self object

 Returns: 
 A string identifying this sequence

=cut

sub ident {
	my $self = shift;
	return "$self->{sourceid}" . 
		   " : $self->{five_prime_pos} - ".
 		   "$self->{three_prime_pos} / " .
		   " $self->{strand} " .
		   ";" . ($self->{description} || "") .
			 to_json ( Serialization::Serializable::to_hash( $self->{source}, 1 ), {allow_blessed => 1} );	
}

=head2 Merge with another sequence object

 Parameters:
 $self : a sequence
 $other : another sequence

 Returns: 
 The two sequences, merged into $self

=cut

sub merge {
	my $self = shift;
	my @rhss = @_;
	
	my %all_anns = ();
	foreach my $ann ( @{ $self->{annotation} } ) {
		my $ann_ident = $ann->ident;
		$all_anns{$ann_ident} = $ann;			
	}
	my $merge = Hash::Merge->new('LEFT_PRECEDENT');

	my $rhs = shift @rhss;
	while (defined ($rhs)) {
		## we can only merge stuff when it comes from the same 
		## sequence base
		if ($self->ident ne $rhs->ident ){
			die "Cannot merge different sequences."
		}
		
		if ( $rhs->seq ne $self->seq ) {
			die "Data consistency error: tried to merge two matching sequence objects, but their sequence data isn't the same.";
		}
		
		$rhs = $rhs->cloneme();	## to make sure we have all in copies
		
		foreach my $ann ( @{ $rhs->{annotation} } ) {
			my $ann_ident = $ann->ident;
			if (!defined ($all_anns{$ann_ident})) {
				$all_anns{$ann_ident} = $ann;			
			} else {
				$all_anns{$ann_ident}->merge ($ann);
			}
		}
		
		$self->{metadata} = $merge->merge ($self->{metadata}, $rhs->{metadata});
			
		$rhs = shift @rhss;
	}
}

=head2 Serialization override

=cut

sub to_hash {
	my $self = shift;

	my $spec = $self->species;

	if ( UNIVERSAL::isa( $spec, "Bio::Species" ) ) {
		$spec = $spec->genus() . "_" . $spec->species();
	}
	
	if (!defined ($spec)) {
		$spec = "unknown species";
		if(UNIVERSAL::isa($self->{source}, 'Serialization::Serializable')) {
			$spec = $self->{source}->unique_id();
		}
		warn("Exporting Genomic_Sequence without defined species, inserting species name $spec");
	}

	return {
		five_prime_pos  => $self->{five_prime_pos},
		three_prime_pos => $self->{three_prime_pos},
		sequence_strand => $self->{strand},
		strand          => $self->{strand},
		source    => Serialization::Serializable::to_hash( $self->{source}, 1 ),
		sequence  => $self->seq,
		sourceid  => $self->id(),
		accession => $self->accession_number,
		description     => $self->desc,
		species         => $spec,
		masked_sequence => $self->{masked_sequence} || $self->seq,
		annotation =>
		  Serialization::Serializable::to_hash( $self->{annotation}, 1 ),
		annotationindex => $self->{annotationindex},
		metadata => Serialization::Serializable::to_hash( $self->{metadata} ),
	};
}

=head2 Serialization override

=cut

sub serialization_restore {
	my $oldself = shift;

	my $self = Sequences::Genomic_Sequence->new(
		-id              => $oldself->{sourceid},
		-species         => $oldself->{species},
		-seq             => $oldself->{sequence},
		-five_prime_pos  => $oldself->{five_prime_pos},
		-three_prime_pos => $oldself->{three_prime_pos},
		-strand          => $oldself->{strand},
		-source          => $oldself->{source},
	);

	while ( my ( $k, $v ) = each(%$oldself) ) {
		$self->{$k} = $v;
	}

	$self->accession_number( $oldself->{accession} );

	%$oldself = %$self;

}

1;
