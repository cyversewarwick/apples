#!/usr/bin/perl

use MooseX::Declare;

=head1 Sequence Location Class

This is a proxy to hand self-created sequences (which might not have come
from a database) as sequence locations.

=cut

class Datatypes::Sequence::Local extends Datatypes::Sequence::Location {
    use Data::Dumper;

	with "Datatypes::Roles::Locatable";

	use Runtime;

	use Scalar::Util qw(blessed);
	use Datatypes::Moose_Types;
	use Data::UUID;
	use File::Temp qw(tempfile);

	use Bio::SeqIO;

	require Sequences::Database::Relative_Location;
	
	## local sequence data
	has 'local_sequence' => (
		is => "rw",
		isa => "Sequences::Genomic_Sequence",
		required => 1,
	);

=head2 Create new local "location" from sequence

 Parameters:
 $class  : Datatypes::Sequence::Local
 $data   : a string containing sequence data
 $species : the species of the sequence (or undef)
 $strand : The strand (or undef)
 $five_prime_pos : The five prime end position (or undef)

 Returns:
 a new Datatypes::Sequence::Local

=cut	
	sub new_from_string {
		my $class = shift;
		my $data = shift
			or die "Cannot create Datatypes::Sequence::Local without sequence data.";
		my $species = shift;
		my $strand = shift || "positive";
		my $five_prime_pos = shift || 1;
		my $nstrand = ($strand eq "positive") ? 1 : -1;
		
		my ($fh, $filename) = tempfile();
		print $fh $data;
		close $fh;
		
		my $seqio_object = Bio::SeqIO->new(-file => $filename);
		my $seq          = $seqio_object->next_seq;
		
		if (!defined ($seq)) {
			die "Could not resolve $data to Bio::Seq object.";
		}
		
		my $id = $seq->id || "locally added sequence in $species";
		
		if (!defined ($species)) {
			my $ug = Data::UUID->new();
			$species = "unknown species " . $ug->create_str();
		}
		
		my $seq_obj = Sequences::Genomic_Sequence->new_from_seq (
			$seq,
			-id => $id,
			-species => $species,
			-five_prime_pos => $five_prime_pos,
			-three_prime_pos => $five_prime_pos + $nstrand * (length($data) - 1),
			-strand => $strand,
			);
        
 
		return $class->new_from_sequence_object( $seq_obj );
	}

=head2 Create new local "location" from sequence

 Parameters:
 $class : Datatypes::Sequence::Local
 $sequence : a Sequences::Genomic_Sequence object

 Returns:
 a new Datatypes::Sequence::Local

=cut	
	sub new_from_sequence_object {
		my $class = shift;
		my $sequence = shift;
		
		my $id = $sequence->id;
		
		if (!defined ($id)) {
			my $ug = Data::UUID->new();
			$id = $ug->create_str();			
		}
		
		my $loc = Sequences::Database::Relative_Location->new (
			id => $id,
			length => length ($sequence->seq),
		);
		
		return __PACKAGE__->new (
			db => "none",
			location => $loc,
			local_sequence => $sequence,
		);
	}

=head2 Get the sequence for this location

This overrides Datatypes::Sequence::Location to do nothing.

=cut

	method validate () {
	}

=head2 Get the sequence for this location

 Parameters:
 $self : $self object

 Returns:
 [ a Genomic_Sequence corresponding to this location ]

=cut

	method get_sequence () {
		return [ $self->local_sequence ];
	}

=head2 Get the metadata for this location

 Parameters:
 $self : $self object

 Returns:
 a metadata hash (see Sequences::Database::Metadata)
 corresponding to this location

 There might not be any metadata available, in which
 case we will return {}

=cut

	method get_meta_data () {
		return $self->local_sequence->{metadata} || {};
	}


};
