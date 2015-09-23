#!/usr/bin/perl

=head1 Class Sequence_Database_Genbank

[document class here]
=cut

package Sequences::Database::Sequence::Genbank;

use strict;

use feature ":5.10";

use File::Temp qw (tmpnam);
use Data::Dumper;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::GenBank;
use Bio::DB::EUtilities;

use Runtime;
use Carp;

use Sequences::Database::Sequence;
use Sequences::Genomic_Sequence;
use Sequences::Annotation;

our @ISA = qw(Sequences::Database::Sequence);

use Moose::Role;

with "Sequences::Database::Sequence::Genbank::Metadata";
with "Sequences::Database::Metadata";

=head2 Constructor

Parameters:
 $class : Sequence_Database_Genbank
 
=cut

sub new {
	my ( $class, $query ) = @_;
	$query = $query || "";
	my $self = { query => $query };
	bless $self, $class;
	return $self;
}

=head2 Get a sequence by its accession from the database

Parameters:
$self : A Database
$location : a Relative_Location

Returns : 
Nothing, but dies if location is invalid

=cut

sub validate_location {
	my $self = shift;
	my $loc  = shift;

	## TODO this should be done faster/better/nicer 
	$self->_get_sequence_by_location ($loc);

	return;
}


=head2 Get a single sequence by sequence location

Location objects specify locations of sequences relative to 
one or more accessions.  

 Parameters:
 $self : a Database
 $location : a Location
 
 Returns:
 An ARRAYREF [Genomic_Sequence]

=cut

sub _get_sequence_by_location {
	my $self     = shift;
	my $location = shift;

	my @results = ();
	our $AR;

	my $acc    = $location->identifier;
	my $anchor = $location->anchor;
	my $offset = $location->offset;
	my $length = $location->length;

	## results limitations (stop if query returns more than this number of results)
	## we should really only get one.
	my $limit = 1;

	## query for the database
	my $query = "\"$acc\" " . $self->{query};
	
	my @dbs = qw(gene nucleotide);

	foreach my $db (@dbs) {
		
		## make the search
		my $factory = Bio::DB::EUtilities->new(
			-eutil  => 'esearch',
			-db     => $db,
			-term   => $query,
			-tool   => 'bioperl',
			-retmax => $limit
		);
	
		my $n_results = $factory->get_count;
		my @ids       = $factory->get_ids;
			
		my $summaries = Bio::DB::EUtilities->new(
			-eutil => 'esummary',
			-db    => $db,
			-id    => \@ids
		);

		while ( my $docsum = $summaries->next_DocSum ) {
			print Dumper ($docsum);
			
			## some items in DocSum are also named ChrStart so we pick the genomic
			## information item and get the coordinates from it
			my ($genomic_info) = $docsum->get_Items_by_name('GenomicInfoType');

			## some entries may have no data on genomic coordinates. This condition filters then out
			
			if ( defined ($genomic_info) ) {
				## get coordinates of sequence
				## get_contents_by_name always returns a list

				my ($chr_acc_ver) = $genomic_info->get_contents_by_name("ChrAccVer");
				my ($chr_start)   = $genomic_info->get_contents_by_name("ChrStart");
				my ($chr_stop)    = $genomic_info->get_contents_by_name("ChrStop");

				my $up_end     = $chr_start;
				my $down_start = $chr_stop;

				my $orig_len = abs ($chr_start - $chr_stop) + 1;
				my $u        = 0;
				my $d        = 0;

				unless  ( $anchor eq "5' end"  || $anchor eq "3' end" ) {
					die "Unsupported sequence anchor : " . $anchor;
				}

				my $strand;
				if ( $chr_start < $chr_stop ) {
					$strand    = 1;

					if ( $anchor eq "5' end" ) {
						$chr_start = $chr_start + $offset;
					} else {
						$chr_start = $chr_stop + $offset;
					}
					$chr_stop  = 1 + $chr_start + $length;
				} else {
					$strand    = 2;			

					if ( $anchor eq "5' end" ) {
						$chr_start = $chr_start - $offset;
					} else {
						$chr_start = $chr_stop - $offset;
					}
					$chr_stop  = 1 + $chr_start - $length;
				}
				my $nstrand = ( $strand == 1 ? 1 : -1 );

				next if ( abs ($chr_start - $chr_stop) + 1  < 0); 

				my $obj = $self->_fetch_nucleotide ($chr_acc_ver, $chr_start, $chr_stop, $strand);

				return $self->_make_sequence( $obj, $chr_start, $chr_stop, $nstrand, $up_end, $down_start, $location );		
			}
			
			if ($db eq 'nucleotide') {
				my ($caption) = $docsum->get_Items_by_name("Caption");
				my ($len) = $docsum->get_Items_by_name("Length");
				
				$len = $len->get_content if defined $len;
				$caption = $caption->get_content if defined $caption;

				if( $acc =~ m/^$caption/
				&&  defined $len && $len > 0) {
					my $obj = $self->_fetch_nucleotide ($acc);
					return $self->_make_sequence( $obj, 1, $len, 1, 1, $len, $location );		
				}
			}
		}
	}
	die "Cannot find " . $location->identifier . " in genbank.";
}

sub _fetch_nucleotide {
	my $self = shift;
	my $chr_acc_ver = shift;
	my $chr_start = shift;
	my $chr_stop = shift;
	my $strand = shift;
	
	my $fetcher = Bio::DB::EUtilities->new(-eutil   => 'efetch',
	                                       -email => 'mymail@foo.bar',
	                                       -db      => 'nucleotide',
	                                       -rettype => 'gbwithparts');
	
	my @parameters = ( -id => $chr_acc_ver );
	
	if (defined $chr_start) {
		push @parameters, -seq_start => $chr_start + 1;
	}
	if (defined $chr_stop) {
		push @parameters, -seq_stop => $chr_stop + 1;
	}
	if (defined $strand) {
		push @parameters, -strand => $strand + 1;
	}
	
	$fetcher->set_parameters( @parameters );

	use File::Temp qw(tempfile);
	my ($fh, $filename) = tempfile();
	close $fh;
	my $gbdata = $fetcher->get_Response(file => $filename);
	my $seqin = Bio::SeqIO->new(-file   => $filename,
	                            -format => 'genbank');
	my $obj = undef;
	while (my $seq = $seqin->next_seq) {
		if (!defined ($obj)) {
			$obj = $seq;
		} else {
			die "More than one sequence was returned.";
		}
	}
	
	if ( (defined $chr_start && defined ($chr_stop)) &&  
		length ($obj->seq) != abs ($chr_start - $chr_stop) + 1  ) {
		die ("Returned sequence lengths don't match : " . (length ($obj->seq)) . " != " .( abs ($chr_start - $chr_stop) + 1));
	}
	
	return $obj;
}

sub _make_sequence {
	my $self = shift;
	my $obj = shift;
	my $chr_start = shift;
	my $chr_stop = shift;
	my $nstrand = shift;
	my $up_end = shift;
	my $down_start = shift;
	my $location = shift;

	my $chr_acc_ver   = $location->identifier;
	
	
	my $gs = Sequences::Genomic_Sequence->new_from_seq(
		$obj,
		(
			-strand          => $nstrand,
			-source          => { db => $self, location => $location },
			-five_prime_pos  => $chr_start,
			-three_prime_pos => $chr_stop,
		)
	);
	
	my $gann = Sequences::Annotation->new(
		$chr_acc_ver,
		"gene",
		$obj->desc(),
		$chr_acc_ver,
		$up_end + 1,
		$down_start + 1,
		$nstrand < 0 ? "negative" : "positive",
	);
	
	$gs->add_annotation($gann);
	
	return $gs;
}

=head2 Getter/Setter for query member

 Parameters:
 $self : $self object
 $val : new value (optional)

=cut

sub query {
	my $self = shift;
	my $val  = shift;

	if ( defined($val) ) {
		$self->{query} = $val;
	}
	return $self->{query};
}

=head2 Overridden from Sequence_Database
=cut

sub name {
	my $self = shift;

	return "Genbank nucleotide"
	  . ( ( $self->{'query'} ? "('$self->{query}')" : "" ) );
}

=head2 Serialization override

=cut

sub to_hash {
	my $self = shift;
	return { 'query' => $self->{query}, };
}

=head2 Serialization override

=cut

sub serialization_restore {
	my $oldself = shift;

	my $self = Sequences::Database::Sequence::Genbank->new( $oldself->{query} );
	%$oldself = %$self;
}

1;
