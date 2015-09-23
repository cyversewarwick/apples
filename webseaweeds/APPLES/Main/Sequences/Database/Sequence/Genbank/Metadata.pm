#!/usr/bin/perl
package Sequences::Database::Sequence::Genbank::Metadata;

use strict;

use feature ':5.10';
use strict;

our @ISA = qw(Sequences::Database::Metadata);
require Sequences::Database::Metadata;

use Moose::Role;
use MooseX::Declare;


use Runtime;
use Carp;

use Serialization::Serializable;

=head2 Get all data available for an identifier

 Parameters:
 $self : a self object
 $id   : an id/accession
 
 Returns:
 {
 	url => an url (if exists),
 	common_name => common_name
 	type => 
 	...
 	
 }

=cut

sub _get_meta_data {
	my $self = shift;
	my $id   = shift;

	my $v = { 
		id => $id,
		url => "http://www.ncbi.nlm.nih.gov/gquery/?$id", 
	};

	my $query = "\"$id\" " . $self->{query};

	my $factory = Bio::DB::EUtilities->new(-eutil => 'esearch',
	                                       -email => 'mymail@foo.bar',
	                                       -db    => 'gene',
	                                       -term  => $query);
	my @ids       = $factory->get_ids;

	my $summaries = Bio::DB::EUtilities->new(
		-eutil => 'esummary',
		-db    => 'gene',
		-id    => \@ids
	);

	while ( my $docsum = $summaries->next_DocSum ) {
		## some items in DocSum are also named ChrStart so we pick the genomic
		## information item and get the coordinates from it
		my ($genomic_info) = $docsum->get_Items_by_name('GenomicInfoType');

		## some entries may have no data on genomic coordinates. This condition filters then out
		if ( !$genomic_info ) {
			## found no genomic coordinates data
			next;
		}

		## get coordinates of sequence
		## get_contents_by_name always returns a list

		my ($chr_acc_ver) = $genomic_info->get_contents_by_name("ChrAccVer");
		my ($chr_start)   = $genomic_info->get_contents_by_name("ChrStart");
		my ($chr_stop)    = $genomic_info->get_contents_by_name("ChrStop");
		
		$v->{acc_ver} = $chr_acc_ver;
		$v->{start} = $chr_start +1;
		$v->{end} = $chr_stop + 1;
	}
		 
	while (my $ds = $factory->next_DocSum) {
	    print "ID: ",$ds->get_id,"\n";
	    # flattened mode
	    while (my $item = $ds->next_Item('flattened'))  {
	        $v->{$item->get_name} = $item->get_content
	        	if $item->get_content;
	    }
	    print "\n";
	}

	my $object_type = "";

	given ($object_type) {
		when (/Gene/) {
#			$v->{nearest_genes} = [ $id ];
#			$v->{url} = "http://www.ensembl.org/Danio_rerio/Gene/Summary?g=$id";
#			$v->{geneid} = $id;
#
#			my $gene_adaptor =
#			  $registry->get_adaptor( $species, $db_type, 'Gene' );
#			my $gene = $gene_adaptor->fetch_by_stable_id($id);
#
#			$v->{biotype}              = $gene->biotype();
#			$v->{canonical_annotation} = $gene->canonical_annotation();
#			$v->{canonical_transcript} =
#			  $gene->canonical_transcript()->stable_id();
#			$v->{description} = $gene->description();
#			$v->{display_id}  = $gene->display_id();
#			$v->{source}      = $gene->source();
#			$v->{status}      = $gene->status();
#			$v->{version}     = $gene->version();
#
#			$v->{start}  = $gene->start();
#			$v->{end}    = $gene->end();
#			$v->{strand} = $gene->strand() < 0 ? "negative" : "positive";
#			$v->{length} = $gene->length();
#			$v->{transcripts} = [];
#			
#			my $transcripts = $gene->get_all_Transcripts();
#			while ( my $transcript = shift @{$transcripts} ) {
#				push @{$v->{transcripts}}, $transcript->stable_id(); 
#			}
		}
		when (/Transcript/) {
#			$v->{url} =
#			  "http://www.ensembl.org/Danio_rerio/Transcript/Summary?t=$id";
#
#			my $transcript_adaptor =
#			  $registry->get_adaptor( $species, $db_type, 'Transcript' );
#			my $gene = $transcript_adaptor->fetch_by_stable_id($id);
#
#			$v->{biotype}             = $gene->biotype();
#			$v->{coding_region_start} = $gene->coding_region_start();
#			$v->{coding_region_end}   = $gene->coding_region_end();
#			$v->{cdna_coding_start}   = $gene->cdna_coding_start();
#			$v->{cdna_coding_end}     = $gene->cdna_coding_end();
#
#			$v->{description} = $gene->description();
#			$v->{display_id}  = $gene->display_id();
#			$v->{seqname}     = $gene->seqname();
#			$v->{length}      = $gene->length();
#			$v->{version}     = $gene->version();
#			$v->{status}      = $gene->status();
#
#			$v->{start} = $gene->start();
#			$v->{end}   = $gene->end();
#
#			$v->{is_canonical} = $gene->is_canonical();
#			$v->{strand} = $gene->strand() < 0 ? "negative" : "positive";
#
#			my $fiv_utr = $gene->five_prime_utr();
#			my $thr_utr = $gene->three_prime_utr();
#			if (defined $fiv_utr) {
#				$v->{five_prime_utr} = $fiv_utr->seq();
#			}
#			if (defined $thr_utr) {
#				$v->{three_prime_utr} = $thr_utr->seq();
#			}
#			## TODO might want to retrieve this
#			#public Bio::Seq 	translate ()
#			#public Text 	translateable_seq ()
#			#public Bio::EnsEMBL::Translation 	translation ()
#
#			$v->{exons} = [];
#	        foreach my $exon ( @{ $gene->get_all_Exons() } ) {
#	        	push @{$v->{exons}}, $exon->stable_id(); 
#	        }
#
		}
		default {
#			$v->{url} = "http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q="
#			  . $id;
#
#			my $slice = $self->get_slice($id);
#			$v->{start}  = $slice->start();
#			$v->{end}    = $slice->end();
#			$v->{length} = $slice->length();
#			$v->{strand} = $slice->strand() < 0 ? "negative" : "positive";
		}
	}

#	if ( $object_type =~ m/Transcript/ || $object_type =~ m/Exon/ ) {
#		my $adaptor =
#		  $registry->get_adaptor( $species, $db_type, $object_type );
#		my $g  = $adaptor->fetch_by_stable_id($id);
#		my $ng = $g->get_nearest_Gene();
#		if ( defined($ng) ) {
#			if ( UNIVERSAL::isa( $ng, "Bio::EnsEMBL::Gene" ) ) {
#				push @{ $v->{nearest_genes} }, $ng->stable_id();
#			} elsif ( ref($ng) eq "ARRAY" ) {
#				push @{ $v->{nearest_genes} }, $_->stable_id() foreach @$ng;
#			} else {
#				warn ("Don't know what to do with nearest gene ref : $ng");
#			}
#		}
#	}

	bless $v, "Serialization::Serializable";

	return $v;
}

1;