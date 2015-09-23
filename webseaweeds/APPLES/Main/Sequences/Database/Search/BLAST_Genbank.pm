#/usr/bin/perl

=head1 Blast search database abstraction

We use a local installation of BLAST, and retrieve Sequences from Genbank.

=cut

package Sequences::Database::Search::BLAST_Genbank;

use strict;
no warnings;
use feature ":5.10";

our @ISA =
  qw(Sequences::Database::Sequence::Genbank Sequences::Database::Search::BLAST_Fasta Serialization::Serializable);
require Sequences::Database::Sequence::Genbank;
require Sequences::Database::Search::BLAST_Fasta;
require Serialization::Serializable;

use File::Temp qw (tmpnam);
use Data::Dumper;
use Cwd qw(abs_path);

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::DB::GenBank;

use Runtime;
use Carp;

=head2 Constructor

 Parameters:
 $class : Sequence_Database_BLAST_Genbank
 @args : BioPerl friendly args hash, containing:
 
  db : the database to search (default : refseq_rna)
  entrezquery : an entrez query to restrict results

 and optionally (defaults from APPLES.dat):
  fastapath : path to FASTA library 
  program : the BLAST program to use
  dbpath : the Blast DB path (set to '-remote' for remote Blasting)

 also, optionally:
  entrezquery : an entrez query to restrict the search results

=cut

sub new {

	my ( $class, @args ) = @_;

	# Bio::Root::RootI::_rearrange does not use this.
	my $dummy = undef;

	my ( $db, $entrezquery ) =
	  &Bio::Root::RootI::_rearrange( $dummy, [qw(DB ENTREZQUERY)], @args );

	my $db          = $db          || 'refseq_rna';
	my $entrezquery = $entrezquery || undef;

	my $self = Sequences::Database::Search::BLAST_Fasta->new(@args);

	bless $self, "Sequences::Database::Sequence::Genbank";
	bless $self, $class;
	$self->{db}          = $db;
	$self->{entrezquery} = $entrezquery;

	return $self;
}

=head2 Search the database for a Bio::Seq

See L<Alignment::BlastDB> 

=cut

sub _search {
	my $self   = shift;
	my $seq    = shift;
	my $params = shift || {};
	my $parameters;
	%{$parameters} = %${params};

	my $tmpfile = tmpnam();
	my $out     = Bio::SeqIO->new(
		-file   => ">$tmpfile",
		-format => 'Fasta'
	);

	$out->write_seq($seq);

	my $tmpfile2 = tmpnam();

	$parameters->{-query}    = $tmpfile;
	$parameters->{-db}       = $self->{db};
	$parameters->{-out}      = $tmpfile2;
	$parameters->{-outfmt}   = "5";
	$parameters->{"-remote"} = "";
	if ( exists( $self->{entrezquery} )
		&& defined( $self->{entrezquery} ) )
	{
		$parameters->{'-entrez_query'} = $self->{entrezquery};
	}

	$self->_run( $self->{program}, $parameters );

	my $result =
	  Sequences::Database::Search::BLAST_Fasta::_read_blast_result($tmpfile2)
	  ;

	unlink($tmpfile);
	unlink($tmpfile2);
	return $result;
}

=head2 Serialization override

=cut

sub to_hash {
	my $self = shift;

	my $db          = $self->{db};
	my $program     = $self->{program};
	my $entrezquery = $self->{entrezquery};
	my $dbpath      = $self->{dbpath};

	if ( $dbpath eq abs_path( config->{blast_db} ) ) {
		$dbpath = undef;
	}

	my $fastapath = $self->{fastapath};

	return {
		db          => $db,
		program     => $program,
		entrezquery => $entrezquery,
		dbpath      => $dbpath,
	  }

}

=head2 Serialization override

=cut

sub serialization_restore {
	my $self = shift;
	%$self = %{ Sequences::Database::Search::BLAST_Genbank->new(
			-db          => $self->{db},
			-program     => $self->{program},
			-dbpath      => $self->{dbpath},
			-entrezquery => $self->{entrezquery},
		)
	  };
}

1;
