#/usr/bin/perl

=head1 Blast search database abstraction

We use a local installation of BLAST.

=cut

package Sequences::Database::Search::BLAST_Fasta;

use strict;
no warnings;
use feature ":5.10";

use Runtime;
use Carp;

our @ISA =
  qw(Sequences::Database::Search Sequences::Database::Sequence::FASTA Serialization::Serializable);
require Sequences::Database::Search;
require Sequences::Database::Sequence::FASTA;
require Serialization::Serializable;

use File::Which;
use File::Spec;
use File::Temp qw( :POSIX );
use Cwd qw(getcwd abs_path);

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::DB::GenBank;

=head2 Constructor

 Parameters:
 $class : Alignment::DB
 @args : BioPerl friendly args hash, containing:

 fastaname : name of the underlying FASTA file
 
 and optionally (defaults from APPLES.dat):
  fastapath : path to FASTA library 
  program : the BLAST program to use
  dbpath : the Blast DB path

=cut

sub new {
	my ( $class, @args ) = @_;

	# Bio::Root::RootI::_rearrange does not use this.
	my $dummy = undef;

	my ( $fastaname, $program, $dbpath, $fastapath ) =
	  &Bio::Root::RootI::_rearrange( $dummy,
		[qw(FASTANAME PROGRAM DBPATH FASTAPATH)], @args );

	$program   = $program   || 'blastn';
	$dbpath    = $dbpath    || config->{blast_db};
	$fastapath = $fastapath || config->{fastapath};

	my $self;
	if(!defined ($fastaname)) {
		$dbpath = "-remote";
	}

	if ( $dbpath ne "-remote" ) {
		$self =
		  Sequences::Database::Sequence::FASTA->new(
			File::Spec->catfile( $fastapath, $fastaname ) );
		$self->{dbpath}         = abs_path($dbpath);
		$self->{fastapath}      = abs_path($fastapath);
		$self->{shortfastaname} = $fastaname;
	
		my $db_file_path =
		  File::Spec->catfile( $self->{dbpath}, $self->{shortfastaname} );
	
		$self->{db_file_path} = $db_file_path;
	} else {
		$self = { dbpath => '-remote', };
	}

	$self->{blastdir}       = _find_blast();
	$self->{program}        = $program;

	bless $self, "Sequences::Database::Search";
	bless $self, $class;

	$self->_check_create_database();

	return $self;
}

=head2 Search the database for a Bio::Seq

Parameters:
$self : a LocalBlastDB
$seq : a Bio::Seq or Bio::PrimarySeq
$parameters : additional parameters to pass to BLAST

Returns:
a Bio::SearchIO result

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

	$parameters->{-query}         = $tmpfile;
	$parameters->{-db}            = $self->{db_file_path};
	$parameters->{-out}           = $tmpfile2;
	$parameters->{-outfmt}        = "5";
	$parameters->{-num_threads}       = "4";

	$self->_run( $self->{program}, $parameters );

	my $result = _read_blast_result ($tmpfile2);

	unlink($tmpfile);
	unlink($tmpfile2);
	return $result;
}

=head2 Read results using Bio::SearchIO

 Parameters:
 $filename : the name of the result file to read
 
 Returns:
 the results in an ARRAYREF
 
=cut

sub _read_blast_result {
	my $filename = shift;

	my @rdata = ();

	my $in = Bio::SearchIO->new(
		-format => 'blastxml',
		-file   => $filename
	);

	while ( my $res = $in->next_result ) {
		## $result is a Bio::Search::Result::ResultI compliant object
		my @hits = ();
		while ( my $hit = $res->next_hit ) {
			eval {
			my @in = $hit->seq_inds( 'query', 'cons', 1 );
			my $inds = \@in || [];
			my @int = $hit->seq_inds( 'sbjct', 'cons', 1 );
			my $indst = \@int || [];
			my $h_hit = {
				hit_name     => $hit->name() || "unknown",
				description  => $hit->description() || "unknown",
				length       => $hit->length() || "unknown",
				algorithm    => $hit->algorithm() || "unknown",
				score        => $hit->raw_score() || "unknown",
				significance => $hit->significance() || "unknown",
				rank         => $hit->rank() || "unknown",
				accession    => $hit->accession() || "unknown",
				seq_inds     => $inds,
				seq_inds_t   => $indst,
				strand       => $hit->strand() || "unknown",
				hsps         => [],
			};
			## $hit is a Bio::Search::Hit::HitI compliant object
			while ( my $hsp = $hit->next_hsp ) {
				## $hsp is a Bio::Search::HSP::HSPI compliant object
				push @{ $h_hit->{hsps} },
				  {
					r_type          => $hsp->algorithm() || "unknown",
					pvalue          => $hsp->pvalue() || "0",
					evalue          => $hsp->evalue() || "unknown",
					frac_id_query   => $hsp->frac_identical('query') || "unknown",
					frac_id_hit     => $hsp->frac_identical('hit') || "unknown",
					frac_id_total   => $hsp->frac_identical('total') || "unknown",
					frac_cons_query => $hsp->frac_conserved('query') || "unknown",
					frac_cons_hit   => $hsp->frac_conserved('hit') || "unknown",
					frac_cons_total => $hsp->frac_conserved('total') || "unknown",
					gaps_query      => $hsp->gaps('query') || "unknown",
					gaps_hit        => $hsp->gaps('hit') || "unknown",
					gaps_total      => $hsp->gaps('total') || "unknown",
					qseq            => $hsp->query_string() || "unknown",
					hseq            => $hsp->hit_string() || "unknown",
					homology_string => $hsp->homology_string() || "unknown",
					length_query    => $hsp->length('query') || "unknown",
					length_hit      => $hsp->length('hit') || "unknown",
					length_total    => $hsp->length('total') || "unknown",
					rank            => $hsp->rank,
				  };
			}

			push @hits, $h_hit;
			};
			if($@) {
				warn("ignoring a hit: $@");
			}
		}

		my $res_data = {
			id          => $res->query_name(),
			description => $res->query_description(),
			dbname      => $res->database_name(),
			size        => $res->database_letters(),
			num_entries => $res->database_entries(),
			parameters  => {},
			statistics  => {},
			hits        => \@hits,
		};

		#			foreach ( @{ $res->available_parameters } ) {
		#				$res_data->{parameters}->{$_} = $res->get_parameter($_);
		#			}
		#
		#			foreach ( @{ $res->available_statistics } ) {
		#				$res_data->{statistics}->{$_} = $res->get_statistic($_);
		#			}
		push @rdata, $res_data;
	}
	
	return \@rdata;
}

=head2 Check or create if we need to create the database 

=cut

sub _check_create_database {
	my $self = shift;

	my $dbtype = "nucl";

	given ( $self->{program} ) {
		when (/tblastx/) {
			$dbtype = "nucl";
		}
		when (/tblastn/) {
			$dbtype = "nucl";
		}
		when (/blastn/) {
			$dbtype = "nucl";
		}
		when (/blastp/) {
			$dbtype = "prot";
		}
		when (/blastx/) {
			$dbtype = "prot";
		}
	}

	unless ($self->{dbpath} eq "-remote" || -f "$self->{db_file_path}_db" ) {
		$self->_run(
			'makeblastdb',
			{
				'-in'           => $self->{fastaname},
				'-out'          => $self->{db_file_path},
				'-parse_seqids' => "",
				'-title'        => $self->{fastaname},
				'-hash_index'   => "",
				'-dbtype'       => $dbtype,
			}
		);
		qx (touch $self->{db_file_path}_db);
	}
}

=head2 Find the blast installation

Returns:

The path to the folder that contains blastn

=cut

sub _find_blast {
	my $loc = which('blastn') || 'blastn';
	if ( -x $loc ) {
		$loc =~ s/blastn$//;
		return $loc;
	}
	my $bd = config->{blast} || $ENV{BLASTDIR};
	if ($bd) {
		if ( -x File::Spec->catfile( $bd, "blastn" ) ) {
			return $bd;
		}
		if ( -x File::Spec->catfile( $bd, "bin", "blastn" ) ) {
			return File::Spec->catdir( $bd, "bin" );
		}
	}

	confess(
		"Blast installation could not be found. Check APPLES.dat : blast\n");
}

=head2 Run a BLAST executable

Parameters:
$self : a LocalBlastDB hash
$program : the program
%$parameters : the hash of parameters

Returns: 
the program output

=cut

sub _run {
	my $self       = shift;
	my $program    = shift;
	my $parameters = shift;

	my $paramstr = "";

	my $cwd = getcwd();

	chdir( $self->{dbpath} );

	my @args = ();

	while ( my ( $k, $v ) = each(%$parameters) ) {
		unless ( $k =~ m/^\-/ ) {
			$k = "-" . $k;
		}
		push @args, $k;
		push @args, $v
		  unless $v eq "";
	}

	my $prog = File::Spec->catfile( "$self->{blastdir}", $program );

	debug("Running: $prog  @args");
	my $result = system $prog, @args;

	chdir($cwd);

	return $result;
}

=head2 Get a name/unique identfier for this database

 Returns:
 A name (needs to be overloaded by implementing classes)

=cut

sub name {
	my $self = shift;
	return $self->{shortfastaname};
}

=head2 Serialization override

=cut

sub to_hash {
	my $self = shift;

	my $fastaname = $self->{shortfastaname};
	my $program   = $self->{program};

	my $dbpath = $self->{dbpath};

	if ( $dbpath eq abs_path( config->{blast_db} ) ) {
		$dbpath = undef;
	}

	my $fastapath = $self->{fastapath};

	if ( $fastapath eq abs_path( config->{fastapath} ) ) {
		$fastapath = undef;
	}

	return {
		fastaname => $fastaname,
		program   => $program,
		dbpath    => $dbpath,
		fastapath => $fastapath,
	};

}

=head2 Serialization override

=cut

sub serialization_restore {
	my $self = shift;

	%$self = %{ Sequences::Database::Search::BLAST_Fasta->new(
			-fastaname => $self->{fastaname},
			-program   => $self->{program},
			-dbpath    => $self->{dbpath},
			-fastaname => $self->{fastaname},
		)
	  };
}

1;
