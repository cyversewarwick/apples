#!/usr/bin/perl 

=head1 Interface definition for a Blast-type sequence searchable database

=cut

package Sequences::Database::Search;

use strict;
no warnings;
use feature ':5.10';

use POSIX;
use Digest::SHA;
use JSON;

use Bio::Seq;

use Runtime;
use Carp;

use Serialization::Serializable;
use Serialization::Serializable_Array;

=head2 Search the database for a Bio::Seq, including cache lookup

 Parameters:
 $self : a BlastDB
 $seq : a Genomic_Sequence
 $parameters : (optional) additional BLAST parameters

 Returns:
 an ARRAYREF of results, hits, hsps...:
 
 [
 
 {
	id          => $res->query_name(),
	description => $res->query_description(),
	dbname      => $res->database_name(),
	size        => $res->database_letters(),
	num_entries => $res->database_entries(),
	parameters  => { ... },
	statistics  => { ... },
	hits        =>  [
 	{
		hit_name => $hit->name(),
		description => $hit->description(),
		length => $hit->length(),
    	algorithm => $hit->algorithm(),
		score => $hit->raw_score(),
		significance => $hit->significance(),
		accession    => $hit->accession(),
    	rank => $hit->rank(),
    	hsps => [
    		{
			r_type => $hsp->algorithm(),
			pvalue => $hsp->pvalue(),
			evalue = $hsp->evalue(),
    		frac_id_query => $hsp->frac_identical( 'query' )
    		frac_id_hit => $hsp->frac_identical( 'hit' )
    		frac_id_total => $hsp->frac_identical( 'total' )
	    	frac_cons_query => $hsp->frac_conserved( 'query' );
	    	frac_cons_hit => $hsp->frac_conserved( 'hit' );
	    	frac_cons_total => $hsp->frac_conserved( 'total' );
	    	gaps_query => $hsp->gaps( 'query' );
	    	gaps_hit => $hsp->gaps( 'hit' );
	    	gaps_total => $hsp->gaps( 'total' );
			qseq => $hsp->query_string();
			hseq => $hsp->hit_string();
			homology_string => $hsp->homology_string();
    		length_query => $hsp->length( 'query' );
    		length_hit => $hsp->length( 'hit' );
    		length_total => $hsp->length( 'total' );
			rank => $hsp->rank;   
    		}
    	]
 	}, ...
 ]
 
 ],

=cut

sub search {
	my $self       = shift;
	my $seq        = shift;
	my $parameters = shift || {};

	my $cache = cache;

	unless ( UNIVERSAL::isa( $seq, "Serialization::Serializable" ) ) {
		confess("Type check failed.");
	}

	my $seq_json = to_json(
		Serialization::Serializable::to_hash($seq),
		{ canonical => 1, allow_blessed => 1 }
	);
	my $digest = Digest::SHA->new (512);
	$digest->add($seq_json);

	my $cache_parameters = bless {
		"searchdb"         => $self,
		"blast_parameters" => $parameters,
		"sequence_hash"    => $digest->hexdigest(),
	  },
	  "Serialization::Serializable";

	# optionally, we might save the search sequence by its hash here.
	# This is a bit of a waste of space though with no obvious applications.

	my $result = $cache->cache_get($cache_parameters);

	if ( !defined($result) ) {
		info( "Computing BLAST hits for " . $seq->id . " in " . $self->name 
			. " " . to_json($cache_parameters, {allow_blessed=>1}) );
		my $rdata = $self->_search( $seq, $parameters );

		my $hitnum = 0;

		foreach my $res (@$rdata) {
			$hitnum += scalar @{ $res->{hits} };
		}

		$result = Serialization::Serializable_Array->new;
		$result->data($rdata);

		if ( int($hitnum) == 0 ) {
			warn("No blast hits found -- that is most likely an error.");
		} else {
			## if no blast hits were found, this can have multiple reasons
			## one of them is that BLAST got aborted

			## we don't want empty results in the cache.
			$cache->cache_put( $cache_parameters, $result );
		}

	} else {
		info(
			"Using cached BLAST hits for " . $seq->id . " in " . $self->name . " " 
				. to_json($cache_parameters, {allow_blessed=>1}) );
	}

	return $result->data;
}

=head2 Search the database for a Bio::Seq

This should be implemented by derived classes.

 Parameters:
 $self : a BlastDB
 $seq : a Bio::Seq or Bio::PrimarySeq
 $parameters : additional parameters to pass to BLAST

Returns:
a Bio::SearchIO result

=cut

sub _search {
	confess("Abstract method called.");
}

=head2 Get a name/unique identfier for this database

 Returns:
 A name (needs to be overloaded by implementing classes)

=cut

sub name {
	confess("Abstract method called");
}

=head2 Get the n best hits for a Bio::Seq in a Search_Database

 Parameters:
 $sdb : a Search_Database
 $seq : a Genomic_Sequence
 $top_x : number of hits to return
 $exclude : coderef to exclude e.g. same species matches
           by default, we accept everything that is
           returned by BLAST. The sub will be passed two 
           parameters: 
           	$inputseq : the input sequence
           	$hit : the BLAST hit
           It must return:
              1 to reject the hit
              0 to accept it
 $parameters : additional parameters to pass to BLAST

 Returns:
 An arrayref containing hash entries like:
 {
	blasthit => Bio::...::Hit,
	source_sequence => Bio::Seq,
 }

=cut

sub get_n_best_hits {
	my $sdb     = shift;
	my $seq     = shift;
	my $top_x   = shift;
	my $exclude = shift || sub {
		## my ($inputseq, $hit) = @_;
	};
	my $parameters = shift || {};

	if ( !defined($top_x) ) {
		$top_x = 250;    ## BLAST default
	}

	## in case we have a function that excludes stuff,
	## we might want a few more than top_x
	my $num_al = $top_x + 5;
	## Don't recompute all alignments for small changes
	## of top_x
	##
	if($num_al < 20) {
		$num_al = 20;
	}
	$parameters->{"-num_alignments"} = $num_al;

	my $results =
	  $sdb->search( $seq, $parameters );
	my @top_hits_sec = ();

	foreach my $result (@$results) {
		foreach my $hit ( @{ $result->{hits} } ) {
			if ( $hit->{rank} <= $top_x ) {
				if ( $exclude->( $seq, $hit ) ) {
					++$top_x;
					next;
				}
				push @top_hits_sec, $hit;
			}
		}
	}

	return \@top_hits_sec;
}

1;
