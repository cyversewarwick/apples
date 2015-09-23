
=head1 Class Search_RBH

Class for searching for reciprocal best hit orthologues

=cut

package Orthology::Search::RBH;

use strict;

use Links::Links_Database;

use Datatypes::Sequence::Location;
use Orthology::Similarity_Link;
use Orthology::RBH_Link;

use Runtime;
use Carp;
use Data::Dumper;

=head2 Constructor

 Parameters:
 $class : Search_RBH
 $db1 : first database to search
 $db2 : second database to search
 $cn  : (optional) the candidate network 

=cut

sub new {
	my ( $class, $db1, $db2 ) = @_;
	my $self = {
				 db1 => $db1,
				 db2 => $db2,
	};
	bless $self, $class;
	return $self;
}

=head2 Perform a reciprocal best hit search

 Parameters:
 $self  : a RBH search object
 $acc   : the accession in the primary database to match 
 $top_x : the number of top hits to include
 $rec   : the level of recursive searches to perform
          when building the candidate network (default is 2)
 
 Returns:
 { 
 	rbhs => <a list of reciprocal best hits>
 	network => a Links::Links_Database with a candidate network.
 } 
 
 
=cut

sub search {
	my $self  = shift;
	my $acc   = shift;
	my $top_x = shift || 1;
	my $rec   = shift || 2;

    #Work on this
	my $candidate_network = Links::Links_Database->new;

	my $pdb = get_search_database( $self->{db1} );
	my $sdb = get_search_database( $self->{db2} );

	my %acc_to_node = ();
	my $nds         = $candidate_network->get_nodes();
    
	foreach my $nd (@$nds) {
		if ( UNIVERSAL::isa( $nd, "Datatypes::Sequence::Location" ) ) {
			$acc_to_node{ $nd->accession } = $nd;
		}
	}

	my $start_node = $acc_to_node{$acc}
	  || Datatypes::Sequence::Location->new_simple_location( $self->{db1},
															  $acc );
	$acc_to_node{$acc} = $start_node;
    
    

	my @nodes_to_search = (
		{
		   node => $start_node,
		   ## FIXME let's hope at this point that each ID only returns
		   ## a single sequence. otherwise, we might want to add
		   ## more nodes here
		   seq      => $pdb->get_sequence($acc)->[0],
		   dbn      => $self->{db2},
		   reclevel => 1,
		},
	);

    #Search nodes recursively
	while ( scalar @nodes_to_search > 0 ) {
		debug( "Number of nodes to search: " . ( scalar @nodes_to_search ) );
		my $search = shift @nodes_to_search;

        #Ths is the top hits
		my $ths = get_search_database( $search->{dbn} )->get_n_best_hits(
			$search->{seq},
			$top_x,
			sub {
				## here we test if we can accept a hit
				my ( $seq, $hit ) = @_;

				## we rule out
				## the same accession coming from the same database
				if ( $seq->{source}->{db}
					    ->equals( get_search_database( $search->{dbn} ) )
					 && $seq->id() eq $hit->{accession} )
				{
					return 1;
				}

				return 0;
			},
			get_search_database_info( $search->{dbn} )->{parameters},
		);
        
        #print "\n" . $search->{"seq"}->{source}->{location}->{identifier};

        #Now we add each hit to the nodes to search, unless recursion limit is met
		# see how/if they map back
		foreach my $th (@$ths) {
			## create network link
            #print "\nHit: " . $th->{hit_name} . ": " . $th->{significance};
			my $th_link = Orthology::Similarity_Link->new( hit => $th );

			my $th_node = undef;
			if ( !exists( $acc_to_node{ $th->{accession} } ) ) {
				$th_node =
				  Datatypes::Sequence::Location
				  ->new_simple_location(
												$search->{dbn}, $th->{accession}
				  );
				$acc_to_node{ $th->{accession} } = $th_node;
			} else {
				$th_node = $acc_to_node{ $th->{accession} };
			}

			my $hit_edge = $search->{node}->link_to($th_node);
			$hit_edge->data( 'BLAST', $th_link );

			push @nodes_to_search,
			  {
				node => $th_node,
				seq  => get_search_database( $search->{dbn} )
				  ->get_sequence( $th->{accession} )->[0],
				dbn => ( $search->{dbn} eq $self->{db1} )
				? $self->{db2}
				: $self->{db1},
				reclevel => $search->{reclevel} + 1,
			  }
			  unless $search->{reclevel} > $rec;
		}
	}
    
    #print "\nTest";
    #print Dumper(keys %acc_to_node);

	my @rbhs = ();
	foreach my $node ( values %acc_to_node ) {
		$candidate_network->add_node($node);

		## filter links to only use the ones discovered in this run
		my $ndb_source = $node->db;
		next unless $ndb_source eq $self->{db1};
		my $ndb_target = $self->{db2};

		## check if we have a reciprocal best hit
		my $nrbh = $self->get_best_hit( $node, $ndb_target );

		my $bnrbh = undef;
		if ( defined($nrbh) ) {
			$bnrbh = $self->get_best_hit( $nrbh, $ndb_source );
		}

		## if the best hit maps back, it's a reciprocal best hit
		if (    defined($nrbh)
			 && defined($bnrbh)
			 && $bnrbh == $node )
		{
            #print "\nRBH!";
            #print "\nBetween: " . $node->{'location'}->{"identifier"} . " and " . $nrbh->{"location"}->{"identifier"};
            #print "\nNode before: " . $node->{"location"}->{"identifier"};
            
            
			my $w1 = $node->link_to($nrbh)->data('BLAST')->weight;
			my $w2 = $nrbh->link_to($node)->data('BLAST')->weight;

			$node->link_bidirectional($nrbh)->data(
				'RBH',
				Orthology::RBH_Link->new(
					## we use the geometric mean.
					## maybe some other way might be better?
					significance => sqrt( $w1 * $w2 ),
				),
			);

			push @rbhs, $nrbh;
		}
	}
    
	return $candidate_network;
}

=head2 sub to get the best hit given a 
Orthology::Gene_Identifer node  

 Parameters:
 $self : a self object
 $node : Orthology::Gene_Identifer
 $dbn  : a search database name the accession needs to be in

 Returns: 
 a Orthology::Gene_Identifer node that is the best
 (connected by Links::Relation with 
  Orthology::Similarity_Link node to give similarity 
  score ) match (smallest significance score)

=cut

sub get_best_hit {
	my $self = shift;
	my $node = shift;
	my $dbn  = shift;

	my $links = $node->get_links;

	my $best_score = undef;
	my $best_link  = undef;

	foreach my $l (@$links) {

		my $d = $l->data('BLAST');

		eval {
			get_search_database($dbn)->get_sequence( $d->hit->{accession} );
		};
		my $right_db = $@ ? 0 : 1;

		if (    defined($d)
			 && UNIVERSAL::isa( $d, "Orthology::Similarity_Link" )
			 && $right_db )
		{
			if ( !defined($best_score) || $d->weight < $best_score ) {
				$best_score = $d->weight;
				$best_link  = $l;
			}
		}
	}

	if ( defined($best_link) ) {
		return $best_link->target;
	} else {
		return undef;
	}
}

1;
