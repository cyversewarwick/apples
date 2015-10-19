### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

=head1 Sequence_Database implementation for Ensembl

=cut

package Sequences::Database::Sequence::Ensembl;

use strict;
use feature ":5.10";

use Moose;

use Runtime;
use Carp;

use Bio::EnsEMBL::Registry;
use List::Util qw(min max);
use Data::Dumper;

use Configuration::AppleSeeds;
use Sequences::Genomic_Sequence;
use Sequences::Annotation;
use Sequences::Database::Relative_Location;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);


with "Sequences::Database::Sequence::Ensembl::Metadata";
with "Sequences::Database::Metadata";

our $registry = 'Bio::EnsEMBL::Registry';

use Sequences::Database::Sequence;
our @ISA = qw(Sequences::Database::Sequence);

## we store registry objects globally for each type of configuration
our %registries = ();

=head2 constructor.

 Parameters:
 $class : Sequence_Database_Ensembl
 $registry_location : the Ensembl registry location, being one of:
 	ensembl, ensemblgenomes, or local

=cut

sub new {
	my $class = shift;
	my $registry_location = shift || "ensembl";

	my $self = {
		registry          => undef,
		registry_location => $registry_location,
	};
	bless $self, $class;
	return $self;
}

=head2 Initialize the registry when a sequence is accessed

 Parameters:
 $self : a Database

=cut

sub _initialize_registry {
	my $self = shift;
	our (%registries);

	my $registry          = 'Bio::EnsEMBL::Registry';
	my $registry_location = $self->{registry_location};

	## check if we have to create this particular registry
	unless ( exists( $registries{$registry_location} )
		&& defined( $registries{$registry_location} ) )
	{
#				-host => 'ensembldb.ensembl.org',
#				-user => 'anonymous'
		if ( $registry_location eq 'ensembl' ) {
			$registry->load_registry_from_db(
				-host => get_config_key('ensembl_db_host') || 'ensembldb.ensembl.org',
				-user => get_config_key('ensembl_db_user') || 'anonymous',
				-pass => get_config_key('ensembl_db_pass') || undef,
				-port => get_config_key('ensembl_db_port') || 5306,

				  #-verbose => 1 # optional (gives verbose output, useful for
				  # checking what DBs are loaded)
			);
		} elsif ( $registry_location eq 'ensemblgenomes' ) {
			$registry->load_registry_from_db(
				# -host         => 'mysql.ebi.ac.uk',
				# -port         => 4157,
				# -user         => 'anonymous',
				# -wait_timeout => 86400,
				# -host         => 'genome-mysql.cse.ucsc.edu',
				# -port         => 3306,
				# -user         => 'genome',
				# -wait_timeout => 86400,
				-host         => 'mysql-eg-publicsql.ebi.ac.uk',
				-port         => 4157,
				-user         => 'anonymous',
				-wait_timeout => 86400,
			);    # testing 24-hour timeout (default is 8 hours)
			      # ONLY loads the databases from ensemblgenomes that are
			      # available with your installed API software version
			$registry->set_disconnect_when_inactive();
		} elsif ( $registry_location eq 'local' ) {
			my $ensembl_registry_conf = config->{ensembl_registry_conf};
			# print Dumper(config);
			# confess "\nEnsembl.pm 113\n";
            #print Dumper($ensembl_registry_conf);
            #die();
			unless ( defined($ensembl_registry_conf)
				&& -e $ensembl_registry_conf )
			{
				confess(
"Ensembl local registry configuration $ensembl_registry_conf not found."
				);
			}
			$registry->load_all( $ensembl_registry_conf, 1 );

		} else {
			confess('unknown ensembl registry location $registry_location\n');
		}
		$registries{$registry_location} = $registry;
	}

	$self->{registry} = $registries{$registry_location};
}

=head2 Get the ensembl registry object

...initialize it if necessary.

 Parameters:
 $self : a Database

 Returns:
 the registry

=cut

sub registry {
	my $self = shift;
	if ( !defined( $self->{registry} ) ) {
		$self->_initialize_registry();
	}
	return $self->{registry};
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

	my $registry = $self->registry;

	my $ident = $location->identifier;

	# find out what we're actually about to retrieve
	my ( $species, $object_type, $db_type );

	( $species, $object_type, $db_type ) =
	  $self->find_species_and_object_type($ident);

	info("Getting data for $ident in $species\n");
	print "\nEnsemble.pm line 182. _get_sequence_by_location, call get_slice.\n";
	my $slice = $self->get_slice( $ident, $location );
	print "\nEnsemble.pm line 184. _get_sequence_by_location, returned from get_slice.\n";
	## get the gene + the upstream/downstream bits (we might get more actually, but we'll
	## add annotations.)
	if ( !defined($slice) ) {
		confess("Cannot create slice for identifier: "
			  . $ident
			  . " and query "
			  . Dumper($location) );
	}
    
	my $gs = Sequences::Genomic_Sequence->new(
		-id              => $ident,
		-species         => $species,
		-five_prime_pos  => $slice->start(),
		-three_prime_pos => $slice->end(),
		-seq             => $slice->seq(),
		-strand          => $slice->strand(),
		-source          => { db => $self, location => $location },
	);

	my $repeat_masked_slice = $slice->get_repeatmasked_seq();

	if ( defined $repeat_masked_slice ) {
		$gs->{masked_sequence} = $repeat_masked_slice->seq();
	} else {
		$gs->{masked_sequence} = $slice->seq();
	}
    
	## expand a little to get neighbours
	$slice->expand (500, 500);

	$self->_annotate_slice_repeats( $gs, $slice );

	my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
	$self->_annotate_slice_genes( $gs, $slice, $gene_adaptor );

	if ( $object_type =~ m/Gene/i ) {

		## the id might not be a gene, so we can return here if it isn't
		if ( !defined $gene_adaptor ) {
			confess("Gene $ident could not be retrieved as ensembl gene.");
		}

		my $ggene = $gene_adaptor->fetch_by_stable_id($ident);
		if ( !defined($ggene) ) {
			confess("Gene $ident could not be retrieved as ensembl gene.");
		}
		$self->_annotate_slice_gene_transcripts( $gs, $ggene );

		my $gann = Sequences::Annotation->new(
			$ggene->stable_id(),
			"gene",
			$ggene->display_id(),
			$ggene->stable_id(),
			$ggene->start(),
			$ggene->end(),
			$ggene->strand() < 0 ? "negative" : "positive"
		);

		$gs->add_annotation($gann);
	}
    


	my $strand = $self->get_strand($ident);
	if ($strand != $slice->strand()) {
		$gs->reverse_strand();
	}
	print "\nEnsemble.pm line 252. _get_sequence_by_location, exiting.\n";
	return $gs;
}

=head2 Annotate all genes in a slice to a Genomic_Sequence

 Fetch all genes and add them to Genomic_Sequence object as location
 annotations


 Parameters:
 $self : a Sequence_Database_Ensembl
 $gs : a Genomic_Sequence
 $slice : an Ensembl Slice
 $gene_adaptor : gene adaptor for getting absolute positions

 Returns:
 Nothing
=cut

sub _annotate_slice_genes {
	my $self  = shift;
	my $gs    = shift;
	my $slice = shift;
	my $gene_adaptor = shift;

	## find all genes contained in the slice.
	my $registry = $self->registry;
	my $genes = $slice->get_all_Genes();
	foreach my $agene ( @{$genes} ) {
		my $ggene = $gene_adaptor->fetch_by_stable_id($agene->stable_id);
		if ($ggene) {
			debug ("have slice gene: " . $ggene->stable_id . " " . $ggene->start . " " . $ggene->end);
			my $gann = Sequences::Annotation->new(
				$ggene->stable_id(),
				"gene",
				$ggene->display_id,
				$ggene->stable_id,
				$ggene->start,
				$ggene->end,
				$ggene->strand() < 0 ? "negative" : "positive"
			);
			$gs->add_annotation($gann);
		}
	}
}

=head2 Annotate all repeats in a slice to a Genomic_Sequence

 Parameters:
 $self : a Sequence_Database_Ensembl
 $gs : a Genomic_Sequence
 $slice : an Ensembl Slice

 Returns:
 Nothing
=cut

sub _annotate_slice_repeats {
	my $self  = shift;
	my $gs    = shift;
	my $slice = shift;

	my @repeats = @{ $slice->get_all_RepeatFeatures() };

	my $rcounter = 1;
	foreach my $repeat (@repeats) {

		# Slice genes may have relative locations.
		my $start = $repeat->start();
		if ( $start < $slice->start() || $start > $slice->end() ) {
			$start = $slice->start() - 1 + $repeat->start();
		}
		my $end = $repeat->end();
		if ( $end < $slice->start() || $end > $slice->end() ) {
			$end = $slice->start() - 1 + $repeat->end();
		}

		my $gann = Sequences::Annotation->new(
			"Repeat_$rcounter",
			"repeat",
			$repeat->display_id(),
			$repeat->display_id(),
			$start,
			$end,
			$repeat->strand() < 0 ? "negative" : "positive"
		);
		$gs->add_annotation($gann);
		++$rcounter;
	}
}

=head2 Annotate all transcripts in a Gene to a Genomic_Sequence

 Fetch all genes and add them to Genomic_Sequence object as location
 annotations


 Parameters:
 $self : a Sequence_Database_Ensembl
 $gs : a Genomic_Sequence
 $gene : the gene

 Returns:
 Nothing
=cut

sub _annotate_slice_gene_transcripts {
	my $self  = shift;
	my $gs    = shift;
	my $ggene = shift;

	my $icounter = 1;
	my $ident    = $ggene->stable_id();
	my $transcripts = $ggene->get_all_Transcripts ();
	foreach my $transcript ( @{$transcripts} ) {
		my $tsid = $transcript->stable_id();
		$gs->add_annotation(
			Sequences::Annotation->new(
				$tsid,
				"transcript",
				"$ident ($tsid.)",
				$transcript->stable_id(),
				$transcript->start(),
				$transcript->end(),
				$transcript->strand() < 0 ? "negative" : "positive",
				$transcript
			)
		);

		foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
			$gs->add_annotation(
				Sequences::Annotation->new(
					$tsid . "exon_$icounter",
					"exon",
					"$ident ($tsid.) : " . $exon->display_id(),
					$exon->stable_id(),
					$exon->start(),
					$exon->end(),
					$exon->strand() < 0 ? "negative" : "positive",
					$exon
				)
			);
			++$icounter;
		}
		foreach my $intron ( @{ $transcript->get_all_Introns() } ) {
			$gs->add_annotation(
				Sequences::Annotation->new(
					$tsid . "intron_$icounter",
					"intron",
					"$ident ($tsid) : intron $icounter",
					$tsid,
					$intron->start(),
					$intron->end(),
					$intron->strand() < 0 ? "negative" : "positive",
					$intron
				)
			);
			++$icounter;
		}
	}
}

=head2 Get the name of this database

=cut

sub name {
	my $self = shift;
	return "Ensembl ($self->{registry_location})";
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

	$self->find_species_and_object_type( $loc->identifier );

	return;
}

=head2 Replacement for find_species_and_object_type in ensembl that also works
       for contigs and such.

 Parameters:
 $self : a self object
 $id   : an id/accession

 Returns:
 ( $species, $object_type, $db_type, $slice )
=cut

sub find_species_and_object_type {
	my $self = shift;
	my $id   = shift;

	my $registry = $self->registry;

	my $type    = undef;
	my $db_type = undef;
	my $species = undef;

	my $cid = bless {
		what              => 'Ensembl::find_species_and_object_type',
		id                => $id,
		registry_location => $self->{registry_location},
	  },
	  "Serialization::Serializable";

	my $result = cache->cache_get($cid);

	if ( defined($result) ) {
		$species = $result->{species};
		$db_type = $result->{db_type};
		$type    = $result->{type};
	}
    
    #Fix
    #if(defined($species))
    #{
    #    $species =~ s/_/ /g;
    #}
    
	## return cached result
	if ( defined($species) && defined($type) && defined($db_type) ) {
		cache->cache_put(
			$cid,
			bless {
				species => $species,
				type    => $type,
				db_type => $db_type,
			},
			"Serialization::Serializable"
		);

		return ( $species, $type, $db_type );
	}

	( $species, $type, $db_type ) = $registry->get_species_and_object_type($id);

	if ( !defined($species) || !defined($type) || !defined($db_type) ) {

		my @dba = @{ $registry->get_all_DBAdaptors() };

		$db_type = "Core";
		foreach my $d (@dba) {
			next if $d->species() eq 'multi';

			my $slice;
			my $slice_adaptor =
			  $registry->get_adaptor( $d->species(), $db_type, "slice" );
			if ( defined($slice_adaptor) ) {
				eval {
					$slice = $slice_adaptor->fetch_by_region( 'toplevel', $id );
				};
				if ( defined($slice) ) {
					$type    = "Slice:toplevel";
					$species = $d->species();
					last;
				}
				eval {
					$slice = $slice_adaptor->fetch_by_region( 'contig', $id );
				};
				if ( defined($slice) ) {
					$type    = "Slice:contig";
					$species = $d->species();
					last;
				}
				eval {
					$slice =
					  $slice_adaptor->fetch_by_region( 'supercontig', $id );
				};
				if ( defined($slice) ) {
					$type    = "Slice:supercontig";
					$species = $d->species();
					last;
				}
				eval {
					$slice = $slice_adaptor->fetch_by_region( 'clone', $id );
				};
				if ( defined($slice) ) {
					$species = $d->species();
					$type    = "Slice:clone";
					last;
				}
			}
		}
	}

	if ( !defined($species) || !defined($type) || !defined($db_type) ) {
		confess "Cannot find information for identifier '$id'";
	}

	cache->cache_put(
		$cid,
		bless {
			species => $species,
			type    => $type,
			db_type => $db_type,
		},
		"Serialization::Serializable"
	);

	return ( $species, $type, $db_type );
}

=head2 Get the strand of an id

 Parameters:
 $self : a self object
 $id   : the object id

 Returns:
 ( 1 or -1 )

=cut

sub get_strand {
	my $self = shift;
	my $id = shift;
	my $strand = 1;

	my ( $species, $type, $db_type ) = $self->find_species_and_object_type($id);

	my $registry = $self->registry;
    
	my $oa = $registry->get_adaptor( $species, $db_type, "Gene" );
	my $gene = undef;
	if ($type =~ m/Gene/i) {
		$gene = $oa->fetch_by_stable_id ($id)
	}
	if ($type =~ m/Exon/i) {
		$gene = $oa->fetch_by_exon_stable_id ($id)
	}
	if ($type =~ m/Transcript/i) {
		$gene = $oa->fetch_by_transcript_stable_id ($id)
	}

	if ($gene) {
		$strand = $gene->strand();
	}

	return $strand;
}

=head2 Get a slice for an id

 Parameters:
 $self : a self object
 $id   : the object id

 Returns:
 a slice

=cut

sub _get_id_slice {
	my $self = shift;
	my $id   = shift;

	my ( $species, $type, $db_type ) = $self->find_species_and_object_type($id);

	my $slice         = undef;
	my $errs          = "";

    my $registry = $self->registry;
	my $slice_adaptor = $registry->get_adaptor( $species, $db_type, 'Slice' );
	given ($type) {
		when (/Gene/i) {
			eval {
				$slice = $slice_adaptor->fetch_by_gene_stable_id($id);
			};
			$errs .= $@;
		}
		when (/Transcript/i) {
			eval {
				$slice = $slice_adaptor->fetch_by_transcript_stable_id($id);
			};
			$errs .= $@;
		}
		when (/Exon/i) {
			eval {
				$slice = $slice_adaptor->fetch_by_exon_stable_id($id);
			};
			$errs .= $@;
		}
		when (/Slice/i) {
			if ( $type =~ m/Slice\:(.*)/ ) {

				my $csys = $1;
				eval { $slice = $slice_adaptor->fetch_by_region( $csys, $id ); };
			}
			$errs .= $@;
		}
	}

	if ( !defined($slice) ) {
		die "Could not get slice for id '$id' $errs";
	}

	return $slice;
}

=head2 Get a slice for an id

 Parameters:
 $self : a self object
 $id   : the object id
 $location : (optional) a Relative_Location

 Returns:
 a slice

=cut

sub get_slice {
	my $self = shift;
	my $id   = shift;
	my $loc  = shift
    || Sequences::Database::Relative_Location->new( identifier => $id );
    
	my ( $species, $type, $db_type ) = $self->find_species_and_object_type($id);
    
    #Get the gene sequence
	my $strand = $self->get_strand ($id);
	my $slice         = $self->_get_id_slice ($id);
        
	## all relative to 5' end?
	my $orig_len = $slice->length();
    
	## if the strand of our object is different from
	## the slice strand, we mirror our relative location
    ## (i.e. an offset of 1 becomes an offset of -1
	my $strand_was_mirrored = 0;
	if ($slice->strand() != $strand) {
		$loc->mirror_strand();
		$strand_was_mirrored = 1;
	}
    
    ## Get the location start + end
	my ( $offset, $length, $stop_at_neighbouring_gene ) =
    ( $loc->offset, $loc->length, $loc->stop_at_neighbours );
    
	my $start;
	my $end;
    
	if ( $loc->anchor eq "5' end" ) {
		## go forward from 5' end
		$start = $slice->start + $offset;
		$end   = $slice->start + $offset + $length;
        
		## all relative to 3' end?
	} elsif ( $loc->anchor eq "3' end" ) {
		$start = $slice->end + $offset;
		$end   = $slice->end + $offset + $length;
	} else {
		die "Unsupported sequence anchor : " . $loc->anchor;
	}
    
	my $u = 0;
	my $d = 0;
    
    #Flip starts and ends
	if ($start > $end) {
		my $tmp = $start;
		$start = $end;
		$end = $tmp;
	}
    
    #upstream and downstream
    #how far do we expand or cut upstream/downstream?
	$u = $slice->start () - $start;
	$d = $end - $slice->end () ;
    
	my $inset_u = 0;
	my $inset_d = 0;
    
	if ( $u < 0 ) {
		$inset_u = -$u;
		$u = 0;
	}
	if ( $d < 0 ) {
		$inset_d = -$d;
		$d = 0;
	}
    
	if ( $u > 0 || $d > 0 ) {
        $slice = $slice->expand( $u, $d );
		if (! defined ($slice) ) {
			die "Could not expand slice! Most likely, there is not enough data on Ensembl to do so.";
		}
	}
    
	if ($inset_u > 0 || $inset_d > 0) {
		my $ostart = $slice->start;
		my $oend = $slice->end;
		$slice = $slice->sub_Slice( 1 + $inset_u, $slice->length() - $inset_d );
		if (! defined ($slice) ) {
			die "Could not sub_Slice to fix length: $start / $ostart  $end / $oend  $inset_u $inset_d. ";
		}
	}
    
    #Resize because of neighbours
	if ($stop_at_neighbouring_gene) {
		my $md    = $self->_get_meta_data($id);
		my $egenes = $slice->get_all_Genes();
        
		my @genes = ();
		my $start = undef;
		my $end   = undef;
        
		foreach my $agene ( @{$egenes} ) {
			## ignore the current gene (if there is one)
			if (scalar(grep { $_ eq $agene->stable_id } @{ $md->{nearest_genes} }) > 0)
			{
				if ( !defined($start) || $agene->start() < $start ) {
					$start = $agene->start();
				}
				if ( !defined($start) || $agene->end() < $start ) {
					$start = $agene->end();
				}
				if ( !defined($end) || $agene->start() > $end ) {
					$end = $agene->start();
				}
				if ( !defined($end) || $agene->end() > $end ) {
					$end = $agene->end();
				}
			} else {
				push @genes,
                {
					id    => $agene->stable_id(),
					start => min( $agene->start(), $agene->end() ),
					end   => max( $agene->start(), $agene->end() ),
                };
			}
		}
        
		my $nearest_up   = undef;
		my $nearest_down = undef;
        
		foreach my $g (@genes) {
			debug ("have neighbour: $g->{id} $g->{start} $g->{end}");
            if ( $g->{end} < $start ) {
                if ( !defined($nearest_up) ) {
                    $nearest_up = $g->{end};
                } else {
                    my $dist1 = $start - $nearest_up;
                    my $dist2 = $start - $g->{end};
                    if ( $dist2 < $dist1 ) {
                        $nearest_up = $g->{end};
                    }
                }
            }
            if ( $g->{start} > $end ) {
                if ( !defined($nearest_down) ) {
                    $nearest_down = $g->{start};
                } else {
                    my $dist1 = $nearest_down - $end;
                    my $dist2 = $g->{start} - $end;
                    if ( $dist2 < $dist1 ) {
                        $nearest_down = $g->{start};
                    }
                }
            }
		}
		my $len = abs ( $slice->start - $slice->end);
		debug ("Length: $len strand " . $slice->strand);
		if ( defined($nearest_down) || defined($nearest_up) ) {
			if ( !defined($nearest_up) ) {
				$nearest_up = 0;
			}
			if ( !defined($nearest_down) ) {
				$nearest_down = $len+1;
			}
            
			if ($nearest_up < 0) {
				$nearest_up = 0;
			}
			if ($nearest_down > $len+1) {
				$nearest_down = $len+1;
			}
            debug ("resizing to accommodate neighbours : $nearest_up $nearest_down");
            
            if($strand_was_mirrored)
            {
                #if there is nothig to take (i.e. empty slice), make it a 1bp slice to avoid errors
                if($nearest_up+1 > $nearest_down)
                {
                    $nearest_down = $nearest_up+1;
                }
                $slice = $slice->sub_Slice( $nearest_up+1, $nearest_down);
            }
            else
            {
                #if there is nothig to take (i.e. empty slice), make it a 1bp slice to avoid errors
                if($nearest_up + 1 > $nearest_down - 1)
                {
                    $nearest_down = $nearest_up+2;
                }
                $slice = $slice->sub_Slice( $nearest_up+1, $nearest_down-1);
            }
            
			if (! defined ($slice) ) {
				die "Could not sub_Slice to fix $id neighbour gene distance : $nearest_up, $nearest_down.";
			}
		}
        
        
	}
    
    #/stop neighbours
    
	## enforce length limit
	if ( abs($length) > 0 && $slice->length > abs($length) )
    {
        if($strand_was_mirrored)
        {
            $slice = $slice->sub_Slice( 2, abs($length)+1 );
            if (! defined ($slice) ) {
                die "Could not sub_Slice to fix length.";
            }
        }
        else
        {
            $slice = $slice->sub_Slice( 1, abs($length) );
            if (! defined ($slice) ) {
                die "Could not sub_Slice to fix length.";
            }
        }
    }
    
	debug ("Resized slice for " . $loc->identifier . " location " . $slice->start . ":" . $slice->end);
    
    #------ Move this function -------#
    #Is this slice inside a gene?
    #[ND] - Addition
    #If so, return a slice of 1 bp
    if ($stop_at_neighbouring_gene) {
        my $egenes = $slice->get_all_Genes();
        foreach my $agene ( @{$egenes} ) {
            if( $agene->start() <= 1 && $agene->end() >= $slice->length() )
            {
                $slice = $slice->sub_Slice(1, 1);
            }
        }
    }
    debug "\nIs inside gene. Resized to 1bp";
    #------- Move --------#
    
	if ($strand_was_mirrored) {
		$loc->mirror_strand();
	}
    
    return $slice;
}


=head2 Serialization override

=cut

sub to_hash {
	use Bio::EnsEMBL::ApiVersion;
	my $self = shift;
	return {
		'registry_location'   => $self->{registry_location},
		'ensembl_api_version' => software_version(),
	};
}

=head2 Get a list of all accessions in this database

 Parameters:
 $self : a self object
 $species : the species

 Returns:
 an ARRAYREF containing all the accessions

=cut

sub get_all_accessions {
	my $self = shift;
	my $species = shift;

	my $registry = $self->registry;

    my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );

    
    my @geneIDlist  = @{ $gene_adaptor->list_stable_ids() };
    return \@geneIDlist;
}

#Gets the canonical transcript sequence of a gene by accession
sub get_canonical_transcript_sequence_by_accession {
    my $self = shift;
	my $species = shift;
    my $accession = shift;
    
    my $registry = $self->registry;
    
    my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
    
    my $gene = $gene_adaptor->fetch_by_stable_id($accession);
    
    my $final_seq = "";

    #Go through all exons
    foreach my $exon ( @{ $gene->canonical_transcript->get_all_Exons } )
    {
        
        #If this exon contains a coding sequence
        if(defined($exon->coding_region_start($gene->canonical_transcript)))
        {
            my $seq;
            
            if($exon->strand == "1")
            {
                $seq = $exon->seq->subseq( ($exon->coding_region_start($gene->canonical_transcript) - $exon->start) + 1, ($exon->coding_region_end($gene->canonical_transcript) - $exon->start) + 1 );
            }
            else
            {
                $seq = $exon->seq->subseq( ($exon->end - $exon->coding_region_end($gene->canonical_transcript)) + 1, ($exon->end - $exon->coding_region_start($gene->canonical_transcript)) + 1);
            }
            
            $final_seq .= $seq;
        }
    }
    #print "\nGene strand: " . $gene->strand;
    #print "\n$final_seq";

    return $final_seq;
}

sub get_strand_by_accession {
    my $self = shift;
	my $species = shift;
    my $accession = shift;
    
    my $registry = $self->registry;
    
    my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
    
    my $gene = $gene_adaptor->fetch_by_stable_id($accession);

    print "\n" . $gene->strand . " : " . $gene->start . " : " . $gene->end;
}

#Gets the sequence of a gene by accession
sub get_gene_sequence_by_accession {
	my $self = shift;
	my $species = shift;
    my $accession = shift;
    
	my $registry = $self->registry;
    # print "\n Get Adaptor...";
    my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
    # print "\n Get Gene...";
    my $gene = $gene_adaptor->fetch_by_stable_id($accession);
    # print Dumper($gene);
    print "\nEnsemble.pm line 1024. get_gene_sequence_by_accession, Accession: " . $accession . "\n End:" . $gene->{end} . "\n Start:" . $gene->{start} . "\n";
    my $gene_length = $gene->{end} - $gene->{start};
    print "\nEnsemble.pm line 1026. get_gene_sequence_by_accession, Get Location...\n";
    my $location =
    Sequences::Database::Relative_Location->new( identifier => $accession ,
                                                length => $gene_length );
    print "\nEnsemble.pm line 1030. get_gene_sequence_by_accession, exiting, call get_sequence_by_location.\n";
	return $self->get_sequence_by_location($location);

}

#Fill this in and we are golden
sub get_gene_information_by_accession   {
    my $self = shift;
	my $species = shift;
    my $accession = shift;
    
	my $registry = $self->registry;
    
    my $gene_adaptor = $registry->get_adaptor( $species, 'Core', 'Gene' );
    
    my $gene = $gene_adaptor->fetch_by_stable_id($accession);
    
    my $strand = "";
    if($gene->strand == -1)
    {
        $strand = "negative";
    }
    else
    {
        $strand = "positive";
    }
    return ($gene->seq_region_name, $gene->coord_system_name);
}

=head2 Serialization override

=cut

sub serialization_restore {
	my $oldself = shift;

	my $self = Sequences::Database::Sequence::Ensembl->new(
		$oldself->{registry_location} );
	%$oldself = %$self;
}

1;
