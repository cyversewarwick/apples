#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut

use MooseX::Declare;

class Jobs::Syntheny_Job extends Jobs::Job {

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

    use Datatypes::Concepts::Gene;

	use Datatypes::Sequence::Search;
	use Datatypes::Sequence::Location;
    use Datatypes::Sequence::Set;

    use Jobs::Subtasks::Sequence_Retrieval;

    use Orthology::Database::Ensembl;
    use Links::Links_Database;

	## used by run to specify the
	## place where the sequences go
	has 'sequenceset' => (
							   is            => 'rw',
							   isa           => 'Str',
							   documentation => 'O1. Output Sequence Set Name',
							   default       => sub { return "Sequences"; },
	);

   	## used by run to specify the
	## place where the sequences go
	has 'orthologyset' => (
							   is            => 'rw',
							   isa           => 'Str',
							   documentation => 'O2. Output Orthology Dataset Name',
							   default       => sub { return "Orthology Information"; },
	);

    ## Search query string
	has 'search' => (
					  is            => 'rw',
					  isa           => 'Sequence_Search',
					  documentation => '1. Sequence IDs|You should specify a larger area that actually includes a few neighbours here.',
					  default       => sub { return "" },
	);

    ## Specify the compara database to use
    has 'species' => (
        is => 'rw',
        isa => 'ArrayRef[Str]',
        documentation => "2. Which species to use",
        default => sub { return []; },
    );

    ## Specify the compara database to use
    has 'which_compara' => (
        is => 'rw',
        isa => 'Str',
        documentation => "3. Location of compara database|Use 'ensembl', 'ensemblgenomes' or 'local'",
        default => sub { return "ensembl"; },
    );


=head2 Overloaded validation method

=cut

	method validate () {
		if ( scalar @{ $self->species } < 1 ) {
			die "You need to select at least one species to search for orthologues in.";
		}

		my $s         = Datatypes::Sequence::Search->new;
		my $locations = $s->parse_search( $self->search );

        if (scalar @$locations <= 0) {
            die "You must specify at least one sequence location.";
        }

		$self->SUPER::validate();
	}

=head2 Run retrieval job

 Parameters:
 None

 Returns:
 A serializable array of Genomic_Sequence objects

=cut

	method _run () {
		## get the sequence, store it in the cache.

		my $s         = Datatypes::Sequence::Search->new;
		my $locations = $s->parse_search( $self->search );
        my $species = $self->species;

        my @seqs = ();

        my $compara = Orthology::Database::Ensembl->new($self->which_compara);

        my $links = Links::Links_Database->new;
        my %accs_done = ();

        while (scalar @$locations > 0) {
            my $loc = shift @$locations;
            debug ("Processing " . $loc->accession);

            $accs_done{$loc->accession} = 1;

            my $n1 = $links->get_or_add_node(
                Datatypes::Concepts::Gene->new_gene(
                    $loc->db,
                    $loc->accession
                )
            );

            my $rseqs = $loc->get_sequence;
            foreach my $seq (@$rseqs) {
                push @seqs, $seq;

                my $seq_genes = $seq->get_annotations_by_type('gene');

                foreach my $sgene (@$seq_genes) {
                    my $sacc = $sgene->{sourceid};
                    $sgene->{metadata} = $compara->get_meta_data($sacc);
                    push @$locations, Datatypes::Sequence::Location->new_simple_location ( $loc->db, $sacc )
                       unless defined ($accs_done{$sacc});
                }

                my $others = $compara->fetch_orthologous_ids($n1, $species);
                foreach my $o (@$others) {
                    $links->get_or_add_node(
                        Datatypes::Concepts::Gene->new_gene(
                            $self->which_compara,
                            $o->accession
                        )
                    )->link_bidirectional($n1)->data("compara", 1);

                    push @$locations, Datatypes::Sequence::Location->new_simple_location ( $loc->db, $o->accession )
                       unless defined ($accs_done{$o->accession});
                }
            }
        }

        ## Save and merge output.
        my $sset = Datatypes::Sequence::Set->new (sequences => \@seqs);

        $self->put_dataitem( $sset, $self->sequenceset );
        $self->put_dataitem( $links, $self->orthologyset );

	  	$self->merge_dataitem_versions( $self->sequenceset,
	  	 	\&Datatypes::Sequence::Set::merged_sets
		);
    	$self->merge_dataitem_versions( $self->orthologyset,
		  \&Links::Links_Database::merge_nodes );

        ## Return nothing (we don't cache and just output)
        my $res = {};
        bless $res, "Serialization::Serializable";

		return $res;
	}

=head2 Expiry time (disable cache)

(all caching is handled by the created BiFa subtasks)

=cut

	method expiry_time() {
		return 0;
	}

};
