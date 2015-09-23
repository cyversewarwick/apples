#!/usr/bin/perl

use MooseX::Declare;

=head1 Compute reciprocal best BLAST hits

This job returns a network of BLAST scores, together with reciprocal
best hit information.

=cut

class Jobs::Subtasks::Orthology_Cleanup extends Jobs::Job {

	with "Jobs::Roles::Dynamic_Data_Loader";
	with "Jobs::Roles::Links_Job";

	use Runtime;

	use Datatypes::Moose_Types;
	use Links::Links_Database;

	use Jobs::Subtasks::RBH;
	require Datatypes::Sequence::Location;
	require Sequences::Database::Relative_Location;

	require Datatypes::Concepts::Gene;
	require Datatypes::Concepts::Transcript;

	require Datatypes::Roles::Locatable;

=head2 Overridden run method

=cut

	method _run () {
		my @all_ids =
		  @{ $self->list_nodes_by_type("Datatypes::Sequence::Location") };

		my %ids_to_keep = ();

		## 1. remove nodes not connected to any others via a RBH
		foreach my $n (@all_ids) {

			my @l = @{ $n->get_links };

			my $isolated = 1;

			foreach my $ln (@l) {
				if ( defined( $ln->data('RBH') ) ) {
					$isolated = 0;
				}
			}

			if ( !$isolated ) {
				$ids_to_keep{ $n->accession } = { old => $n, };
			}
		}

		## 2. Make genes/transcripts
		foreach my $id ( keys %ids_to_keep ) {
			my $old = $ids_to_keep{$id}->{old};
			my $new = undef;

			my $db = $old->db;

			## TODO there must be a less hard-coded way to do this
			if (   $db ne 'ensembl'
				&& $db ne 'ensemblgenomes'
				&& $db ne 'genbank' )
			{
				$db = Datatypes::Roles::Locatable::_guess_db_by_identifier(
					$old->accession );
			}

			my $md =
			  get_sequence_database($db)->get_meta_data( $old->accession );

			if ( defined( $md->{type} ) ) {
				if ( $md->{type} =~ m/Gene/i ) {
					$new = Datatypes::Concepts::Gene->new_gene( $db,
						$old->accession );
				} elsif ( $md->{type} =~ m/Transcript/i ) {
					$new = Datatypes::Concepts::Transcript->new_transcript( $db,
						$old->accession );
				}
			}

			if ( !defined($new) ) {
				$new =
				  Datatypes::Sequence::Location->new_simple_location( $old->db,
					$db );
			}

			$ids_to_keep{$id}->{new}  = $new;
			$ids_to_keep{$id}->{meta} = $md;
		}

		## 3. Make links
		foreach my $id ( keys %ids_to_keep ) {
			my $old = $ids_to_keep{$id}->{old};
			my $new = $ids_to_keep{$id}->{new};
			my $md  = $ids_to_keep{$id}->{meta};

			if ( defined( $md->{nearest_genes} ) ) {
				foreach my $iid ( @{ $md->{nearest_genes} } ) {
					next if $iid eq $id;
					if ( !defined( $ids_to_keep{$iid} ) ) {
						$ids_to_keep{$iid} = {
							'new' => Datatypes::Concepts::Gene->new_gene(
								$new->db, $iid
							)
						};

						$ids_to_keep{$iid}->{old} = $ids_to_keep{$iid}->{new};
						$ids_to_keep{$iid}->{meta} =
						  $ids_to_keep{$iid}->{new}->get_meta_data;
					}
					$new->link_bidirectional( $ids_to_keep{$iid}->{new} )
					  ->data( 'via', $new->db );
				}
			}

			my @l = @{ $old->get_links };

			foreach my $ln (@l) {
				my $rbh = $ln->data('RBH');
				if ( defined($rbh) ) {
					my $ids = $ln->source->accession;
					if ( $ids eq $id ) {
						$ids = $ln->target->accession;
					}

					## we should have kept all the nodes which have RBH links
					## therefore, $ids_to_keep must have a key $ids
					$new->link_bidirectional( $ids_to_keep{$ids}->{new} )
					  ->data( 'RBH', $rbh );
				}
			}
		}

		my $new_ldb = Links::Links_Database->new;
		foreach my $id ( keys %ids_to_keep ) {
			$new_ldb->add_node( $ids_to_keep{$id}->{new} );
		}

		return $new_ldb;
	}

=head2 This method will be run by the cache after the result has
been retrieved

 Parameters:
 $result : the computation result (cached or just computed)

 Returns:
 a postprocessed version of the result

=cut

	method postprocess ( Serialization::Serializable $result ) {
		$self->put_dataitem( $result, 'Orthology Search Result' );
		return $result;
	}

=head2 Expiry time (disable cache)

=cut

	method expiry_time() {
		return 0;
	}

};
