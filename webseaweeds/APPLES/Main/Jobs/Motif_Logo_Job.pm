#!/usr/bin/perl

=head1 PSSM Clustering Job

Clusters PSSMs in selected sets.

=cut

use MooseX::Declare;

class Jobs::Motif_Logo_Job extends Jobs::Job {
	use Runtime;

    require Serialization::Serializable_Array;

=head2 Run Clustering

 Parameters:
 None

 Returns:
 a PSSM set

=cut

	method _run() {
        my $pssms = $self->selected_objects_of_kind("Datatypes::Motifs::PSSM");
        my $pset = Datatypes::Motifs::PSSM_Set->new (
            pssms => $pssms,
        );

		if ((scalar @$pssms) <= 0) {
			die "You need to select some PSSMs for this to work.";
		}

        for my $pssm (@$pssms) {
            $pssm->get_png();
        }

        return Serialization::Serializable_Array->new( data => [] );
    }

=head2 Expiry time (disable cache)

(all caching is handled by the created BiFa subtasks)

=cut

	method expiry_time() {
		return 0;
	}

    with "Jobs::Roles::Dataitem_Selection_Job";
};