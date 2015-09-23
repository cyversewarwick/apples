#!/usr/bin/perl

=head1 PSSM Clustering Job

Clusters PSSMs in selected sets.

=cut

use MooseX::Declare;

class Jobs::PSSM_Clustering_Job extends Jobs::Job {
	use Runtime;

    require Datatypes::Motifs::PSSM_Set;
	require Jobs::Subtasks::PSSM_Clustering;

   	## used by postprocess to specify the
	## place where the retrieved PSSMs go
	## this parameter is not cached since
	## it is only used in postprocess
	has 'nc__pssmset' => (
		is            => 'rw',
		isa           => 'Str',
		documentation => '1. Output PSSM set dataitem name',
		default       => sub { return "PSSMs"; },
	);

    ## Distance metric (hellinger or Kullbackleibler)
    has 'metric' => (
        is => 'rw',
        isa => 'Str',
        default  => sub { return "hellinger"; },
        required => 1,
		documentation => '2. The distance metric to use.|Supported: hellinger or kullback'
    );

    ## Cutoff height
    has 'cutoff_height' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 1.5; },
        required => 1,
		documentation => '3. The tree cutoff height.',
    );


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

        my $cluster_matrices = Jobs::Subtasks::PSSM_Clustering->new(
            pssms => $pset,
            metric => $self->metric,
            cutoff_height => $self->cutoff_height,
            nc__pssmset => $self->nc__pssmset
        );
        return $cluster_matrices->run();
    }


=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		$self->put_dataitem(
			Datatypes::Motifs::PSSM_Set->new (pssms => $result->data()),
			$self->nc__pssmset()
		);
		return $result;
	}

=head2 Expiry time (disable cache)

(all caching is handled by the created BiFa subtasks)

=cut

	method expiry_time() {
		return 0;
	}

    with "Jobs::Roles::Dataitem_Selection_Job";
};