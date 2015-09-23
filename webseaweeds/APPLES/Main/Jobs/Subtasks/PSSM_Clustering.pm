#!/usr/bin/perl

=head1 PSSM Clustering using the PSSM Tool and R

=cut

use MooseX::Declare;

class Jobs::Subtasks::PSSM_Clustering extends Jobs::Job {
	use Runtime;
	use Configuration::AppleSeeds;

	use Scalar::Util qw(blessed);
	use JSON;
    use File::Temp qw(tempfile);

	require Serialization::Serializable_Array;
	require Datatypes::Motifs::PSSM;

	## the input PSSM set
	has 'pssms' => (
		is => 'rw',
		isa => 'Datatypes::Motifs::PSSM_Set'
	);

    ## Distance metric (hellinger or Kullbackleibler)
    has 'metric' => (
        is => 'rw',
        isa => 'Str',
        default  => sub { return "hellinger"; },
        required => 1,
		documentation => 'The distance metric to use.|Supported: hellinger or kullback'
    );

    ## Cutoff height
    has 'cutoff_height' => (
        is => 'rw',
        isa => 'Num',
        default  => sub { return 1.5; },
        required => 1,
		documentation => 'The tree cutoff height.',
    );

=head2 Run Clustering

 Parameters:
 None

 Returns:
 a Serializable_Array of PSSM objects

=cut

	method _run () {
		my $pssmtool_executable = find_executable("pssmtool");

		if ( !-e $pssmtool_executable ) {
			die "PSSM Tool executable $pssmtool_executable not found.";
		}

        my %matrices_by_name = ();

        ## write matrices to file
		my @matrices = @{ $self->pssms->pssms() };

		die "Must use more than 1 PSSM for clustering."
			if scalar @matrices <= 1;

        my ($fh, $filename) = tempfile();
        foreach my $matrix ( @matrices ) {
            $matrices_by_name{$matrix->name()} = $matrix;
            print $fh to_json (Serialization::Serializable::to_hash($matrix), {allow_blessed => 1}) . "\n";
        }
        close ($fh);

        ## Run distances and info tool

        my $method = -1;

        $method = 0 if $self->metric =~ m/hellinger/i;
        $method = 1 if $self->metric =~ m/kullback/i;

        die "Unknown distance metric for clustering."
            if $method < 0;

        my $dists = `$pssmtool_executable -a distances -f $filename -t $method`;

        my ($fh2, $filename2) = tempfile ();
        print $fh2 $dists;
        close($fh2);

        my %mat_info = ();

        my $info = `$pssmtool_executable -a info -f $filename`;
        my @inf_lines = split/[\r\n]+/, $info;
        foreach my $il (@inf_lines) {
            my @ilt = split /\t/, $il;
            next if scalar @ilt < 2;
            die "Invalid output from PSSM info command" if scalar @ilt != 3;
            $mat_info{$ilt[0]} = {
                length => $ilt[1],
                entropy => $ilt[2],
            };
        }

        my $R_filename = $filename2;
        $R_filename =~ s/\\/\//g;
        my $cutoff = $self->cutoff_height();

        eval_R(<<END
data = read.delim('$R_filename')
x = as.matrix(data[,(1:dim(data)[2]-1)])
clrin = cutree( hclust(as.dist(x), method="ward"), h=$cutoff )
write.table(as.matrix(clrin), '$R_filename', sep="\t", quote=FALSE)
END
);
        my %mats = ();
        open F, "<", $filename2;
        while (<F>) {
            s/(.*)[\r\n]+/$1/;
            my $l = $_;
            my @vals = split /\t/, $l;
            next if scalar @vals < 2;

            my $name = $vals[0];
            my $cluster = $vals[1];
            my $info = $mat_info{$name};
            my $m = $matrices_by_name{$name};

            die "No cluster for matrix $name" if (!defined $cluster);
            die "No info for matrix $name" if (!defined $info);

            if (exists ($mats{$cluster}) ) {
                my $on = $mats{$cluster}->name();
                my $oni = $mat_info{$on};
                die "No info for matrix $on" if (!defined $oni);

                my $q1 = $oni->{entropy} / $oni->{length};
                my $q2 = $info->{entropy} / $info->{length};

                if ($q1 < $q2) {
                    $mats{$cluster} = $m;
                }
            } else {
                $mats{$cluster} = $m;
            }
        }
        close F;
        return Serialization::Serializable_Array->new ( values %mats );
    }
}
