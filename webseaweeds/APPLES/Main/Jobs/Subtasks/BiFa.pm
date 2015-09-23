#!/usr/bin/perl

=head1 Sequence Retrieval Job

This job retrieves a set of sequences from sequence databases into a dataset.

=cut
use MooseX::Declare;

class Jobs::Subtasks::BiFa extends Jobs::Job {

	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	use Motifs::BiFa_Server_Interface;
    
    use Data::Dumper;


	## Sequence length limit.
	## this is to limit the length of sequences
	## submitted. You can change this at initialisation
	## for using longer sequences, mainly it serves the
	## purpose of preventing web users from submitting
	## huge jobs.
	has 'length_limit' => (
		is       => "rw",
		isa      => "Int",
		default  => sub { return 20000; },
		required => 1,
	);

	## Sequence to analyze: must be set
	has 'sequence' => (
		is       => "rw",
		isa      => "Str",
		required => 1,
	);

	## a title to show at the back of the SVG file
	has 'title' => (
		is       => "rw",
		isa      => "Str",
		default  => sub { return "BiFa analysis"; },
		required => 1,
	);

	## related sequences to filter output
	has 'phylogeneticly_related_sequences' => (
		is       => "rw",
		isa      => "ArrayRef[Str]",
		default  => sub { return []; },
		required => 1,
	);

	## BiFa p-value threshold
	has 'threshold' => (
		is       => "rw",
		isa      => "Num",
		default  => sub { return 0.0001; },
		required => 1,
	);

	## BiFa phylo p-value threshold
	has 'phylo_threshold' => (
		is       => "rw",
		isa      => "Num",
		default  => sub { return 0.05; },
		required => 1,
	);

	## BiFa's 'old algorithm' option. whatever that might mean.
	has 'use_old_algorithm' => (
		is       => "rw",
		isa      => "Bool",
		default  => sub { return 0; },
		required => 1,
	);

	## ALG_BAYESIAN or ALG_OTT
	has 'bifa_algorithm' => (
		is       => "rw",
		isa      => "Int",
		default  => sub { return Motifs::BiFa_Server_Interface::ALG_BAYESIAN; },
		required => 1,
	);

	## selection of PSSM sets to pass to BiFa.
	## from the original BiFa test:
	##	Alternative for specific pssm sets
	##	my @pssmSets = ["transfac"];
	has 'pssmsets' => (
		is       => "rw",
		isa      => "ArrayRef[Str]",
		default  => sub { return []; },
		required => 1,
	);

	## pssm filter. apparently, "." means no filtering
	has 'pssm_filter' => (
		is       => "rw",
		isa      => "Str",
		default  => sub { return "."; },
		required => 1,
	);

	## selection of species specific PSSMS ( PVI or any combination of these )
	has 'species' => (
		is       => "rw",
		isa      => "Str",
		default  => sub { return "PVI"; },
		required => 1,
	);

	## use consensus sequences? We use 1 as default, like the BiFa example script
	has 'use_consensus' => (
		is       => "rw",
		isa      => "Bool",
		default  => sub { return 1; },
		required => 1,
	);

	## dummy to make sure we use the same data version through all bifa jobs and also save it for caching
	## this should be set upon creation, see Jobs::BiFa_Job
	has 'data_and_tf_version' => (
		is      => "rw",
		isa     => "Str",
		default => sub {
			my $soapIf = Motifs::BiFa_Server_Interface->new;
			if ( !$soapIf->init() ) {
				die "Error connecting to BiFa Server.";
			}

			#	Get the transfac version
			my $tv = $soapIf->TransfacVersion();

			#	Get the custom PSSM version
			my $cv = $soapIf->CustomPssmVersion();

			# Construct the data version which is used as a parameter
			# in subsequent requests
			my $dv = $tv . '.' . $cv;
            
			return $dv;
            
		},
	);

	## if you set this to one, the output will contain a
	## svg item, containing the BiFa SVG output.
	## This might become quite large, so by default
	## we don't use it.
	has 'make_svg' => (
		is       => "rw",
		isa      => "Bool",
		default  => sub { return 0; },
		required => 1,
	);


=head2 Overloaded validation

=cut

	method validate () {
		$self->SUPER::validate();

		if (length ( $self->sequence ) > $self->length_limit ) {
			die "Input sequence is too long. Limit is " . $self->length_limit;
		}
		
		foreach my $ps (@{$self->phylogeneticly_related_sequences()}) {
			if (length ( $ps ) > $self->length_limit ) {
				die "Input phylo-sequence is too long. Limit is " . $self->length_limit;
			}
		}			
	}

=head2 Run BiFa analysis job

 Parameters:
 None

 Returns:
 A Serializable object with one key 'hits' containing the BiFa hits.

=cut

	method _run () {
		## get the sequence, store it in the cache.
		my $result = bless {}, "Serialization::Serializable";

		#	Initialise the soap interface to the server

		my $algorithm = Motifs::BiFa_Server_Interface::ALG_BAYESIAN;

		my $soapIf = Motifs::BiFa_Server_Interface->new;
		if ( !$soapIf->init() ) {
			die "Error connecting to BiFa Server.";
		}
        

		#	Get the transfac version
		my $tv = $soapIf->TransfacVersion();

		#	Get the custom PSSM version
		my $cv = $soapIf->CustomPssmVersion();

		# Construct the data version which is used as a parameter
		# in subsequent requests
		my $dv = $tv . '.' . $cv;

		if (   $self->data_and_tf_version ne ""
			&& $dv ne $self->data_and_tf_version )
		{
			die "Data version mismatch : $dv vs " . $self->data_and_tf_version;
		}

		my $cumulative = 0;
		if ( $self->bifa_algorithm() == Motifs::BiFa_Server_Interface::ALG_OTT )
		{
			$cumulative = 1;
		}

		#	get the information out as an SVG file
		if ( $self->make_svg ) {
			$result->{svg} = $soapIf->bifa(
				$self->sequence(),
				$self->threshold(),
				$self->title(),
				$self->bifa_algorithm(),
				0,    ## showLabels = 0
				$self->phylogeneticly_related_sequences(),
				$self->use_consensus(),
				$self->species(),
				$self->phylo_threshold(),
				$self->pssm_filter(),
				$self->use_old_algorithm(),
				$cumulative,
				$self->pssmsets(),
			);
		}
        
        my $specs = "PVINFB";
        my @pssmNames = $soapIf->PssmSetNames();
        push(@pssmNames, "transfac");
        

		#	get the information out as a series of hits
		my @hits = $soapIf->bifaHits(
			$self->sequence(),
			$self->threshold(),
			$self->bifa_algorithm(),
			$self->phylogeneticly_related_sequences(),
			$self->use_consensus(),
			$specs,
			$self->phylo_threshold(),
			$self->pssm_filter(),
			$cumulative,
			$dv,
			\@pssmNames
		);


		$result->{hits} = \@hits;

		return $result;
	}

};
