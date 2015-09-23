#!/usr/bin/perl

=head1 PSSM Retrieval Job

This jobs retrieves a set of PSSMs from BiFa

=cut
use MooseX::Declare;

class Jobs::Subtasks::BiFa_PSSM_Retrieval extends Jobs::Job {
	use Runtime;

	use Scalar::Util qw(blessed);
	use JSON;

	use Motifs::BiFa_Server_Interface;

	require Serialization::Serializable_Array;
	require Datatypes::Motifs::PSSM;
	require Datatypes::Concepts::Transcription_Factor;
	require Links::Links_Database;
	require Datatypes::Motifs::PSSM_Set;
    use Data::Dumper;

	## used by postprocess to specify the
	## place where the retrieved PSSMs go
	## this parameter is not cached since
	## it is only used in postprocess
	has 'nc__pssmset' => (
		is            => 'rw',
		isa           => 'Str',
		documentation => 'The output PSSM Set Name',
		default       => sub { return "PSSMs"; },
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
		documentation => "1. Non-transfac PSSM sets.|"
		  . "Examples are: PLACE, jaspar, ...",
	);

	## pssm filter. apparently, "." means no filtering
	has 'pssm_filter' => (
		is       => "rw",
		isa      => "Str",
		default  => sub { return "."; },
		required => 1,
		documentation => "2. PSSM Regex filter",
	);

	## selection of species specific PSSMS ( PVI or any combination of these )
	has 'species' => (
		is       => "rw",
		isa      => "Str",
		default  => sub { return "PVINFB"; },
		required => 1,
		documentation => "3. Select the species to use.|"
		  . "Possible values: plants, vertebrates, insects, nematodes, fungi, bacteria, or any combination of these, separated by spaces).",
	);

	## use consensus sequences? We use 1 as default, like the BiFa example script
	has 'use_consensus' => (
		is       => "rw",
		isa      => "Bool",
		default  => sub { return 1; },
		required => 1,
		documentation =>
		  "4. Use consensus sequences as well as score matrices.",
	);

	## dummy to make sure we use the same data version through all bifa jobs and also save it for caching
	## this should be set upon creation, see Jobs::BiFa_Job
	has 'data_and_tf_version' => (
		is      => "rw",
		isa     => "Str",
		documentation => "_hidden.",
		default => sub {
			my $dv = "";
			eval {
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
			$dv = $tv . '.' . $cv;
			};
			if ($@) {
				warn $@;
			}
			return $dv;
		},
	);

=head2 Overloaded validation

=cut

	method validate () {
		$self->SUPER::validate();
	}

=head2 Run BiFa PSSM Retrieval

 Parameters:
 None

 Returns:
 A serializable array of PSSM objects

=cut

	method _run () {
		## get the sequence, store it in the cache.
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
        
        my @pssmNames = $soapIf->PssmSetNames();
        push(@pssmNames, "transfac");

		#	get the matrices out in an array
		my @matrices =  uniq($soapIf->Pssms(
			$self->use_consensus(), 
			$self->species(), 
			$self->pssm_filter(), 
			$dv,
			\@pssmNames,
		));
        
        print "\nNo. matrices: " . scalar(@matrices);

		my @pssms = ();
		
        my $no = 0;
		foreach my $matrix ( @matrices ) {
            $no++;
            print "\nDoing $no";
			my $cache_param = {
				id => $matrix,
				type => "BIFA_MATRIX_RECORD",
			};
						
			my $cache_value = cache->cache_get_hash ( $cache_param );
			my $info;
			my @freqs;
			my @counts;
			
			if (!defined ($cache_value)) {
				$info = $soapIf-> PssmInfo ($matrix, $dv);
				@freqs = $soapIf-> PssmFreqs ($matrix, $dv);
				@counts = $soapIf-> PssmCounts ($matrix, $dv);
				$cache_value = {
					info => $info,
					freqs => \@freqs,
					counts => \@counts,
				};
				cache->cache_put_hash ($cache_param, $cache_value);
			} else {
				$info = $cache_value->{info};
				@freqs = @{$cache_value->{freqs}};
				@counts = @{$cache_value->{counts}};
			}

			my $pssm = Datatypes::Motifs::PSSM->new (
				name => $info->{name} || $matrix,
				accession => $matrix,
				metainfo => $info,
				pseudocount => $info->{pseudo_count} || 0.25,
				pathway => $info->{pathway} || "unknown",
				nsites => $info->{sites} || 0,
			);
			my $pos = 0;
			foreach my $frqs (@freqs) {
				my $char = 0;
				foreach my $freq (@{$frqs}) {
					$pssm->set_probability ($pos, $char, $freq);
					++$char;
				}
				++$pos;
			}
			
			if (!Serialization::Serializable::is_hash($info)) {
				$info = {
					info => $info,
				};
			}
			$info->{counts} = \@counts;

			if (defined ($info->{factors})
			&&  defined ($info->{factors}->{string})) {
				foreach my $fac (@{$info->{factors}->{string}}) {
					my @inf = split /\;/, $fac;
					my $id = trim ($inf[0]);
					
					my $species = $inf[2];
					$species =~ s/^Species\: //;
					$species =~ s/.*,(.+)\./$1/;
					$species = lc (trim ($species));
					$species =~ s/\s/_/g;
					my $tf = Datatypes::Concepts::Transcription_Factor->new (
						id => $id,
						name => trim ($inf[1]),
						species => $species,
					);
					$pssm->add_factor($tf);
				}
			}
			
			## we also put this PSSM into the cache
			my $cacheid = {
				what => "BIFA_PSSM_BY_NAME",
				acc => $pssm->name,
			};
			bless $cacheid, "Serialization::Serializable";
			cache->cache_put($cacheid, $pssm);
			$cacheid = {
				what => "BIFA_PSSM_BY_ACCESSION",
				acc => $pssm->accession,
			};
			bless $cacheid, "Serialization::Serializable";
			cache->cache_put($cacheid, $pssm);
			
			## and make a default PNG for it.
			$pssm->get_png;			

			push @pssms, $pssm;
		}
		print "\nData";
		return Serialization::Serializable_Array->new (@pssms);
	}
	
	
=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		my $return = Datatypes::Motifs::PSSM_Set->new (
			pssms => $result->data
		);
		
		$self->put_dataitem($return, $self->nc__pssmset);
		$self->merge_dataitem_versions($self->nc__pssmset, \&Datatypes::Motifs::PSSM_Set::merged_sets);
		
		return $return;
	}
};
