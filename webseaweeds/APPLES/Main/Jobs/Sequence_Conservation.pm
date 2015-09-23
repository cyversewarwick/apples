#!/usr/bin/perl

use MooseX::Declare;

class Jobs::Sequence_Conservation extends Jobs::Job {

	use Runtime;
	use Datatypes::Moose_Types;

	use Serialization::Serializable_Array;
	use Datatypes::Sequence::Search;
	use Jobs::Subtasks::Seaweed_Job;

	has 'primary_seq' => (
						   is            => "rw",
						   isa           => "Sequence_Search",
						   documentation => '1. Primary sequence',
						   default       => sub { return '' },
	);

	has 'secondary_seq' => (
					   is            => "rw",
					   isa           => "Sequence_Search",
					   documentation => '2. Phylogenetically related sequences',
					   default       => sub { return '' },
	);

	has 'windowsize' => (
						  is            => "rw",
						  isa           => "Int",
						  documentation => 'P1: Window Size',
						  default       => sub { return 60; },
	);

	has 'masked' => (
					  is            => 'rw',
					  isa           => 'Bool',
					  documentation => 'P2: Use repeatmasked sequences',
					  default       => sub { return 0; },
	);

=head2 Overloaded validate method

=cut

	method validate () {
		$self->SUPER::validate();

		if ( $self->windowsize < 10 || $self->windowsize > 500 ) {
			die "Window size must be between 10 and 500";
		}
		my $s                   = Datatypes::Sequence::Search->new;
		my $primary_locations   = $s->parse_search( $self->primary_seq, 0 );
		my $secondary_locations = $s->parse_search( $self->secondary_seq, 0 );

		my $plocs = scalar @{$primary_locations};
		my $slocs = scalar @{$secondary_locations};

		if ( $plocs != 1 ) {
			die
"You should specify exactly one primary sequence location. Your search string gives me $plocs.";
		}
		if ( $slocs < 1 ) {
			die "You should specify at least one secondary sequence location.";
		}
	}

=head2 Overloaded run method

=cut

	method _run () {
		my $s = Datatypes::Sequence::Search->new;

		my $primary_locations   = $s->parse_search( $self->primary_seq,   1 );
		my $secondary_locations = $s->parse_search( $self->secondary_seq, 1 );

		my @arr         = ();
		my $queued_jobs = 0;
		foreach my $sl (@$secondary_locations) {
			eval {
				$main::use_scheduler = 1;
				my $rj = Jobs::Subtasks::Seaweed_Job->new(
					'output_dataset' => $self->output_dataset,
					"sequence_1"     => $primary_locations->[0],
					"sequence_2"     => $sl,
					"windowsize"     => $self->windowsize,
					"masked"         => $self->masked,
				);
				push @arr, $rj->run;
			};
			my $err = $@;
			if ($err) {
				$main::use_scheduler = 0;
				if ( !UNIVERSAL::isa( $err, "Runtime::Job_Queued_Exception" ) )
				{
					confess $err;
				} else {
					$queued_jobs++;
				}
			}
		}

		$main::use_scheduler = 0;
		## don't have all the results? return that we have queued the job
		if ( $queued_jobs > 0 ) {
			die Runtime::Job_Queued_Exception->new( job => $self );
		}
		
		return Serialization::Serializable_Array->new( @arr );
	}

=head2 Overloaded postprocess method

=cut

	method postprocess (Serialization::Serializable $result) {
		foreach my $plot (@{$result->data()}) {
			my $name = $self->primary_seq . "/" . $plot->sequence_2->id();
			if ($self->masked) {
				$name.= " (masked)";
			}
			$name .= ", window:" . $self->windowsize();
			$self->put_dataitem($plot, $name );
		}
		return $result;
	}
}
