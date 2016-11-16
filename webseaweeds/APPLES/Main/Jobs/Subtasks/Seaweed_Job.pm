#!/usr/bin/perl



use MooseX::Declare;

=head2 Two sequence input seaweed job.

=cut

class Jobs::Subtasks::Seaweed_Job extends Jobs::Job {

	use Configuration::AppleSeeds;

	use Runtime;

	use Sequences::Genomic_Sequence;

	use File::Temp qw/ tempfile /;
	use Cwd;
	use List::Util qw(min max);

	use POSIX;# ":sys_wait_h";
	use Time::HiRes qw(sleep);
	use Test::More 'no_plan';
	use DBI;

	use Data::Dumper;

	require Datatypes::Results::Alignment_Plot;

	has 'sequence_1' => (
						  is       => "rw",
						  isa      => "Datatypes::Sequence::Location",
						  required => 1,
	);

	has 'sequence_2' => (
						  is       => "rw",
						  isa      => "Datatypes::Sequence::Location",
						  required => 1,
	);

	has 'windowsize' => (
						  is            => "rw",
						  isa           => "Int",
						  documentation => 'P1: Window Size',
						  default       => sub { return 60; },
	);

	has 'dbhs' => (
					is 				=> 'rw',
					traits  		=> ['Array'],
					isa				=> 'ArrayRef[DBI::db]',
					documentation	=> 'Array of DBI handles',
					required 		=> 0,
					default			=> sub { [] },
					handles 		=> {
							add_handle=>'push',
						},
	);

	has 'masked' => (
					  is            => 'rw',
					  isa           => 'Bool',
					  documentation => 'P0: Use repeatmasked sequences',
					  default       => sub { return 1; },
	);

	has 'step_1' => (
					  is            => "rw",
					  isa           => "Int",
					  documentation => 'P2: Step size in sequence 1',
					  default       => sub { return 1; },
	);

	has 'step_2' => (
					  is            => "rw",
					  isa           => "Int",
					  documentation => 'P2: Step size in sequence 2',
					  default       => sub { return 1; },
	);

=head2 Overloaded validate method

=cut

	method validate () {
		$self->SUPER::validate();

		if ( $self->windowsize < 10 || $self->windowsize > 500 ) {
			die "Window size must be between 10 and 500";
		}
	}

=head2 Overloaded run method

=cut

	method _run () {
		my $seaweed_executable = find_executable("windowalignment");

		if ( !-e $seaweed_executable ) {
			die "Seaweed Code executable $seaweed_executable not found.";
		}

		my $primary   = $self->sequence_1->get_sequence()->[0];
		my $secondary = $self->sequence_2->get_sequence()->[0];

		## set seaweed algorithm parameters
		my $primary_stepsize   = $self->step_1();
		my $secondary_stepsize = $self->step_2();
		my $windowlength       = $self->windowsize();

		my $primary_size   = length( $primary->seq );
		my $secondary_size = length( $secondary->seq );

		## check size
		my $size_limit_1 = get_config_key("windowalignment_size_limit_1") || 20000;
		my $size_limit_2 = get_config_key("windowalignment_size_limit_2") || 20000;

		if ( $primary_size * $secondary_size > $size_limit_1 * $size_limit_2 ) {
			die(   "FAILED: This job ($primary_size x $secondary_size = "
				 . $primary_size * $secondary_size
				 . ") is too big! The limit is $size_limit_1 * $size_limit_2."
			);
		}
        
        my $directory = getcwd . "/tempfiles/";

		my ( $primary_sequence_file, $primary_sequence_file_name ) =
		  tempfile( DIR => $directory, UNLINK => 0 );
		my ( $primary_reversed_sequence_file,
			 $primary_reversed_sequence_file_name )
		  = tempfile( DIR => $directory, UNLINK => 0 );
		my ( $secondary_sequence_file, $secondary_sequence_file_name ) =
		  tempfile( DIR => $directory, UNLINK => 0 );
		my ( $secondary_reversed_sequence_file,
			 $secondary_reversed_sequence_file_name )
		  = tempfile( DIR => $directory, UNLINK => 0 );

		my $prseq = $primary->seq;
		my $scseq = $secondary->seq;

		if ( $self->masked ) {
			$prseq = $primary->{masked_sequence}   || $primary->seq;
			$scseq = $secondary->{masked_sequence} || $secondary->seq;
		}

		my $rprseq =
		  Sequences::Genomic_Sequence::reverse_dna_sequence_strand($prseq);
		my $rscseq =
		  Sequences::Genomic_Sequence::reverse_dna_sequence_strand($scseq);

		if ( length($prseq) == 0 ) {
			die "Primary sequence has length zero.";
		}

		if ( length($scseq) == 0 ) {
			die "Secondary sequence has length zero.";
		}

		if ( length($rprseq) == 0 ) {
			die "Reversed primary sequence has length zero.";
		}

		if ( length($rscseq) == 0 ) {
			die "Reversed secondary sequence has length zero.";
		}

		print $primary_sequence_file $prseq;
		print $secondary_sequence_file $scseq;
		print $primary_reversed_sequence_file $rprseq;
		print $secondary_reversed_sequence_file $rscseq;
        
        #print "\nPrim: " . $primary_sequence_file_name;
        #print "\nPrim: " . $secondary_sequence_file_name;
        #print "\nPrim: " . $primary_reversed_sequence_file_name;
        #print "\nPrim: " . $secondary_reversed_sequence_file_name;

		close $primary_sequence_file;
		close $primary_reversed_sequence_file;
		close $secondary_sequence_file;
		close $secondary_reversed_sequence_file;

		# my $code_parameters =
		#     "-w $windowlength "
		#   . "-s $primary_stepsize "
		#   . "-t $secondary_stepsize " . "-V";

		my $code_parameters =
		    "-w $windowlength "
		  . "-m seaweednw ";

		## compare both forward and reverse strands.
		## the only direction we don't need to compare are
		## the two reverse strands, this is equivalent
		## to comparing the forward strands

		my $prefix = $primary_reversed_sequence_file_name;

        #print "\nPref: ${prefix}\n";
		# construct the call to the seaweed code
		# old-style calling
		# my @plot_exec = (
		# 				  "$seaweed_executable "
		# 					. "$primary_sequence_file_name "
		# 					. "$secondary_sequence_file_name "
		# 					. "${prefix}_result_nn.dat "
		# 					. "${prefix}_pp_nn.dat "
		# 					. "${prefix}_sp_nn.dat "
		# 					. $code_parameters,
		# 				  "$seaweed_executable "
		# 					. "$primary_sequence_file_name "
		# 					. "$secondary_reversed_sequence_file_name "
		# 					. "${prefix}_result_nr.dat "
		# 					. "${prefix}_pp_nr.dat "
		# 					. "${prefix}_sp_nr.dat "
		# 					. $code_parameters,
		# );

		# new-style calling
		my @plot_exec = ( "$seaweed_executable " 
							. "--first-sequence $primary_sequence_file_name "
							. "--second-sequence $secondary_sequence_file_name "
							. "--output-file ${prefix}_nn "
							. $code_parameters ,
						  "$seaweed_executable "
							. "--first-sequence $primary_sequence_file_name "
							. "--second-sequence $secondary_reversed_sequence_file_name "
							. "--output-file ${prefix}_nr "
							. $code_parameters	);

		# foreach my $pe (@plot_exec) {
		# 	# here stuff gets run!
		# 	debug("Running command: $pe\n");
		# 	print "Running command: $pe\n";

			# qx($pe 2>&1);
		# }

		my @file_suffix = ("nn","nr");
		foreach my $ii (0..1) { # iterate commands
			# here stuff gets run!

			foreach my $jj (1..3){ # iterate attempts

				eval{
					local $SIG{ALRM} = sub{ die "alarm_sound\n"};
					alarm(10*$jj); #seconds
					# eval{
						print("Running command (attempt no. $jj): $plot_exec[$ii]\n");
						qx($plot_exec[$ii] 2>&1);
					# };
					alarm(0);
				};
				alarm(0); # eliminate race condition, http://docstore.mik.ua/orelly/perl4/cook/ch16_22.htm

				if ($@) {
					die "Unexpected exit from seaweed" unless $@ eq "alarm_sound\n";
					print "Timeout (10*$jj secs) called on seaweed\n";
				}

				if (-f "${prefix}_$file_suffix[$ii]\_profile_1" &&
					-f "${prefix}_$file_suffix[$ii]\_profile_2" &&
					-f "${prefix}_$file_suffix[$ii]\_result") {
					print "'$file_suffix[$ii]' output files are ready.\n";
					last; # iterate attempts
				} else {
					print "'$file_suffix[$ii]' output files are not found.\n";
				}

			}

		}

		# my @file_suffix = ( "nn" , "nr"); # Todo in the next commit

		# my $timeout = 10; # seconds, max to wait for seaweed binary to complete
		# # my @db_handles = $self->dbhs();

		# for (my $i=1; $i<=3; $i++) {

		# 	print "Seaweed_Job: Running command 1/2, ($i/3): \n$plot_exec[0]\n";

		# 	my $start_time = time();
		# 	my $child_pid = fork();
		# 	die "Could not fork.\n" if not defined $child_pid;

		# 	if (not $child_pid) {
		# 		# In child
		# 		# print "Seaweed_Job: db_handles: ". Dumper(@db_handles)."\n";
		# 		# print Dumper($db_handles[1]);
		# 		# isa_ok($db_handles[1],'DBI::db');
		# 		# $db_handles[1]->{"InactiveDestroy"} = 1;
		# 		# $db_handles[2]->{"InactiveDestroy"} = 1;
		# 		foreach my $dbh (@{ $self->dbhs }) {
		# 			print "Seaweed_Job: dbh: " . Dumper($dbh) ."\n";
		# 			isa_ok($dbh, 'DBI::db');
		# 			$dbh->{"InactiveDestroy"} = 1;
		# 		}

		# 	print qx($plot_exec[0] 2>&1);
		# 	print "Seaweed_Job: Child process completed.\n";
		# 	exit 3; # Child completed.

		# 	} else {
		# 		# In parent
		# 		while (1) {

		# 			my $res = waitpid($child_pid, WNOHANG);
		# 			sleep(0.1);

		# 			if ($res == -1) {
		# 				print "Error from seaweed binary: ". ($? >> 8) ."\n";
		# 				last;
		# 			}

		# 			if ($res) {
		# 				last;
		# 			}

		# 			if (time() - $start_time > $timeout) {
		# 				print "Timeout, killing child pid: $child_pid.\n";
		# 				kill (SIGKILL, $child_pid); # SIGKILL = 15
		# 				last;
		# 			}
		# 		}
		# 	}

		# 	sleep $i-1;

		# 	if (-f "${prefix}_nn_profile_1" &&
		# 			-f "${prefix}_nn_profile_2" &&
		# 			-f "${prefix}_nn_result") {
		# 		print "nn output files are ready.\n";
		# 		last;
		# 	} else {
		# 		print "No nn output files, run command again.\n";
		# 	}
		# }

		# for (my $i=1; $i<=3; $i++) {

		# 	print "Seaweed_Job: Running command 2/2, ($i/3): \n$plot_exec[1]\n";
			
		# 	my $start_time = time();
		# 	my $child_pid = fork();
		# 	die "Could not fork.\n" if not defined $child_pid;

		# 	if (not $child_pid) {
		# 		# In child
		# 		foreach my $dbh (@{ $self->dbhs }) {
		# 			print "Seaweed_Job: dbh: " . Dumper($dbh) ."\n";
		# 			isa_ok($dbh, 'DBI::db');
		# 			$dbh->{"InactiveDestroy"} = 1;
		# 		}
				
		# 	print qx($plot_exec[1] 2>&1);
		# 	exit 3; # Child completed.

		# 	} else {
		# 		# In parent
		# 		while (1) {

		# 			my $res = waitpid($child_pid, WNOHANG);
		# 			sleep (0.1);

		# 			if ($res == -1) {
		# 				print "Error from seaweed binary: ". ($? >> 8) ."\n";
		# 				last;
		# 			}
		# 			if ($res) {
		# 				last;
		# 			}
		# 			if (time() - $start_time > $timeout) {
		# 				print "Timeout, killing child pid: $child_pid.\n";
		# 				kill (SIGKILL, $child_pid); # SIGKILL = 15
		# 				last;
		# 			}
		# 		}
		# 	}

		# 	sleep $i-1;

		# 	if (-f "${prefix}_nr_profile_1" &&
		# 			-f "${prefix}_nr_profile_2" &&
		# 			-f "${prefix}_nr_result") {
		# 		print "nr output files are ready.\n";
		# 		last;
		# 	} else {
		# 		print "No nr output files, run command again.\n";
		# 	}
		# }

        #print "\nReading PROFILE from : ${prefix}";
		# read profiles and plots
		my $profile1 = $self->_read_profile("${prefix}_nn_profile_1");
        #print "\nRead one";
		my $profile2 = $self->_read_profile("${prefix}_nn_profile_2");
		my $plot     = $self->_read_plot("${prefix}_nn_result");

		my $profile1_a = $self->_read_profile("${prefix}_nr_profile_1");
		my $profile2_a = $self->_read_profile("${prefix}_nr_profile_2");
		my $plot_a     = $self->_read_plot("${prefix}_nr_result");

		my $xlen = length($prseq) - $windowlength + 1;
		my $ylen = length($scseq) - $windowlength + 1;
		foreach my $pos ( keys %$plot_a ) {
			if ( $pos =~ m/([0-9]+)\_([0-9]+)/ ) {
				my $xp  = $1;
				my $yp  = $2;
				my $ypr = $ylen - $yp - 1;
				my $val = $plot->{"${xp}_${ypr}"} || 0;
				$plot->{"${xp}_${ypr}"} = max( $val, $plot_a->{$pos} );
			} else {
				warn("invalid plot key.");
			}
		}

		## read and merge histograms
		my @histograms = (
						   $self->_read_profile(
										   "${prefix}_nn_histogram"),
						   $self->_read_profile(
										   "${prefix}_nr_histogram"),
		);

		## histograms should have the same length
		my $i    = 0;
		# my @new_hist = (0) x ($windowlength+1);

		my @new_hist = ();
		for($i = 0; $i < $windowlength + 1; ++$i) {
			$new_hist[$i] = 0;
		}

		my $ii = 0;
		foreach my $h (@histograms) {
			for($i = 0; $i < (scalar @$h); ++$i) {
				$new_hist[($windowlength+1)*$i/(scalar @$h)] += $h->[$i];
				# print "hist $ii : $i is mapped to ", int(($windowlength+1)*$i/(scalar @$h)), "\n";
			}
			++$ii;
		}
		# dump_json(\@new_hist);

		for ( $i = 0 ; $i < scalar @{$profile1} ; ++$i ) {
			$profile1->[$i] = max( $profile1->[$i], $profile1_a->[$i] );
		}

		my $plen = scalar @{$profile2};
		if (scalar @{$profile2_a} != $plen) {
			die "Output profiles don't have the same length! : $plen != ". scalar @{$profile2_a};
		}
		for ( $i = 0 ; $i < $plen; ++$i ) {
			$profile2->[$i] = max( $profile2->[$i], $profile2_a->[$plen - $i - 1] );
		}
        
        #Remove all the temp files
        #system("rm ${prefix}* 2>/dev/null");
        #system("rm $primary_reversed_sequence_file_name 2>/dev/null");
        #system("rm $secondary_reversed_sequence_file_name 2>/dev/null");
        #system("rm $primary_sequence_file_name 2>/dev/null");
        #system("rm $secondary_sequence_file_name 2>/dev/null");
        #Remove all the temp files
        
		return
		  Datatypes::Results::Alignment_Plot->new( job        => $self,
												   histogram  => \@new_hist,
												   profile_1  => $profile1,
												   profile_2  => $profile2,
												   plot       => $plot,
												   sequence_1 => $primary,
												   sequence_2 => $secondary,
		  );

	}


=head2 Read a conservation profile

 Parameters:
 $filename : the file name

 Returns:
 an ArrayRef[Num] with all profile values
=cut

	method _read_profile(Str $filename) {
		my @profile = ();

		open P1, "<$filename" or do {
			die("Profile computation: profile $filename not found.");
		};
		while (<P1>) {
			my $v = $_;
			if ($v =~ m/^\s*\#.*/) {
				next;
			}

			$v =~ s/[\s\n]//g;

			## some kind of numbers should be in there.
			if (m/[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?/) {
				push @profile, $v;
			}
		}
		close P1;

		return \@profile;
	}

=head2 Read a conservation plot

 Parameters:
 $filename : the file name

 Returns:
 ( ArrayRef[[x,y,score]] with all values > $threshold, $threshold )
=cut

	method _read_plot(Str $filename) {
		my %plot = ();
		open P1, "<$filename" or do {
			warn("Alignment plot computation: $filename not found.");
		};

		my $threshold = undef;
		while (<P1>) {
			my $line = $_;
			if ( !defined($threshold) ) {
				$threshold = int($line);
				if ( !defined($threshold) ) {
					die "Output file $filename has invalid file format.";
				}
			} else {
				if ( $line =~ m/^([0-9]+)\s+([0-9]+)\s+([0-9]+.*)$/ ) {
					my $xp = $1;
					my $yp = $2;
					$plot{"${xp}_${yp}"} = $3;
				} else {
					warn("Output line ignored: $line");
				}
			}
		}
		close P1;
		return \%plot;
	}
}
