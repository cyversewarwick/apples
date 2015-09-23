
use MooseX::Declare;

=head2 Recursive RBH search

=cut

class Jobs::Recursive_RBH extends Jobs::Job {

	with "Jobs::Roles::Dynamic_Data_Loader";
	with "Jobs::Roles::Links_Job";

	use Runtime;

	use Datatypes::Moose_Types;
	use Links::Links_Database;

	use Jobs::Subtasks::RBH;
	require Datatypes::Sequence::Location;
	require Sequences::Database::Relative_Location;

	has "target_species" => (
			is            => "rw",
			isa           => "Str",
			documentation => "1. Target species (separated by semicolon/comma)",
			default       => sub { return "" },
			required      => 1,
	);

	has "max_hits" => (
						is            => "rw",
						isa           => "Int",
						default       => 5,
						documentation => "2. Breadth of search",
	);

	has "startid" => (
		   is      => "rw",
		   isa     => "Str",
		   default => sub { return "ALL CURRENT"; },
		   documentation =>
			 "3. Start at identifiers (optional, separated by semicolon/comma)",
	);

	has "searchtool" => (
						  is            => "rw",
						  isa           => "BLASTTool",
						  default       => default_value("BLASTTool"),
						  documentation => "4. Search tool to use",
	);

=head2 Overridden run method

=cut

	method _run () {
		my @startids = ();

		foreach my $n ( @{ $self->list_selected_nodes } ) {
			if ( UNIVERSAL::isa( $n, "Datatypes::Sequence::Location" ) ) {
				push @startids, $n;
			}
		}

		if ( $self->startid =~ m/[;,]?\s*ALL\s*CURRENT[;,]?/i ) {
			foreach my $n ( @{ $self->list_nodes_by_type ("Datatypes::Sequence::Location") } ) {
				push @startids, $n;
			}
		}


		if ( $self->startid ne "" ) {
			my @ids = split /[,;]/, $self->startid;
			my @sdbs = @{ get_search_database_ids() };
			## for each start id, find all the FASTA databases it matches.
			foreach my $sdb (@sdbs) {
				my $db = get_search_database($sdb);

				foreach my $id (@ids) {
					my $sloc =
					  Sequences::Database::Relative_Location->new(
															  identifier => $id,
															  length     => -1
					  );
					if ( UNIVERSAL::can( $db, 'validate_location' ) ) {
						eval { $db->validate_location($sloc); };
						if ( !$@ ) {
							my $si =
							  Datatypes::Sequence::Location
							  ->new_simple_location( $sdb,
													 $id );
							push @startids, $si;
						}
					}
				}
			}
		}

		if (scalar @startids == 0) {
			die "No start ids selected.";
		}

		my @targetdbs      = split /[,;]/, $self->target_species;
		my @nets           = ();
		my $jobs_submitted = 0;
		my $errors         = "";
		my $queued_jobs    = 0;

		foreach my $tdb (@targetdbs) {

			eval { get_search_database( $tdb . " " . $self->searchtool ); };
			if ($@) {
				warn("Search species $tdb not found.");
				$errors .= "; " unless $errors eq "";
				$errors .= "Search species $tdb not found ($@).";
			}

			foreach my $sid (@startids) {
				debug("Trying target db $tdb for $sid.");

				eval {
#					$main::use_scheduler = 1;
					my $j =
					  Jobs::Subtasks::RBH->new(
								   'output_dataset' => $self->output_dataset,
								   'startID'        => $sid,
								   'targetDB' => $tdb . " " . $self->searchtool,
								   'max_hits' => $self->max_hits,
					  );
					++$jobs_submitted;
					push @nets, $j->run;
				};
				my $err = $@;
				$main::use_scheduler = 0;
				if ($err) {
					if (
						 !UNIVERSAL::isa( $err, "Runtime::Job_Queued_Exception"
						 )
					  )
					{
						confess $err;
					} else {
						$queued_jobs++;
					}
				}
			}
		}

		$main::use_scheduler = 0;

		## don't have all the results? return that we have queued the job
		if ( $queued_jobs > 0 ) {
			die Runtime::Job_Queued_Exception->new( job => $self );
		}

		if ( $jobs_submitted == 0 ) {
			die "No valid new species/id pairs could be added to the graph.";
		}

		my $all_cn = Links::Links_Database->new;
		foreach my $cn (@nets) {
			$all_cn = Links::Links_Database::merge_nodes( $all_cn, $cn );
		}

		return $all_cn;
	}

=head2 This method will be run by the cache after the result has
been retrieved

 Parameters:
 $result : the computation result (cached or just computed)

 Returns:
 a postprocessed version of the result

=cut

	method postprocess ( Serialization::Serializable $result ) {
		$self->merge_links_network($result);
		return $result;
	}

=head2 Expiry time (disable cache)

=cut
	method expiry_time() {
		return 0;
	}

}
