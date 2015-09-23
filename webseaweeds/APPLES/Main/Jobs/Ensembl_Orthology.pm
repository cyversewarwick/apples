
use MooseX::Declare;

=head2 Recursive RBH search

=cut

class Jobs::Ensembl_Orthology extends Jobs::Job {

	with "Jobs::Roles::Dynamic_Data_Loader";
	with "Jobs::Roles::Links_Job";

	use Runtime;

	use Datatypes::Moose_Types;
	use Links::Links_Database;

	require Datatypes::Concepts::Gene;
	require Orthology::Database::Ensembl;

	has "target_species" => (
			is            => "rw",
			isa           => "Str",
			documentation => "1. Target species (separated by semicolon/comma)",
			default       => sub { return "" },
			required      => 1,
	);

	has "startid" => (
		   is      => "rw",
		   isa     => "Str",
		   default => sub { return "ALL CURRENT"; },
		   documentation =>
			 "2. Start at identifiers (optional, separated by semicolon/comma)",
	);

	has "item_name" => (
						 is            => "rw",
						 isa           => "Str",
						 default       => sub { return ""; },
						 documentation => "3. Data item name",
	);

=head2 Overridden run method

=cut

	method _run () {
		my @startids = ();

		my $all_cn = $self->get_dataitem( $self->item_name );

		if ( $self->startid =~ m/[;,]?\s*ALL\s*CURRENT[;,]?/i ) {
			foreach my $n ( @{ $all_cn->get_nodes } ) {
				push @startids, $n
				  if UNIVERSAL::isa( $n, "Datatypes::Concepts::Gene" )
				;
			}
		}

		if ( $self->startid ne "" ) {
			my @ids = split /[,;]/, $self->startid;

			foreach my $id (@ids) {
				my $db = "ensembl";
				if ( $id =~ m/(.*)\sin\s+(.*)$/ ) {
					$id = $1;
					$db = $2;
				}
				push @startids, Datatypes::Concepts::Gene->new( $db, $id );
			}
		}

		if ( scalar @startids == 0 ) {
			die "No start ids selected.";
		}

		my @targets = split /[,;]/, $self->target_species;

		my %odbs = ();

		foreach my $sid (@startids) {
			my $ortho_db = $odbs{ $sid->db };

			if ( !defined($ortho_db) ) {
				$odbs{ $sid->db } =
				  Orthology::Database::Ensembl->new( $odbs{ $sid->db } );
				$ortho_db = $odbs{ $sid->db };
			}
			my @ids = @{ $ortho_db->fetch_orthologous_ids( $sid, \@targets ) };

			$sid = $all_cn->get_or_add_node($sid);
			foreach my $i (@ids) {
				$i = $all_cn->get_or_add_node($i); 
				$sid->link_bidirectional($i)->data( 'compara', '1' );
			}
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
		$self->put_dataitem( $result, $self->item_name );
		return $result;
	}

=head2 Expiry time (disable cache)

=cut

	method expiry_time() {
		return 0;
	}

}
