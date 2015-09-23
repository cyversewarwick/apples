
=head1 Sequence location input parser

Parses strings into sequence locations

=cut

package Datatypes::Sequence::Search;

use strict;

use Parse::RecDescent;
use Data::Dumper;

use Configuration::AppleSeeds;
use Runtime;

use Datatypes::Sequence::Location;
use Datatypes::Roles::Locatable;
use Sequences::Database::Relative_Location;

=head2 Constructor

 Parameters:
 $class : the class 
 
 Returns: 
 a new parser object

=cut

sub new {
	my ($class) = @_;
	my $parser = new Parse::RecDescent(
		q{
	startrule: len location(s?) geneids {
		$return = {
			"length" => $item[1],
			location => $item[2],
			genes => $item[3],
		};
	}
		
	geneids: /.+/ { 
		$return = $item[1];
	}
	
	len: position {
		$return = $item[1];
	} | "" {
		$return = 1000;
	}
	
	location: loc_end position loc_direction /of/i loc_anchor(?) {
		if (defined ($item[5][0]) ) {			 
			$return = {
				$item[1] => {
					offset => $item[2] * $item[3],
					anchor => $item[5][0],
				},
			};
		} else {
			if ($item[3] < 0) {
				$return = {
					$item[1] => {
						offset => $item[2] * $item[3],
						anchor => "5' end",
					},
				};
			} else {
				$return = {
					$item[1] => {
						offset => $item[2] * $item[3],
						anchor => "3' end",
					},
				};
			}
		}
	} | loc_direction /of/i loc_anchor {
		if ($item[1] < 0) {
			$return = {
				'to' => {
					offset => 0,
					anchor => $item[3],
				},
			};
		} else {
			$return = {
				'from' => {
					offset => 0,
					anchor => $item[3],
				},
			};
		}
	} | loc_direction /of/i {
		if ($item[1] < 0) { ##�upstream of -> meaning upstream of identifier
			$return = {
				'to' => {
					offset => 0,
					anchor =>  "5' end",
				},
			};
		} else {  ## downstream of -> downstream of idenitifier
			$return = {
				'from' => {
					offset => 0,
					anchor =>  "3' end",
				},
			};
		}
	} | loc_end loc_anchor {
		$return = {
			$item[1] => {
				offset => 0,
				anchor => $item[2],
			},
		};		
	} | /((ex)|(in))c(lud(e|(ing)))?/i /neighbou?rs/i {
		$return = {
			stop_at_neighbours => (substr ($item[1], 0, 1) eq "e" ? 1 : 0),
		};		
	} | /in/i quoted_word {
		$return = {
			db => $item[2],
		}
	}
	
	loc_end : (/from/i|/to/i) {
		$return = lc ($item[1]);
	}
	
	loc_direction: (/upstream/i|/downstream/i) {
		$return = ({
			upstream => -1,
			downstream => 1,
		})->{lc ($item[1])};
	} 
	
	loc_anchor: ("5'"|"3'"|/start/i) (/end/i|"") {
		$return = ({
			"" => "5' end",
			"5'" => "5' end",
			"start" => "5' end",
			"3'" => "3' end",
		})->{$item[1]};
	}
	
	quoted_word: word {
		$return = $item[1];
	} | "'" dqword(s?) "'" {
		$return= (join " ", @{$item[2]});
	}  | '"' sqword(s?) '"' {
		$return= (join " ", @{$item[2]});
	}
	
	word : /[^\s'"]+/ {
		$return = $item[1];
	}

	dqword : /[^\s']+/ {
		$return = $item[1];
	}
	
	sqword : /[^\s"]+/ {
		$return = $item[1];
	}

	position: /([\+\-]?[0-9]+)[ \t]*([mMkK]?)[ \t]*([bB]?[Pp]?|bases)/ {
		if ($item[1] =~ m/([\+\-]?[0-9]+)[ \t]*([mMkK]?)[ \t]*([bB]?[Pp]?|bases)/) {
			$return = int($1) * (lc($2) eq 'm' ? 1e6 : (lc($2) eq 'k' ? 1000 : 1) ); 
		} else {
			$return = int($item[1]);
		}
	}
	}
	);

	my $self = { 'seq_parser' => $parser, };

	bless $self, $class;
	return $self;
}

=head2 Parse a search string

 Parameters:
 $self : a Datatypes::Sequence::Search object
 $search : the search string
 $lookup : (optional) try to determine the search database

 Returns:
 an ArrayRef [ Datatypes::Sequence::Location ]

=cut

sub parse_search {
	my ( $self, $search, $lookup ) = @_;

	if ( !defined($lookup) ) {
		$lookup = 1;
	}
	$lookup = 0;

	my $s = $self->{'seq_parser'}->startrule($search);

	my @gids =
	  grep { $_ ne "" }
	  ( grep ( !/^(and)|(of)|(for)$/i, ( split /[\s\,\;\:]/, $s->{genes} ) ) );

	my $loc = {};

	if ( defined( $s->{location} ) ) {
		foreach my $l ( @{ $s->{location} } ) {
			foreach my $k ( keys %$l ) {
				$loc->{$k} = $l->{$k};
			}
		}
	}

	$s->{location} = $loc;

	my @locs = ();
	foreach my $gid (@gids) {
		my $dbi = $loc->{db};
		$dbi = Datatypes::Roles::Locatable::_find_db_identifier( $loc->{db},
																 $lookup );

		if ( !defined($dbi) ) {
			if ( !$lookup ) {
				$dbi =
				  Datatypes::Roles::Locatable::_guess_db_by_identifier($gid);
			} else {
				$dbi = Datatypes::Roles::Locatable::_find_db_identifier(
					Sequences::Database::Relative_Location->new(
															  identifier => $gid
					  )
				);
			}

		}

		if ( !$dbi ) {
			my $message =
"Could not determine sequence database for identifiers @gids (lookup: $lookup)";
			if ( defined( $loc->{db} ) ) {
				$message .= " in $loc->{db}";
			}
			die $message;
		}

		if (    defined( $loc->{from} )
			 && defined( $loc->{to} ) )
		{
			die "You cannot specify both \"from [...]\" and \"to [...]\" ";
		}

		my $offset = -$s->{length};
		my $anchor = "5' end";

		if ( defined( $loc->{from} ) ) {
			$offset = $loc->{from}->{offset};
			$anchor = $loc->{from}->{anchor};
		} elsif ( defined( $loc->{to} ) ) {
			$offset = $loc->{to}->{offset} - $s->{length};
			$anchor = $loc->{to}->{anchor};
		} else {
			## use defaults : length bo from 5' end
		}

		push @locs,
		  Datatypes::Sequence::Location->new(
					  db => $dbi,
					  location =>
						Sequences::Database::Relative_Location->new(
						  identifier         => $gid,
						  length             => $s->{length} || 1000,
						  stop_at_neighbours => $loc->{stop_at_neighbours} || 0,
						  offset             => $offset,
						  anchor             => $anchor
						)
		  );

	}

	return \@locs;
}

=head2 Run a search as a test

 Parameters:
 $self : a self object
 $search : a search string
 
 Returns: 
 Nothing (everything is dumped to stdout)

=cut

sub test_search {
	$::RD_HINT  = 1;
	$::RD_TRACE = 1;

	my $self   = shift;
	my $search = shift;

	my $parser = $self->{'seq_parser'};

	print Dumper( $parser->startrule($search) );

	$::RD_HINT  = 0;
	$::RD_TRACE = 0;
}

1;
