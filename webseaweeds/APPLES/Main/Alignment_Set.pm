### (c) copyright University of Warwick 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Alignment_Set Class ###
use MooseX::Declare;
class Alignment_Set {
  use Data::Dumper;

  has 'matches' => (is => 'rw', isa => 'HashRef', required => 1); # gives array of target-IDs for given query-ID
  
  method give_hits(Str $query_ID) {
    # if query ID is undefined, do what? return empty array?
    my @result = '';  
    my $array_ref = ${$self->matches}{$query_ID};
    if (!defined $array_ref) {
      return @result;
    }
    @result = @{$array_ref};
    return @result; 
  } # give_hits #
  
} # Alignment_Set #

