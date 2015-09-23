### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Memory_Cache Class ###
# used by the Job_Handler to store some of the results obtained from cache in memory to
# avoid repetitive retrieval

use MooseX::Declare;

class Memory_Cache {

  has 'memory_cache_data' => (is => 'rw', isa => 'HashRef');

  my $GU = General_Utilities->new();

  method retrieve (Str $keyfileMD5, Str $method) {
      $self->private_ensure_hash_is_defined();
      my $result;
      my $key = $keyfileMD5.$method;
      $result = ${$self->memory_cache_data}{$key};
      return $result;
  } # retrieve #

  method store (Str $keyfileMD5, Str $method, Any $result) {
      $self->private_ensure_hash_is_defined();
      my $key = $keyfileMD5.$method;
      ${$self->memory_cache_data}{$key} = $result;
  } # store #

  method private_ensure_hash_is_defined () {
      if (!defined $self->memory_cache_data) {
	  my %hash;
	  $self->memory_cache_data(\%hash);
      }
  } # private_ensure_hash_is_defined #

} # Memory_Cache #
