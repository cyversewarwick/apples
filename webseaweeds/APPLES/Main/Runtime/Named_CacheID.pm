#!/usr/bin/perl

=head1 Class Runtime::Named_CacheID

This class allows to create cache entries with a named 
identifier.

=cut

use MooseX::Declare;

class Runtime::Named_CacheID extends Serialization::Serializable {
	
	use Runtime;
	
	## We only use an identifier as a string
	## By default, this is undef, in which case
	## we will not store or retrieve anything
	has 'identifier' => (
		is => "rw",
		isa => "Undef|Str",
	);
	
	##ÊGet the result from the cache
	method get_result () {
		if(!defined $self->{identifier}) {
			return undef;
		}
		return cache->cache_get($self);
	}

	##ÊPut the result into the cache
	method put_result (Serialization::Serializable $result) {
		if(!defined $self->{identifier}) {
			return;
		}
		return cache->cache_put($self, $result);
	}
	
	## check if the cache has the result
	method has_result () {
		if(!defined $self->{identifier}) {
			return 0;
		}
		return cache->cache_test($self);
	}
}


1;