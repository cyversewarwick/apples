#!/usr/bin/perl

=head1 This role allows jobs to load data dynamically from the job 
dataset. 

Jobs implementing this role need to do this:

1. Define an attribute 'nc__param' to indicate the data item source this will
   be read from. This attribute should be of type Str.
   
   TODO : an enhancement would be to allow selection of different datasets/
          versions using this attribute.
   
2. Define a non-serialized attribute 'ns__param' which will receive the data 
   item. You can type-check this attribute, but you should not supply another
   default value than undef, since the value will only be loaded if $self 
   doesn't contain it yet.
   
When a cache id is created, nc__param will be replaced by what is returned 
by sub load_data, which is a checksum of the loaded data item.

=cut

package Jobs::Roles::Dynamic_Data_Loader;

use feature ':5.10';

use Runtime;

use Moose::Role;
use MooseX::Declare;

use Datatypes::Moose_Types;
use Digest::SHA qw(sha512_hex);
use JSON;

use Datatypes::Moose_Types;

## we need a method get_dataitem from job to
requires('get_dataitem');

=head2 Load data dynamically

This is called by Cache.pm before using the object as a cache id

This function reads a parameter dynamically from a job's in/output dataset.

Note that this can make computations validate before running, but then fail
when they are executed on a cluster: the dataset can change in between.

You should make sure that jobs don't leave datasets in an unacceptable state.

 Parameters:
 $self : a self object
 $p    : a parameter name (optional, if undef: will load all)
 
 Returns:
 A SHA2 checksum for the data object that was loaded (this is used in 
 caching to detect if we have a cached version of this computation)

=cut

sub load_data {
	my $self = shift;
	my $p    = shift;

	## if no parameter name is given, we load all of them.
	if ( !defined($p) ) {
		my @attrs = $self->meta->get_all_attributes();

		foreach my $attr (@attrs) {
			if ( $attr->name =~ m/^nc__(.*)/ ) {
				$self->load_data( $attr->name );
			}
		}
		return;
	}

	my $p_source = $p;

	if ( $p =~ m/^nc__(.*)/ ) {
		$p = "ns__$1";
	} else {
		confess( "Can only dynamically load non-cached parameters "
			. "(parameter name must start with 'nc__' -- failed for parameter $p)"
		);
	}

	## prevent multiple loading
	return if ( defined( $self->{$p} ) );

	my $load_or_create = "";

	my $type = undef;

	my @attrs = $self->meta->get_all_attributes();
	foreach my $attr (@attrs) {
		if (    $attr->name eq $p ) {
			my $t = undef;
			
			if ( $attr->has_type_constraint ) {
				$t = $attr->type_constraint->name;
			}
			
			if ( defined ($t) && type_isa_class ( $t, 'UNIVERSAL') ) 
			{			
				$type = $t;
			}
		}
	}

	my $item = undef;
	$load_or_create .=
	  '$item = $self->get_dataitem ( $self->' . $p_source . ');';
	eval $load_or_create;
	if ($@) {
		warn ("Loading of item $p failed: " . $@);
	}
	if ( !defined ($item) ) {
		if ( defined ($type) && $type ne "" ) {
			$load_or_create =
			    'if (!defined ($item)) { $item = ' 
			  . $type
			  . '->new(); }';
			eval $load_or_create;
			if ($@) {
				die $@;
			}		
		}
	}
	
	if (!defined ($item)) {
		die "Data item $p could neither be loaded nor created.";
	}
	
	$load_or_create = '$self->' . $p . '( $item );';
	eval $load_or_create;
	if ($@) {
		die $@;
	}		
	

	if ( !defined( $self->{$p} ) ) {
		confess(
"Failed to load item $p dynamically : data for source $self->{$p_source} was not valid."
		);
	}

	$self->{ "ns__checksum__" . $p_source } =
	  sha512_hex(
				  to_json(
						   Serialization::Serializable::to_hash( $self->{$p} ),
						   { allow_blessed => 1 }
				  )
	  );
}

1;
