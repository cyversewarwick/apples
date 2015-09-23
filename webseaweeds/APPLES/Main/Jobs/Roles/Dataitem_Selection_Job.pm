#!/usr/bin/perl

=head1 Jobs that work on a selection of dataitems in the current dataset

This role gives jobs to request a selection of dataitems in the
current dataset as an input.

=cut

package Jobs::Roles::Dataitem_Selection_Job;

use feature ':5.10';

use Runtime;

use Moose::Role;
use MooseX::Declare;

use Datatypes::Moose_Types;
use Digest::SHA;
use JSON;

## jobs can take selected data items as an input
## Each String gives the source id of a data item
## the job will work on the latest version
has 'selected_dataitems' => (
	is        => 'ro',
 	isa     => 'ArrayRef[Str]',
	default => sub { return [] },
	documentation => "_ automatically filled: selected data items",
);

=head2 Return selected objects of a given kind

 Parameters:
 $self : a self object
 $type : the class of the selected object

 Returns:
 an array containing these objects (or array reference in scalar context)

=cut

sub selected_objects_of_kind {
    my $self = shift;
    my $type = shift || "Serialization::Serializable";

    my @arr = ();

    foreach my $di ( @{ $self->selected_dataitems } ) {
        my $x = $self->get_dataitem($di);

        _recursive_find_and_add(\@arr, $x, $type);
    }

    if (wantarray) {
        return @arr;
    }

    return \@arr;
}

=head2 Recursively search an object for subobjects of a given type

 Parameters:
 $arr : an arrayref to add to
 $obj : the object to look at
 $type : what to look for

 Returns:
 Nothing

=cut

sub _recursive_find_and_add {
    my $arr = shift;
    my $obj = shift;
    my $type = shift;

    if ( UNIVERSAL::isa( $obj, $type ) ) {
        push @$arr, $obj;
    } elsif (UNIVERSAL::isa($obj, "Serialization::Serializable_Array") ) {
        foreach my $x (@{$obj->data}) {
            _recursive_find_and_add($arr, $x, $type);
        }
    } elsif (UNIVERSAL::isa($obj, "Datatypes::Sequence::Set") ) {
        foreach my $x (@{$obj->sequences}) {
            _recursive_find_and_add($arr, $x, $type);
        }
    } elsif (UNIVERSAL::isa($obj, "Datatypes::Motifs::PSSM_Set") ) {
        foreach my $x (@{$obj->pssms}) {
            _recursive_find_and_add($arr, $x, $type);
        }
    }
}

1;