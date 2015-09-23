
use MooseX::Declare;

=head1 APPLES Datatype base class

We are quite demanding on the data types that are exposed
to the users of APPLES (via web-interface/etc.): They must 
be Serializable Moose classes.

The serialization part is used when storing them in the
Cache and ResultDB database.

The Moose part is used for validation and introspection
when generating HTML forms for object creation. 

This class gives a default implementation of an automatic
form generator to create JSON versions of APPLES datatypes. 

=cut

class Datatypes::Datatype extends Serialization::Serializable {

	use feature qw(state);

	use Runtime;
	use Scalar::Util qw(blessed);
	use Serialization::Serializable;
	use Datatypes::Moose_Types;

=head2 Get this computation's parameters and constraints


 Returns:
 a HASHREF as follows:
 {
 	name => {
 		name => <parameter name>,
 		documentation => <parameter documentation>,
 		validate => <a coderef to a sub which can be called 
 		             to validate this parameter, and returns
 		             a message >
 	}
 }

=cut

	method get_parameters () {
		my %attrs = ();
		for my $attr ( $self->meta->get_all_attributes ) {
			next if $attr->name =~ m/^ns__/;
			
			$attrs{ $attr->name } = {
				name          => $attr->name,
				documentation => $attr->documentation,
				validate      => sub {
					my $value = shift;
					$attr->verify_against_type_constraint($value);
				  }
			};
		}
		return \%attrs;
	}
	
=head2 Get a summary of the values in this datatype 

 Returns: 
 A string with a visual summary of this datatype

=cut
	method get_summary () {
		my $val = blessed ($self) . " : {";
		while ( my ($attr, $entry) = (each %{$self->get_parameters()}) ) {
			$val .= "$attr : $self->{$attr}, ";
		}
		$val.= "}";
		return $val;
	}

=head2 Validate all attributes

This method will validate all attributes in this class by 
setting them to their own value stored in the $self HASHREF. 

=cut

	method validate () {
		for my $attr ( $self->meta->get_all_attributes ) {
			next if $attr->name eq "ns__started";
			$attr->set_value( $self, $self->{ $attr->name } );
		}
	}

}
