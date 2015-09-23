=pod

=head1 Class Serializable

Serializable - Class for outputting various data structures to XML for re-use in webservices and other places.

=head1 SYNOPSIS
If a class wishes to be able to export itself and all the data it holds to an XML file, then it needs to extend Serializable.

=head1 DESCRIPTION

=head2 Methods

=over 12

=item C<to_hash>

Takes all class attributes and puts them into hash ready to be written to XML file. N.B. all members starting 
with "ns__" will be skipped and not exported to XML file. 

=item C<from_hash>

NOT IMPLEMENTED. 
Reverse of to_hash(), reads data from the hash and constructs the class. Maybe should be implemented by each class?

=back

=head1 LICENSE

This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package and distributed 
under Academic Non-Commercial Use Licence.

=head1 COPYRIGHT

(c) Copyright University of Warwick 2009-2010

=head1 AUTHOR

=cut

#BEGIN {
#  our $global_print_level = "1";
#  use lib '../';
#}

use MooseX::Declare;

class Serializable {
    use Data::Dumper;

=pod
    
=head2 from_hash()

NOT IMPLEMENTED. 
Reverse of to_hash(), reads data from the hash and constructs the class. Maybe should be implemented by each class?

=cut

    method from_hash() {
        print "NOT IMPLEMENTED!\n";
    }

=pod
    
=head2 to_hash()

Takes all class attributes and puts them into hash ready to be written to XML file. N.B. all members starting 
with "ns__" will be skipped and not exported to XML file.

Returns a hash ready to be written by XML_Writer.

=cut
    
    method to_hash() {
        my $self_obj = shift;
        my $object = {};

    
        while ( my ($key, $value) = each %$self_obj ) {
            if( $key !~ m/^ns__/) {    # filter bits which should not be serialized
                
                if(ref($value) eq 'ARRAY') { #If element is an array then go through elements
                    my @elemArray;
                    foreach my $element (@{$value}) {
                        if($element->isa('Serializable')) {
                            push(@elemArray, $element->to_hash()); #if element is serializable then hash it
                        } else {
                            push(@elemArray, $element); #otherwise store it's value directly
                        }
                    }
                    
                    $object->{$key} = [@elemArray];
                    
                } elsif (ref($value) ne '') { #if it is not an array but also not a class
                    if($value->isa('Serializable')) {
                        $object->{$key} = $value->to_hash; #if it is serializable then hash it
                    } else {
                    	$object->{$key} = $value; #otherwise store it's value
                    }
                } else {
                	$object->{$key} = $value; #if it just plain value then store it directly
                }
                
            }
        }
        
        return $object;     
    }

}    # Serializable
