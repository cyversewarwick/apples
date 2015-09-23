### (c) copyright University of Warwick 2013 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Conserved_Region Class ###
#Notes

use MooseX::Declare;

class WPTools::Region_Match
{
    use Data::Dumper;
    has 'ID_match' => (is => 'rw', isa => 'Int');

    has 'start' => (is => 'rw', isa => 'Num');
    has 'end' => (is => 'rw', isa => 'Num');
    
    method new_from_region_match (WPTools::Region_Match $reg) {
        $self->ID_match($reg->ID_match);
        $self->start($reg->start);
        $self->end($reg->end);
    }
	
} # Region_Match #
