### (c) copyright University of Warwick 2013 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Conserved_Region Class ###
#Notes

use MooseX::Declare;

class WPTools::Conserved_Region
{
    use Data::Dumper;
    use WPTools::Region_Match;

    has 'ID' => (is => 'rw', isa => 'Int');

    has 'seq' => (is => 'rw', isa => 'Str');
    has 'species' => (is => 'rw', isa => 'Str');
    has 'gene_accession' => (is => 'rw', isa => 'Str');
    has 'five_prime_pos' => (is => 'rw', isa => 'Int');

    has 'start' => (is => 'rw', isa => 'Num');
    has 'end' => (is => 'rw', isa => 'Num');
    has 'other_species_matches' => (is => 'rw', isa => 'ArrayRef');
    
    method new_from_conserved_region (WPTools::Conserved_Region $reg) {
        $self->ID($reg->ID);
        $self->seq($reg->seq);
        $self->species($reg->species);
        $self->gene_accession($reg->gene_accession);
        $self->five_prime_pos($reg->five_prime_pos);
        $self->start($reg->start);
        $self->end($reg->end);
        
        foreach my $other (@{$reg->other_species_matches})
        {
            my $new_match = WPTools::Region_Match->new;
            $new_match->new_from_region_match($other);
            push(@{$self->{"other_species_matches"}}, $new_match);
        }
    }
	
} # Conserved_Region #
