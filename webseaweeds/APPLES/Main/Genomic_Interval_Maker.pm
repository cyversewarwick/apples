### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Genomic_Interval_Maker class ###
use MooseX::Declare;
use Genome_DB_Utilities;

class Genomic_Interval_Maker {
	use Genomic_Interval;
	
	method copy_genomic_interval (Genomic_Interval $gi) {
		my $genomic_interval = Genomic_Interval->new(
			'genome_db' => $gi->genome_db,
			'region' => $gi->region,
			'five_prime_pos' => $gi->five_prime_pos,
			'three_prime_pos' => $gi->three_prime_pos,
			'strand' => $gi->strand,
			'working_sequence' => $gi->working_sequence,
		);
		if (defined $gi->coord_sys_name) {
			$genomic_interval->coord_sys_name($gi->coord_sys_name);
		}
		if (defined $gi->type) {
			$genomic_interval->type($gi->type);
		}
		if (defined $gi->label) {
			$genomic_interval->label($gi->label);
		}
		if (defined $gi->sub_interval_annotation) {
			$genomic_interval->sub_interval_annotation($gi->sub_interval_annotation);
		}
		return $genomic_interval;
	} # copy_genomic_interval #
	
} # Genomic_Interval_Maker #