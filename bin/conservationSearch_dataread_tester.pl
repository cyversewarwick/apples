#!/usr/bin/perl

#/Projects/cyverse.tools/APPLES/repo/bin$
# conservationSearch_dataread_tester.pl Niben101 Niben101Ctg00074g00004

use strict;
use warnings;

use lib "../webseaweeds/";

# This is the ORANGES Loader for APPLES. Go fruit.
BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
	1;
}



use Data::Dumper;

use Bio::DB::Fasta;

use constant {
	CHROMID		=> 0, # These define the column positions in the lookup hash of arrays
	SEQSTART	=> 1, # The arrangement of the first 6 columns (0-5) matches BED standard
	SEQEND		=> 2,
	GENEID		=> 3,
	SCORE		=> 4, # << We don't use this column in our calculation
	DIRECTION	=> 5,
	FIVESTART	=> 6, # The last four columns do not follow BED format
	FIVEEND		=> 7,
	THREESTART	=> 8,
	THREEEND	=> 9,
	LASTCOLUMN	=> 5,
};


my $lookup_species = $ARGV[0];
my $lookup_accession = $ARGV[1];


my $species_1 = "Niben 101";
my $species_2 = "Arabidopsis TAIR10";

$species_1 =~ s/\s//g; # remove spaces in name
$species_2 =~ s/\s//g;

my @species_list = ($species_1, $species_2);

my $fn_rbh = "../inputs/rbhSearch_result_$species_1\_$species_2.txt";


# Part 1 Sequence Data

my %sequence_info_lookup; # e.g. $sequence_info_lookup{"Niben101_Niben101Scf10386g00010"}[1]
my %all_geneids; # e.g. $all_geneids{"Niben101"} gives you all the gene IDs of Niben101
my %db_fasta;

foreach my $species (@species_list) {

	my $fn_main = "../inputs/$species.bed";
	my $fn_utr5 = "../inputs/$species\_utr5.bed";
	my $fn_utr3 = "../inputs/$species\_utr3.bed";

	my $fn_fasta = "../inputs/$species.fa";

	my $fh;

	open($fh, '<', $fn_main) or die "Could not open file '$fn_main' \n$!";
	while (my $line = <$fh>) {
		chomp $line;
		my @array = split(" ", $line);
		next if $#array < LASTCOLUMN;

		$sequence_info_lookup{$species . "_" . $array[GENEID]} = \@array;

		push @{ $all_geneids{$species} }, $array[GENEID];

	}
	close $fh;

	open($fh, '<', $fn_utr5) or die "Could not open file '$fn_utr5' \n$!";
	while (my $line = <$fh>) {
		chomp $line;
		my @array = split(" ", $line);
		next if $#array < LASTCOLUMN;

		$sequence_info_lookup{$species . "_" . $array[GENEID]}[FIVESTART] = $array[SEQSTART];
		$sequence_info_lookup{$species . "_" . $array[GENEID]}[FIVEEND] = $array[SEQEND];

	}
	close $fh;

	open($fh, '<', $fn_utr3) or die "Could not open file '$fn_utr3' \n$!";
	while (my $line = <$fh>) {
		chomp $line;
		my @array = split(" ", $line);
		next if $#array < LASTCOLUMN;

		$sequence_info_lookup{$species . "_" . $array[GENEID]}[THREESTART] = $array[SEQSTART];
		$sequence_info_lookup{$species . "_" . $array[GENEID]}[THREEEND] = $array[SEQEND];

	}
	close $fh;

	$db_fasta{$species} = Bio::DB::Fasta->new($fn_fasta);
}


if ( grep(/^$lookup_accession$/, @{$all_geneids{$lookup_species}} ) ) {
	print "I have information on '$lookup_accession' of '$lookup_species'.\n";
}


# print $sequence_info_lookup{"Niben101_Niben101Scf10386g00010"}[1] . "\n";

if (exists $sequence_info_lookup{$lookup_species . "_" .$lookup_accession} ) {
	# print "\nThe 5' location of " . $lookup_species . "_" .$lookup_accession . " is " . $sequence_info_lookup{$lookup_species . "_" .$lookup_accession}[1] . "\n";
	print "[\n";
	print join ("\n", @{$sequence_info_lookup{$lookup_species . "_" .$lookup_accession}});
	print "]\n";
} else {
	print "I have no information on '$lookup_accession' of '$lookup_species'.\n";
}

print "The sequences of '$lookup_accession' of '$lookup_species' is [" . $db_fasta{$lookup_species}->seq($lookup_accession) . "]";


# Part 2 RBH

my %rbhs;

open(my $fh, '<', $fn_rbh)
  or die "Could not open file '$fn_rbh'. \n$!";
while (my $line = <$fh>) {
	chomp $line;
	my @array = split(/\t/, $line);
	next if $#array < 3;

	# There are 4 columns in the rbh output file:
	# 	the 1st and last belong to the first species
	#	the 2nd and 3rd belong to the second species
	#	the 1st and 2nd are protein IDs
	#	the 3rd and last are gene IDs if known, "unknown" otherwise
	# If any of the gene IDs are "unknown", we try and extract it from the protein ID
	# by taking a substring before the first ".", if a "." exists in the protein ID
	foreach my $ii (2..3) {
		if ($array[$ii] eq 'unknown') {
			$array[$ii] = index($array[3-$ii], '.') == -1 ? $array[3-$ii] : substr($array[3-$ii], 0, index($array[3-$ii], '.'));
		}
	}
    # $array[2] = $array[2] eq 'unknown' ? substr($array[1], 0, index($array[1], '.')) : $array[2] ;
    # $array[3] = $array[3] eq 'unknown' ? substr($array[0], 0, index($array[0], '.')) : $array[3] ;

    # print "2: $array[2], 3: $array[3] \n";

    $rbhs{$array[3]} = $array[2];
}
close $fh;

print "rbh of '$lookup_accession' is '$rbhs{$lookup_accession}'.\n";
print "rbh of 'Niben101Ctg00074g00004' is '$rbhs{Niben101Ctg00074g00004}'.\n";
print $rbhs{"Niben101Ctg00074g00004"} . "\n" ;
print "rbh of 'Niben101Scf01922g11013' is '$rbhs{Niben101Scf01922g11013}'.\n";

