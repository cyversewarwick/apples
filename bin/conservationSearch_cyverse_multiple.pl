#!/usr/bin/perl

# derived form Nathaniel's conservation_multiple.pl
# my @species = ("nasonia vitripennis", "apis mellifera", "atta cephalotes", "bombyx mori", "acyrthosiphon pisum", "danaus plexippus", "heliconius melpomene", "tribolium castaneum");
# my @dbs = ("nasonia_vitripennis_core_18_71_1", "apis_mellifera_core_18_71_1", "atta_cephalotes_core_18_71_1", "bombyx_mori_core_18_71_1", "acyrthosiphon_pisum_core_18_71_1", "danaus_plexippus_core_18_71_1", "heliconius_melpomene_core_18_71_1", "tribolium_castaneum_core_18_71_3");

use strict;

# BASEDIR contains path to ORANGES base directory
use lib "../webseaweeds/";

# This is the ORANGES Loader for APPLES. Go fruit.
BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
	1;
}

#<=== Includes ===>#

#General
use JSON;
use NDU::Gen qw(:ALL);
use Data::Dumper;

use File::Spec;
use Try::Tiny;

#Read Commandline Options
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

#Grannysmith
use Runtime;
# use Sequences::Database::Relative_Location;
# use Sequences::Database::Sequence::Ensembl;
# use Sequences::Database::Sequence::Genbank;
use Datatypes::Sequence::Local;
use Jobs::Subtasks::Seaweed_Job;
use Serialization::Serializable;
use List::Util qw(shuffle);

#APPLES
use Star_Bundler;
use Genomic_Interval;
use Genomic_Interval_Set;
use Genomic_Interval_Set_Maker;
use Evolutionary_Tree;
use Partial_Threshold_Matrix_Maker;
use Partial_Threshold_Matrix;
use Parameters;

use constant {
    CHROMID     => 0, # These define the column positions in the lookup hash of arrays
    SEQSTART    => 1, # The arrangement of the first 6 columns (0-5) matches BED standard
    SEQEND      => 2,
    GENEID      => 3,
    SCORE       => 4, # << We don't use this column in our calculation
    DIRECTION   => 5,
    FIVESTART   => 6, # The last four columns do not follow BED format
    FIVEEND     => 7,
    THREESTART  => 8,
    THREEEND    => 9,
    LASTCOLUMN  => 5,
};

use Time::HiRes qw( time ); # Measures runtime for Seaweed and Star_Bundler


#<=== SET PARAMETERS ===>#

# my @species = ("nasonia vitripennis", "apis mellifera", "atta cephalotes", "bombyx mori", "acyrthosiphon pisum", "danaus plexippus", "heliconius melpomene", "tribolium castaneum");
# my @dbs = ("nasonia_vitripennis_core_18_71_1", "apis_mellifera_core_18_71_1", "atta_cephalotes_core_18_71_1", "bombyx_mori_core_18_71_1", "acyrthosiphon_pisum_core_18_71_1", "danaus_plexippus_core_18_71_1", "heliconius_melpomene_core_18_71_1", "tribolium_castaneum_core_18_71_3");
#Rbh file for 1-2, 1-3, 1-4, etc
# my @rbhs = ("/home/grannysmith/data/new_rbh/nasonia_apis_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_atta_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_bombyx_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_acrythosiphon_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_danaus_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_heliconiusmelpomene_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_tribolium_rbh_new.txt");

# Default parameters
my $species_string = "PlantA,PlantB";
my $db_directory = "../inputs/db_dir";
my $output_directory = "../outputs";

# my $sequence_length = 2000;
my $window_size = 60;
my $pseudo_orthologs = 0;

my $threshold_a = 78;
my $threshold_b = 100;

GetOptions(
    'species|s=s' => \$species_string,
    'dbdir|d=s' => \$db_directory,
    'outdir|o=s' => \$output_directory,
    # 'slength|l=i' => \$sequence_length,
    'wsize|w=i' => \$window_size,
    'pseudo|p' => \$pseudo_orthologs,
    'tha|a=i' => \$threshold_a,
    'thb|b=i' => \$threshold_b,
    ) or die "Usage: $0 --slength INTEGER([2000], 500, 5000) --wsize INTEGER([60], 30, 80, 100) --pseudo\n";

if($pseudo_orthologs){
    $output_directory = $output_directory . "/pseudo";
}

my $output_directory_extra = $output_directory . "/extra";

if (not -e $output_directory) {
    mkdir $output_directory, 0755;
}

if (not -e $output_directory_extra) {
    mkdir $output_directory_extra, 0755;
}

my @species_list = split(/,/, $species_string);

my $fn_output = "$output_directory/cnsOut_m_wB_ws" . $window_size ."_" . substr($species_list[0], 0, 6) . "_" . scalar(@species_list) . ( $pseudo_orthologs ? "_pseudo" : "") . ".txt";
my $fn_output2 = $fn_output . ".alignmax";
my $fn_output3 = $fn_output . ".tabular";
my $fn_log = $fn_output . ".log";

open my $outfile, ">$fn_output";
open my $outfile_alignmax, ">$fn_output2";
open my $outfile_tabular, ">$fn_output3";
open my $logfile, ">$fn_log";

print $logfile "Given " . scalar(@species_list) . " species:\n";
for(my $i=0; $i<scalar(@species_list); $i++){
    print $logfile "\t $species_list[$i]\n";
}
print $logfile "DB at $db_directory\n";

# my $outfile_fn = "testing_pseudo.txt";
# open my $outfile, ">$outfile_fn";

#<=== LOAD RBHS ===>#
my @rbh_results = ();

for(my $i=1; $i<scalar(@species_list); $i++)
{
    # open my $rbhs_data, "<$rbhs[$i]", or die "\nError: could not open rbh file '$rbhs[$i]'";
    # $_ = <$rbhs_data>;
    # while(<$rbhs_data>)
    # {
    #     chomp;
    #     my @split = split(/\t/);
        
    #     #if($split[2] ne "none")
    #     #{
    #         push(@{$rbh_results[$i]->{$split[0]}}, $split[1]);
    #     #}
    # }

    # my %rbhs;
    my $fn_rbh = "$db_directory/$species_list[$i]/rbhSearch_result.txt";

    open(my $fh, '<', $fn_rbh) or die "Could not open file '$fn_rbh'. \n$!";
    while (my $line = <$fh>) {
        chomp $line;
        my @array = split(/\t/, $line);
        next if $#array < 3;

        # There are 4 columns in the rbh output file:
        #   the 1st and last belong to the first species
        #   the 2nd and 3rd belong to the second species
        #   the 1st and 2nd are protein IDs
        #   the 3rd and last are gene IDs if known, "unknown" otherwise
        # If any of the gene IDs are "unknown", we try and extract it from the protein ID
        # by taking a substring before the first ".", if a "." exists in the protein ID
        foreach my $ii (2..3) {
            if ($array[$ii] eq 'unknown') {
                $array[$ii] = index($array[3-$ii], '.') == -1 ? $array[3-$ii] : substr($array[3-$ii], 0, index($array[3-$ii], '.'));
            }
        }
        # $array[2] = $array[2] eq 'unknown' ? substr($array[1], 0, index($array[1], '.')) : $array[2] ;
        # $array[3] = $array[3] eq 'unknown' ? substr($array[0], 0, index($array[0], '.')) : $array[3] ;

        # push @{ $rbhs{$array[3]} }, $array[2];
        push @{ $rbh_results[$i]->{ $array[3]} }, $array[2];
        # $rbhs{$array[3]} = $array[2];
    }
    close $fh;

}

# print $rbh_results[1]->{"AT1G01320"} . "\n";
# print $rbh_results[2]->{"AT1G01320"} . "\n";
# print Dumper ($rbh_results[1]);
# print Dumper ($rbh_results[2]);
# exit;

#<== SET PSEUDO ORTHOLOGS ==>#
if($pseudo_orthologs)
{
    for(my $i=1;$i<scalar(@species_list);$i++)
    {
        # my @useful_genes = ();
        # my @useful_rbhs = ();
        # my %rbh_numbers;
        
        # foreach my $rbh (keys %{$rbh_results[$i]})
        # {
        #     push(@useful_genes, $rbh);
        #     my $no_rbhs = 0;
        #     foreach my $ortholog (@{$rbh_results[$i]->{$rbh}})
        #     {
        #         push(@useful_rbhs, $ortholog);
        #         $no_rbhs++;
        #     }
        #     $rbh_numbers{$rbh} = $no_rbhs;
        # }
        
        # $rbh_results[$i] = ();
        
        # @useful_genes = shuffle(@useful_genes);
        # @useful_rbhs = shuffle(@useful_rbhs);
        
        # my $j = 0;
        # for(my $k=0;$k<scalar(@useful_genes);$k++)
        # {
        #     while($rbh_numbers{$useful_genes[$k]})
        #     {
        #         push(@{$rbh_results[$i]->{$useful_genes[$k]}}, $useful_rbhs[$j]);
        #         $j++;
        #         $rbh_numbers{$useful_genes[$k]}--;
        #     }
        # }

        print $logfile "Warning: Running with Pseudo Orthologs\n";

        # my %rbhs = $rbh_results[$i];

        my @useful_genes = ();
        my @useful_rbhs = ();
        
        foreach my $rbh (keys %{$rbh_results[$i]})
        {
            # if($rbh_results[$i]->{$rbh} ne "none")
            # {
                push(@useful_genes, $rbh);
                push(@useful_rbhs, $rbh_results[$i]->{$rbh});
            # }
        }
        
        $rbh_results[$i] = ();
        
        @useful_genes = shuffle(@useful_genes);
        # @useful_rbhs = shuffle(@useful_rbhs);
        
        for(my $j=0;$j<scalar(@useful_genes);$j++)
        {
            # push @{ $rbh_results[$i]->{ $useful_genes[$j] } }, $useful_rbhs[$j];
            $rbh_results[$i]->{ $useful_genes[$j] } = $useful_rbhs[$j];

        }
    }
}

# print Dumper ($rbh_results[1]);
# print Dumper ($rbh_results[2]);
# exit;

#<=== Load Sequence Location Data ===>#

my %sequence_info_lookup; # e.g. $sequence_info_lookup{"Niben101_Niben101Scf10386g00010"}[1]
my %all_geneids; # e.g. $all_geneids{"Niben101"} gives you all the gene IDs of Niben101
my %db_fasta;

foreach my $species (@species_list) {

    my $fn_main = "$db_directory/$species/PlantA.bed";
    my $fn_utr5 = "$db_directory/$species/PlantA\_utr5.bed";
    my $fn_utr3 = "$db_directory/$species/PlantA\_utr3.bed";

    my $fn_fasta = "$db_directory/$species/PlantA.fa";

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

# #<=== GET GENES ===>#
# my $local_db = get_sequence_database("ensembl_local");

# my @genes = ();
# for(my $i=0;$i<scalar(@species_list);$i++)
# {
#     $genes[$i] = $local_db->get_all_accessions($species_list[$i]);
# }

my @species_1_genes = @{$all_geneids{$species_list[0]}};

#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes
my $start_id = "NV10220";
my $begin = 1;

my $total = scalar @species_1_genes;
my $count = 0;

# foreach my $s1_gene_accession (@{$genes[0]})
foreach my $s1_gene_accession (@species_1_genes)
{

    $count++;

    #<-- If we're ready to begin, start with a gene -->#
    # print "\nGene: $s1_gene_accession";
    if($s1_gene_accession eq $start_id)
    {
        $begin = 1;
    }
    
    #<-- If we're beginning -->#
    if($begin)
    {
        #<-- Keep track of the number of comparisons made -->#
        my $no_comparisons = 0;
        my @curgene_comparisons = ();
        
        # #<-- Firstly, set up S1 gene -->#
        # my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species_list[0], $s1_gene_accession);

        # #<-- Positive strand error fixing -->#
        # my $gene_1_start = $sequence_length;
        # my $g1_strand = $gene_1_sequence->[0]->{"strand"};
        # $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
        
        # #<-- Upstream sequence length default -->#
        # my $species_1_length = $sequence_length;
        
        # #If we can't physically take that length of sequence, then we need to cut it down
        # #Don't need to bother doing this for negative strand really
        # if($gene_1_start < $species_1_length)
        # {
        #     $species_1_length = $gene_1_start;
        # }
        
        # #Take the upstream sequence
        # my $species_1_sequence = $local_db->get_sequence_by_location
        # (
        #     Sequences::Database::Relative_Location->new
        #     (
        #         identifier => $s1_gene_accession,
        #         offset => - $species_1_length,
        #         'length' => $species_1_length,
        #         stop_at_neighbours => 1,
        #         'anchor' => "5' end"
        #     )
        # );
        
        my $species_1_sequence = $db_fasta{$species_list[0]}->seq($s1_gene_accession);

        print $logfile "$count/$total";
        print $logfile "\t$s1_gene_accession(" . $sequence_info_lookup{$species_list[0] . "_" . $s1_gene_accession}[CHROMID] . ")";
        print $logfile "[". length($species_1_sequence);

        #If the length is above the window size
        # if(length($species_1_sequence->[0]->seq) >= $window_size)
        if(length($species_1_sequence) >= $window_size)
        {
            # print "\n$species_list[0] sequence ok!";
            #Remove IUPAC codes
            # my $sequence_one = $species_1_sequence->[0]->seq;
            my $sequence_one = $species_1_sequence;
            $sequence_one =~ s/[^(A|T|C|G)]/N/g;
            
            #Create the new sequences
            my $species_1_final = Datatypes::Sequence::Local->new_from_string
            (
                $sequence_one, $species_list[0],
            );
            
            #Set up the stuff we need for the bundler run with gene 1
            # my $rdb1 = Ensembl_Database_Parameters->new("dbname" => $dbs[0],
            #     "alias" => $species_list[0],
            #     "location" => "local");
            
            my $rdb1 = FASTA_Sequence_Database_Parameters->new(
                    "dbname" => $species_list[0],
                    "filename" => "$db_directory/" . $species_list[0] . "/PlantA.fa",
                    "alias" => $species_list[0],
                    "location" => "local",
                    "natural_species_name" => $species_list[0],
                );

            #Get sequence information
            # my ($g1_region, $g1_coordsys) = $local_db->get_gene_information_by_accession($species_list[0], $s1_gene_accession);
            # my $g1_5prime = $species_1_sequence->[0]->{"five_prime_pos"};
            # my $g1_3prime = $species_1_sequence->[0]->{"three_prime_pos"};
            # my $g1_strand = $species_1_sequence->[0]->{"strand"};
            
            my $g1_region = $sequence_info_lookup{$species_list[0] . "_" . $s1_gene_accession}[GENEID]; #CHROMID
            my $g1_coordsys = "not_in_use"; # chromosome, scaffold, contig
            my $g1_5prime = $sequence_info_lookup{$species_list[0] . "_" . $s1_gene_accession}[FIVESTART];
            my $g1_3prime = $sequence_info_lookup{$species_list[0] . "_" . $s1_gene_accession}[THREEEND];
            my $g1_strand = $sequence_info_lookup{$species_list[0] . "_" . $s1_gene_accession}[DIRECTION];
                        
            print $logfile "$g1_strand,$g1_5prime,$g1_3prime]";

            # Dumbing the Bundler
            $g1_strand = '+';
            $g1_3prime = $g1_strand eq '+' ? length($sequence_one) : 1;
            $g1_5prime = $g1_strand eq '+' ? 1 : length($sequence_one);

            print $logfile "<$g1_strand,$g1_5prime,$g1_3prime>\n";

            #Set genomic interval data
            my $s1_ginterval = Genomic_Interval->new("genome_db" => $rdb1,
            "region" => $g1_region,
            "five_prime_pos" => $g1_5prime,
            "three_prime_pos" => $g1_3prime,
            "strand" => $g1_strand eq '+' ? "positive" : "negative",
            "working_sequence" => "ref_sequence",
            "coord_sys_name" => $g1_coordsys);
            
            #Just in case
            # $s1_ginterval->{"gi_sequence"} = $species_1_sequence->[0]->seq;
            # $s1_ginterval->{"gi_sequence_repeatmasked"} = $species_1_sequence->[0]->{"masked_sequence"};
            $s1_ginterval->{"gi_sequence"} = $sequence_one;
            $s1_ginterval->{"gi_sequence_repeatmasked"} = $sequence_one;

            #Set s1 gene comparison information
            $curgene_comparisons[0]->{"species"} = $species_list[0];
            $curgene_comparisons[0]->{"db"} = $species_list[0];
            $curgene_comparisons[0]->{"gi"} = $s1_ginterval;
            # $curgene_comparisons[0]->{"gene_1_start"} = $gene_1_start;
            $curgene_comparisons[0]->{"gene_1_strand"} = $g1_strand;
            
            #Now that species 1 is done, we can do the same for all of the other species
            #Now do the conservation for this gene
            #For each species
            for(my $i=1;$i<scalar(@species_list);$i++)
            {
                #Which species is it?
                my $species_two = $species_list[$i];
                # my $db_two = $dbs[$i];
                
                # print "\nSpecies: " . $species_list[$i];
                
                #Get orthologues for this gene
                my @s2_genes = ();
                if(defined($rbh_results[$i]->{$s1_gene_accession}))
                {
                    @s2_genes = @{$rbh_results[$i]->{$s1_gene_accession}};
                }
                
                print $logfile "-s$i/" . scalar(@s2_genes);
                #Do we have any orthologues for this gene?
                if(scalar(@s2_genes) > 0)
                {
                    #Go through all s2 genes
                    FORS2: foreach my $s2_gene_accession (@s2_genes)
                    {
                        print $logfile "\t$s2_gene_accession(" . $sequence_info_lookup{$species_two . "_" . $s2_gene_accession}[CHROMID] . ")";

                        # print "\n\t$species_two ortholog: $s2_gene_accession";
                        
                        # #Start making sequences
                        # my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_two, $s2_gene_accession);
                        
                        # #Positive strand error fixing (some genes start at the beginning of chromosomes)
                        # my $gene_2_start = $sequence_length;
                        # my $g2_strand = $gene_2_sequence->[0]->{"strand"};
                        # $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
                        
                        # #Upstream sequence length defaults
                        # my $species_2_length = $sequence_length;
                        
                        # #If we can't physically take that length of sequence, then we need to cut it down
                        # #Don't need to bother doing this for negative strand really
                        # if($gene_2_start < $species_2_length)
                        # {
                        #     $species_2_length = $gene_2_start;
                        # }
                        
                        # #Take the upstream sequences
                        # my $species_2_sequence = $local_db->get_sequence_by_location
                        # (
                        # Sequences::Database::Relative_Location->new
                        # (
                        # identifier => $s2_gene_accession,
                        # offset => - $species_2_length,
                        # 'length' => $species_2_length,
                        # stop_at_neighbours => 1,
                        # 'anchor' => "5' end"
                        # )
                        # );
                        my $species_2_sequence = $db_fasta{$species_two}->seq($s2_gene_accession);
                        
                        print $logfile "[". length($species_1_sequence) . "]";
                        # if(length($species_2_sequence->[0]->seq) >= $window_size)
                        if(length($species_2_sequence) >= $window_size)
                        {
                            #Remove IUPAC codes
                            my $sequence_two = $species_2_sequence;
                            $sequence_two =~ s/[^(A|T|C|G)]/N/g;
                            
                            #If we're ready to conduct the conservation analysis, then...
                            #Create the sequences
                            my $species_2_final = Datatypes::Sequence::Local->new_from_string
                            (
                            $sequence_two, $species_two,
                            );
                            
                            #Set up the job
                            my $job = Jobs::Subtasks::Seaweed_Job->new
                            (
                            sequence_1 => $species_1_final,
                            sequence_2 => $species_2_final,
                            windowsize => $window_size,
                            );
                            
                            #And run the job
                            # my $result = $job->run;
                            my $result;
                            try{
                                my $start_seaweed = time();
                                $result = $job->run;
                                my $end_seaweed = time();
                                printf $logfile "[ST%.2f]", $end_seaweed - $start_seaweed ;
                            } catch {
                                print $logfile "[ST_aborted]\t";
                                # print $outfile "Seaweed_aborted\n";
                                # print $outfile_alignmax "-1\n";
                                # print $outfile_tabular "Seaweed_aborted\n";
                                next FORS2;
                            };

                            # In case there are any hanging processes from running the Seaweed executable
                            # We log the kill count at the end of the script because of <defunct> processes after kill
                            system('killall -9 AlignmentPlot_posix_default_release 2>/dev/null');

                            # This removes any core dumps may be generated by seaweed
                            # By adding "ulimit -c 0" in the wrapper script (cyverse),
                            # this should not be necessary, but we keep this here in case
                            # this script is called without the wrapper.
                            if ( glob("core.*") ) { 
                                system('rm core.* 2>/dev/null');
                                print $logfile "c";
                            }

                            #What is the max alignment value for this run of seaweeds?
                            my $alignmax = 0;
                            foreach my $cur_val (values %{$result->{"plot"}})
                            {
                                if($cur_val > $alignmax)
                                {
                                    $alignmax = $cur_val;
                                }
                            }
                            my $alignment_results = convert_gs_seaweed_result_to_ap($result->plot);
                            #Convert alignment results to APPLES format
                            #push(@{$alignment_results}, convert_gs_seaweed_result_to_ap($result->plot));
                            
                            # my $rdb2 = Ensembl_Database_Parameters->new("dbname" => $db_two,
                            # "alias" => $species_two,
                            # "location" => "local");

                            my $rdb2 = FASTA_Sequence_Database_Parameters->new(
                                "dbname" => $species_two,
                                "filename" => "$db_directory/" . $species_two . "/PlantA.fa",
                                "alias" => $species_two,
                                "location" => "local",
                                "natural_species_name" => $species_two,
                            );

                            #Get sequence information for the genomic intervals
                            # my ($g2_region, $g2_coordsys) = $local_db->get_gene_information_by_accession($species_two, $s2_gene_accession);
                            my $g2_region = $sequence_info_lookup{$species_two . "_" . $s2_gene_accession}[GENEID];
                            my $g2_coordsys = "not_in_use";

                            # my $g2_5prime = $species_2_sequence->[0]->{"five_prime_pos"};
                            # my $g2_3prime = $species_2_sequence->[0]->{"three_prime_pos"};
                            # my $g2_strand = $species_2_sequence->[0]->{"strand"};
                            
                            my $g2_5prime = $sequence_info_lookup{$species_two . "_" . $s2_gene_accession}[FIVESTART];
                            my $g2_3prime = $sequence_info_lookup{$species_two . "_" . $s2_gene_accession}[THREEEND];
                            my $g2_strand = $sequence_info_lookup{$species_two . "_" . $s2_gene_accession}[DIRECTION];

                            print $logfile "[$g2_strand,$g2_5prime,$g2_3prime]";

                            $g2_strand = '+';

                            $g2_3prime = $g2_strand eq '+' ? length($sequence_two) : 1;
                            $g2_5prime = $g2_strand eq '+' ? 1 : length($sequence_two);

                            print $logfile "<$g2_strand,$g2_5prime,$g2_3prime>";

                            #Set genomic interval data
                            my $s2_ginterval = Genomic_Interval->new("genome_db" => $rdb2,
                            "region" => $g2_region,
                            "five_prime_pos" => $g2_5prime,
                            "three_prime_pos" => $g2_3prime,
                            "strand" => $g2_strand eq '+' ? "positive" : "negative",
                            "working_sequence" => "ref_sequence",
                            "coord_sys_name" => $g2_coordsys);
                            
                            #Just in case
                            # $s2_ginterval->{"gi_sequence"} = $species_2_sequence->[0]->seq;
                            # $s2_ginterval->{"gi_sequence_repeatmasked"} = $species_2_sequence->[0]->{"masked_sequence"};
                            $s2_ginterval->{"gi_sequence"} = $sequence_two;
                            $s2_ginterval->{"gi_sequence_repeatmasked"} = $sequence_two;

                            $no_comparisons++;
                            
                            #Set s2 gene
                            $curgene_comparisons[$no_comparisons]->{"species"} = $species_two;
                            $curgene_comparisons[$no_comparisons]->{"db"} = $species_two;
                            $curgene_comparisons[$no_comparisons]->{"gi"} = $s2_ginterval;
                            $curgene_comparisons[$no_comparisons]->{"alignment"} = $alignment_results;
                            # $curgene_comparisons[$no_comparisons]->{"gene_2_start"} = $gene_2_start;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_acc"} = $s2_gene_accession;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_strand"} = $g2_strand;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_alignmax"} = $alignmax;
                        }
                    }
                }

                print $logfile "\n";

            }
        } else { # length($species_1_sequence) >= $window_size
            print $logfile "]\n";
        }

        my @arr = ();
        
        #What is going on here
        
        for(my $i=1;$i<scalar(@curgene_comparisons);$i++)
        {
            my $cur_g = $curgene_comparisons[$i];
            push(@arr, $cur_g->{"gi"});
            #print "\n" . $cur_g->{"species"};
            #print "\n" . $cur_g->{"db"};
        }
        
        if(scalar(@arr) > 0)
        {
            #Make an interval set for comparison
            my $interval_set = Genomic_Interval_Set->new("genomic_interval_set" => \@arr);
            
            #Make an evolutionary tree (will not be used yet)
            my $evtree = Evolutionary_Tree->new("root" => $species_list[0],
            "sub_trees" => []);
            
            #Make the parameters
            my $ptm = Partial_Threshold_Matrix->new("evolutionary_tree" => $evtree);
            
            my $wpap = Seaweed_Algorithm_Parameters->new("stepwidth" => 1,
            "windowlength" => $window_size,
            "cutoff_for_uninteresting_alignments" => 57);
            
            my $sbp = Star_Bundler_Parameters->new("overlap_tolerance" => 20,
            "belief_value" => 0.05,
            "threshold_a" => $threshold_a,
            "threshold_b" => $threshold_b,
            "partial_threshold_matrix" => $ptm);
            
            my $sp = Sequence_Parameters->new("region" => "upstream",
            "min_length_to_return" => 50,
            "max_length_to_search" => 10000,
            );
            
            my $parameters = ReMo_Set_Phylogenetic_Constructor_Parameters->new("window_pair_algorithm_parameters" => $wpap,
            "star_bundler_parameters" => $sbp,
            "sequence_parameters" => $sp,
            "sequence_databases_to_use_for_homologs" => []);
            
            my $alignment_results;

            for(my $i=1;$i<scalar(@curgene_comparisons);$i++)
            {
                my $cur_g = $curgene_comparisons[$i];

                push(@{$alignment_results}, $cur_g->{"alignment"});
            }

            my $bundler = Star_Bundler->new();
            my $start_bundler = time();
            my @remo_sets = $bundler->truncated_bundle($curgene_comparisons[0]->{"gi"}, $interval_set, $parameters, $alignment_results);
            my $end_bundler = time();

            printf $logfile "-b\t[BT%.2f]", $end_bundler - $start_bundler;
            
            print $outfile "--PairStart\n";
            print $outfile $s1_gene_accession . "\n"; #file2, # filename-i.fa
            print $outfile_tabular "$s1_gene_accession\t";
            # print $outfile $curgene_comparisons[0]->{"gene_1_start"} . "\n";
            print $outfile $curgene_comparisons[0]->{"gene_1_strand"} . "\n";
            
            for(my $i=1;$i<scalar(@curgene_comparisons);$i++)
            {
                my $cur_g = $curgene_comparisons[$i];
                
                print $outfile $cur_g->{"species"} . "\n";
                print $outfile $cur_g->{"gene_2_acc"} . "\n";
                print $outfile_tabular $cur_g->{"gene_2_acc"} . "\t";
                # print $outfile $cur_g->{"gene_2_start"} . "\n";
                print $outfile $cur_g->{"gene_2_strand"} . "\n";
                print $outfile $cur_g->{"gene_2_alignmax"} . "\n";
                print $outfile_alignmax $cur_g->{"gene_2_alignmax"} . "\n";
            }
            
            my $extra_file_counter = 1;

            # print "open >$output_directory_extra/$s1_gene_accession-$extra_file_counter.fa";
            print $outfile "--RemoStart\n";
            foreach my $remo_set (@remo_sets)
            {
                # $s1_gene_accession-i.fa
                my $s1_gene_accession_filtered =~ s/[^A-Za-z0-9\-\.\_]//g;
                open my $outfile_extra, ">$output_directory_extra/$s1_gene_accession-$extra_file_counter.fa";
                # print ">$output_directory_extra/$s1_gene_accession-$extra_file_counter.fa";
                $extra_file_counter++;

                print $outfile "--InnerStart\n";
                foreach my $remo (@{$remo_set->{"remo_set"}})
                {
                    print $outfile "--Organism: " . $remo->{"genome_db"}->{"alias"} . "\n"; # line1 > species_name
                    print $outfile_extra ">" . $remo->{"genome_db"}->{"alias"} . " " . $remo->{"conservation"} . " " . $remo->{"belief_score"} . "\n";
                    print $outfile $remo->{"gi_sequence"} . "\n"; # file2, #line2
                    print $outfile_extra $remo->{"gi_sequence"} . "\n";
                    print $outfile_tabular $remo->{"gi_sequence"} . "\t";
                    print $outfile $remo->{"strand"} . "\n";
                    print $outfile $remo->{"five_prime_pos"} . "\n";
                    print $outfile $remo->{"three_prime_pos"} . "\n";
                    print $outfile $remo->{"conservation"} . "\n"; # file2,#line1 
                    print $outfile $remo->{"belief_score"} . "\n"; # file2,#line1
                    print $outfile_tabular $remo->{"conservation"} . "\t";
                    print $outfile_tabular $remo->{"belief_score"} . "\t";
                    print $outfile $remo->{"repeat_ratio"} . "\n";
                }
                print $outfile "--InnerEnd\n";

                close($outfile_extra);
                # close
            }
            print $outfile "--RemoEnd\n";
            print $outfile "--PairEnd\n";
            print $outfile_tabular "\n"; #file2 end
            
            foreach my $remo_set (@remo_sets)
            {
                print "\nReMo:";
                foreach my $remo (@{$remo_set->{"remo_set"}})
                {
                    print "\n\t" . $remo->{"belief_score"} . " - " . $remo->{"genome_db"}->{"alias"};
                }
            }
        }
    }
}

my $kill_count = `ps | awk '/AlignmentPlot*/ && /defunct/ && !/awk/' | wc -l`;
print $logfile "KCount: $kill_count";
