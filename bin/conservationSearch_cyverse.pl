#!/usr/bin/perl

# derived form LB_two_species_conservation_updated_1.pl
# additional custom library path

use strict;

# BASEDIR contains path to ORANGES base directory
use lib "../webseaweeds/";

# This is the ORANGES Loader for APPLES. Go fruit.
BEGIN {
    use Configuration::AppleSeeds;
    Configuration::AppleSeeds::load_APPLES();
    1;
}

#General
use JSON;
use NDU::Gen qw(:ALL);
use Data::Dumper;

use File::Spec;

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


#<=== PART 1: Set Parameters ===>#

#Two species that you're going to be comparing
# my $species_1 = "plantA";#"oryza sativa";
# my $species_2 = "plantB";#"arabidopsis thaliana";

my $species_1 = "PlantA";
my $species_2 = "PlantB";

$species_1 =~ s/\s//g; # remove any spaces in name
$species_2 =~ s/\s//g;

my @species_list = ($species_1, $species_2);

# Default Parameters
my $sequence_length = 2000;#How much upstream sequence to take
my $window_size = 60;#Window size for seaweeds algorithm
my $pseudo_orthologs = 0;

GetOptions(
    'slength|l=i' => \$sequence_length,
    'wsize|w=i' => \$window_size,
    'pseudo|p' => \$pseudo_orthologs,
    ) or die "Usage: $0 --slength INTEGER([2000], 500, 5000) --wsize INTEGER([60], 30, 80, 100) --pseudo\n";

my $fn_output = "../outputs/conservation_result_wB_ws" . $window_size ."_" . substr($species_1, 0, 6) . "_" . substr($species_2, 0, 6) . ( $pseudo_orthologs ? "_pseudo" : "") . ".txt";
my $fn_log = $fn_output . ".log";
# my $fn_output = "../output/conservation_result_wB_pd_two_species_plantV_plantA_long.txt";
open my $outfile, ">$fn_output";
open my $logfile, ">$fn_log";

#<=== PART 2: Load Inputs ===>#

#<=== Load RBH Data ===>#
# my $rbh_file = "../output/rbhSearchForked_result_plantV_plantA.txt"; # short version

my %rbhs;
my $fn_rbh = "../inputs/rbhSearch_result_$species_1\_$species_2.txt";

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

    push @{ $rbhs{$array[3]} }, $array[2];
    # $rbhs{$array[3]} = $array[2];
}
close $fh;

#<== SET PSEUDO ORTHOLOGS ==>#
if($pseudo_orthologs)
{
    my @useful_genes = ();
    my @useful_rbhs = ();
    
    foreach my $rbh (keys %rbhs)
    {
        if($rbhs{$rbh} ne "none")
        {
            push(@useful_genes, $rbh);
            push(@useful_rbhs, $rbhs{$rbh});
        }
    }
    
    %rbhs = ();
    
    @useful_genes = shuffle(@useful_genes);
    @useful_rbhs = shuffle(@useful_rbhs);
    
    for(my $i=0;$i<scalar(@useful_genes);$i++)
    {
        $rbhs{$useful_genes[$i]} = $useful_rbhs[$i];
    }
}
# print Dumper (%rbhs);
# exit;

#<=== Load Sequence Location Data ===>#

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

# #<=== GET GENES ===>#
# my $local_db = get_sequence_database("ensembl_local");

# my @species_1_genes = @{$local_db->get_all_accessions($species_1)};
# my @species_2_genes = @{$local_db->get_all_accessions($species_2)};

my @species_1_genes = @{$all_geneids{$species_1}};
my @species_2_genes = @{$all_geneids{$species_2}};


#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes

my $start_id = "Niben101Scf02104g00001";
my $begin = 1; # 0 to begin with the gene specified in $start_id, 1 otherwise.

my $total = scalar @species_1_genes;
my $count = 0;

foreach my $s1_gene_accession (@species_1_genes)
{
    $count++;

    if($s1_gene_accession eq $start_id)
    {
        $begin = 1;
    }
    
    if($begin)
    {
        # if( defined($rbhs{$s1_gene_accession}) && $rbhs{$s1_gene_accession} ne "none" ) {
        if( defined($rbhs{$s1_gene_accession}) ) {
            foreach my $s2_gene_accession ( @{$rbhs{$s1_gene_accession}} ) {
                if(  $s2_gene_accession ne "none" ) {
                
                    print "Progress: ". $count . "/" . $total . " s1_gene_accession:" . $s1_gene_accession . "\n";
                    #Get the RBH in species 2
                    # my $s2_gene_accession = $rbhs{$s1_gene_accession};
                    
                    print $logfile "$count/$total";
                    print $logfile "\t$s1_gene_accession(" . $sequence_info_lookup{$species_1 . "_" . $s1_gene_accession}[CHROMID] . ")";
                    print $logfile "\t$s2_gene_accession(" . $sequence_info_lookup{$species_2 . "_" . $s2_gene_accession}[CHROMID] . ")";
                    
                    #Now we have to check if we can take the sequence we want to take
                    #So, get the sequences
                    # my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species_1, $s1_gene_accession);
                    # my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_2, $s2_gene_accession);

                    # #        print Dumper($gene_2_sequence);
                    # #Where do they start on the chromosome?
                    # #Positive strand error fixing (some genes start at the beginning of chromosomes)
                    # my $gene_1_start = $sequence_length;
                    # my $gene_2_start = $sequence_length;
                    # my $g1_real_start;
                    # my $g2_real_start;
                    # my $g1_strand = $gene_1_sequence->[0]->{"strand"};
                    # my $g2_strand = $gene_2_sequence->[0]->{"strand"};
                    
                    # $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
                    # $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
                    
                    # $g1_real_start = $gene_1_sequence->[0]->{"five_prime_pos"};
                    # $g2_real_start = $gene_2_sequence->[0]->{"five_prime_pos"};
                    
                    # #Upstream sequence length defaults
                    # my $species_1_length = $sequence_length;
                    # my $species_2_length = $sequence_length;
                    
                    # #If we can't physically take that length of sequence, then we need to cut it down
                    # if($gene_1_start < $species_1_length){
                    #     $species_1_length = $gene_1_start;
                    # }
                    # if($gene_2_start < $species_2_length){
                    #     $species_2_length = $gene_2_start;
                    # }
                    
                    # #Take both of the upstream sequences
                    # my $species_1_sequence = $local_db->get_sequence_by_location(
                    #     Sequences::Database::Relative_Location->new(
                    #         identifier => $s1_gene_accession,
                    #         offset => - $species_1_length,
                    #         'length' => $species_1_length,
                    #         stop_at_neighbours => 1,
                    #         'anchor' => "5' end"
                    #     )
                    # );
                    
                    # my $species_2_sequence = $local_db->get_sequence_by_location
                    # (
                    #     Sequences::Database::Relative_Location->new
                    #     (
                    #         identifier => $s2_gene_accession,
                    #         offset => - $species_2_length,
                    #         'length' => $species_2_length,
                    #         stop_at_neighbours => 1,
                    #         'anchor' => "5' end"
                    #     )
                    # );
                    
                    my $species_1_sequence = $db_fasta{$species_1}->seq($s1_gene_accession);
                    my $species_2_sequence = $db_fasta{$species_2}->seq($s2_gene_accession);

                    print $logfile "\tLen-(". length($species_1_sequence) . "," . length($species_2_sequence) . ")";

                    #Both sequences must be above or equal to the window size
                    # if(length($species_1_sequence->[0]->seq) >= $window_size && length($species_2_sequence->[0]->seq) >= $window_size)
                    if( length($species_1_sequence) >= $window_size &&
                        length($species_2_sequence) >= $window_size &&
                        $species_1_sequence ne "" &&
                        $species_2_sequence ne "")
                    {
                        #Remove IUPAC codes
                        my $sequence_one = $species_1_sequence;
                        $sequence_one =~ s/[^ATCG]/N/g;
                        my $sequence_two = $species_2_sequence;
                        $sequence_two =~ s/[^ATCG]/N/g;

                        #If we're ready to conduct the conservation analysis, then...
                        #Create the sequences
                        my $species_1_final = Datatypes::Sequence::Local->new_from_string(
                                $sequence_one,
                                $species_1,
                            );
                        my $species_2_final = Datatypes::Sequence::Local->new_from_string(
                                $sequence_two,
                                $species_2,
                            );
                        
                        
                        #Set up the job
                        my $job = Jobs::Subtasks::Seaweed_Job->new(
                                sequence_1 => $species_1_final,
                                sequence_2 => $species_2_final,
                                windowsize => $window_size,
                                nprocess => 8,
                            );

                        # my @db_handles;
                        # # my $ensembl_registry = 'Bio::EnsEMBL::Registry';
                        # my $ensembl_registry = $local_db->{"registry"};
                        # # print Dumper($local_db);
                        
                        # my @db_adaptors = @{ $ensembl_registry->get_all_DBAdaptors() };
                        # foreach my $db_adaptor (@db_adaptors) {
                        #     # print "db_adaptor: " . Dumper($db_adaptor);
                        #     my $db_connection = $db_adaptor->dbc();
                        #     # push @db_handles, $db_connection->db_handle();
                        #     $job->add_handle($db_connection->db_handle());

                        #     # push @db_handles, $db_connection->db_handle();
                        #     # print "script: $db_connection->db_handle()". Dumper($db_connection->db_handle()) . "\n";
                        # }
                        # # print "db_handles:" . Dumper (@db_handles);

                        
                        #And run the job
                        my $start_seaweed = time();
                        my $result = $job->run;
                        my $end_seaweed = time();
                        printf $logfile "\t[ST%.2f]", $end_seaweed - $start_seaweed ;

                        system('killall -9 AlignmentPlot_posix_default_release');
                        if ( $? != -1) {
                            print $logfile "(kill)";
                        }

                        # print "\nConservation script two species, removing temporary files.\n";
                        # unlink glob "'./tempfiles/*'";
                        # print "\nConservation script two species, seaweed tempfiles removed.\n";
                        
                        #What is the max alignment value for this run of seaweeds?
                        my $alignmax = 0;
                        foreach my $cur_val (values %{$result->{"plot"}}){
                            if($cur_val > $alignmax){
                                $alignmax = $cur_val;
                            }
                        }
                        
                        my $alignment_results;

                        #Convert alignment results to APPLES format
                        push(@{$alignment_results}, convert_gs_seaweed_result_to_ap($result->plot));
                        
                        #<==Bundling==>#
                        # my $rdb1 = Ensembl_Database_Parameters->new(
                        #         "dbname" => $db_1,
                        #         "alias" => $species_1,
                        #         "location" => "local",
                        #     );
                        
                        # my $rdb2 = Ensembl_Database_Parameters->new(
                        #         "dbname" => $db_2,
                        #         "alias" => $species_2,
                        #         "location" => "local",
                        #     );

                        my $rdb1 = FASTA_Sequence_Database_Parameters->new(
                                "dbname" => "no_dbname1",
                                "filename" => File::Spec->rel2abs("../inputs/$species_1\.fa"), #_full.fa
                                "alias" => $species_1,
                                "location" => "local",
                                "natural_species_name" => $species_1,
                            );
                        
                        my $rdb2 = FASTA_Sequence_Database_Parameters->new(
                                "dbname" => "no_dbname2",
                                "filename" => File::Spec->rel2abs("../inputs/$species_2\.fa"), #_full.fa
                                "alias" => $species_2,
                                "location" => "local",
                                "natural_species_name" => $species_2,
                            );
                        
                        #Get sequence information for the genomic intervals
                        # my ($g1_region, $g1_coordsys) = $local_db->get_gene_information_by_accession($species_1, $s1_gene_accession);
                        # my ($g2_region, $g2_coordsys) = $local_db->get_gene_information_by_accession($species_2, $s2_gene_accession);
                        
                        my $g1_region = $sequence_info_lookup{$species_1 . "_" . $s1_gene_accession}[GENEID]; #CHROMID
                        my $g1_coordsys = "not_in_use"; # chromosome, scaffold, contig
                        my $g2_region = $sequence_info_lookup{$species_2 . "_" . $s2_gene_accession}[GENEID]; #CHROMID
                        my $g2_coordsys = "not_in_use";

                        # my $g1_5prime = $species_1_sequence->[0]->{"five_prime_pos"};
                        # my $g1_3prime = $species_1_sequence->[0]->{"three_prime_pos"};
                        # my $g1_strand = $species_1_sequence->[0]->{"strand"};
                        
                        # my $g2_5prime = $species_2_sequence->[0]->{"five_prime_pos"};
                        # my $g2_3prime = $species_2_sequence->[0]->{"three_prime_pos"};
                        # my $g2_strand = $species_2_sequence->[0]->{"strand"};

                        my $g1_5prime = $sequence_info_lookup{$species_1 . "_" . $s1_gene_accession}[FIVESTART];
                        my $g1_3prime = $sequence_info_lookup{$species_1 . "_" . $s1_gene_accession}[THREEEND];
                        my $g1_strand = $sequence_info_lookup{$species_1 . "_" . $s1_gene_accession}[DIRECTION];
                        
                        my $g2_5prime = $sequence_info_lookup{$species_2 . "_" . $s2_gene_accession}[FIVESTART];
                        my $g2_3prime = $sequence_info_lookup{$species_2 . "_" . $s2_gene_accession}[THREEEND];
                        my $g2_strand = $sequence_info_lookup{$species_2 . "_" . $s2_gene_accession}[DIRECTION];

                        print $logfile "\t<$g1_strand/$g2_strand,$g1_5prime/$g2_5prime,$g1_3prime/$g2_3prime,". length($sequence_one) . "/" . length($sequence_two) .">";

                        $g1_strand = '+';
                        $g2_strand = '+';

                        $g1_3prime = $g1_strand eq '+' ? length($sequence_one) : 1;
                        $g1_5prime = $g1_strand eq '+' ? 1 : length($sequence_one);
                        $g2_3prime = $g2_strand eq '+' ? length($sequence_two) : 1;
                        $g2_5prime = $g2_strand eq '+' ? 1 : length($sequence_two);

                        print $logfile "\tg1(". length($sequence_one) .")$g1_strand\[5:$g1_5prime,3:$g1_3prime]";
                        print $logfile "\tg2(". length($sequence_two) .")$g2_strand\[5:$g2_5prime,3:$g2_3prime]";

                        #Set genomic interval data
                        my $s1_ginterval = Genomic_Interval->new(
                                "genome_db" => $rdb1,
                                "region" => $g1_region,
                                "five_prime_pos" => $g1_5prime,
                                "three_prime_pos" => $g1_3prime,
                                "strand" => $g1_strand eq '+' ? "positive" : "negative",
                                "working_sequence" => "ref_sequence",
                                "coord_sys_name" => $g1_coordsys,
                            );
                        
                        #Just in case
                        # $s1_ginterval->{"gi_sequence"} = $species_1_sequence->[0]->seq;
                        # $s1_ginterval->{"gi_sequence_repeatmasked"} = $species_1_sequence->[0]->{"masked_sequence"};

                        $s1_ginterval->{"gi_sequence"} = $sequence_one;
                        $s1_ginterval->{"gi_sequence_repeatmasked"} = $sequence_one;
                        
                        #Set genomic interval data
                        my $s2_ginterval = Genomic_Interval->new(
                                "genome_db" => $rdb2,
                                "region" => $g2_region,
                                "five_prime_pos" => $g2_5prime,
                                "three_prime_pos" => $g2_3prime,
                                "strand" => $g2_strand eq '+' ? "positive" : "negative",
                                "working_sequence" => "ref_sequence",
                                "coord_sys_name" => $g2_coordsys,
                            );
                        
                        #Just in case
                        # $s2_ginterval->{"gi_sequence"} = $species_2_sequence->[0]->seq;
                        # $s2_ginterval->{"gi_sequence_repeatmasked"} = $species_2_sequence->[0]->{"masked_sequence"};

                        $s2_ginterval->{"gi_sequence"} = $sequence_two;
                        $s2_ginterval->{"gi_sequence_repeatmasked"} = $sequence_two;
                        
                        #Make an interval set for comparison
                        my @arr = ($s2_ginterval);
                        my $interval_set = Genomic_Interval_Set->new(
                                "genomic_interval_set" => \@arr,
                            );
                        
                        #Make an evolutionary tree (will not be used yet)
                        my $evtree1 = Evolutionary_Tree->new(
                                "root" => $species_1,
                                "sub_trees" => [],
                            );
                        my $evtree2 = Evolutionary_Tree->new(
                                "root" => $species_2,
                                "sub_trees" => [],
                            );
                        my $evtree3 = Evolutionary_Tree->new(
                                "root" => "",
                                "sub_trees" => [$evtree1, $evtree2],
                            );
                        
                        #Make the parameters
                        my $ptm = Partial_Threshold_Matrix->new(
                                "evolutionary_tree" => $evtree3,
                            );
                        
                        my $wpap = Seaweed_Algorithm_Parameters->new(
                                "stepwidth" => 1,
                                "windowlength" => $window_size,
                                "cutoff_for_uninteresting_alignments" => 57
                            );
                        my $sbp = Star_Bundler_Parameters->new(
                                "overlap_tolerance" => 20,
                                "belief_value" => 0.05,
                                "partial_threshold_matrix" => $ptm,
                            ); # LB - windowlength was set to 50 here, changed it to 60
                        
                        my $sp = Sequence_Parameters->new(
                                "region" => "upstream",
                                "min_length_to_return" => 50,
                                "max_length_to_search" => 10000,
                            );
                        
                        my $parameters = ReMo_Set_Phylogenetic_Constructor_Parameters->new(
                                "window_pair_algorithm_parameters" => $wpap,
                                "star_bundler_parameters" => $sbp,
                                "sequence_parameters" => $sp,
                                "sequence_databases_to_use_for_homologs" => [],
                            );
                        
                        my $bundler = Star_Bundler->new();

                        my $start_bundler = time();
                        my @remo_sets = $bundler->truncated_bundle($s1_ginterval, $interval_set, $parameters, $alignment_results);
                        my $end_bundler = time();

                        printf $logfile "\t[BT%.2f]", $end_bundler - $start_bundler;

                        #Print the appropriate output
                        print $outfile "--PairStart\n";
                        print $outfile $s1_gene_accession . "\n";
                        print $outfile $s2_gene_accession . "\n";
                        print $outfile $g1_strand . "\n";
                        print $outfile $g2_strand . "\n";
                        print $outfile '$g1_real_start' . "\n";
                        print $outfile '$g2_real_start' . "\n";
                        print $outfile $species_1_sequence . "\n";
                        print $outfile $species_2_sequence . "\n";
                        print $outfile $alignmax . "\n";
                        print $outfile "--RemoStart\n";
                        foreach my $remo_set (@remo_sets)
                        {
                            print $outfile "--InnerStart\n";
                            foreach my $remo (@{$remo_set->{"remo_set"}})
                            {
                                print $outfile "--Organism: " . $remo->{"genome_db"}->{"alias"} . "\n";
                                print $outfile $remo->{"gi_sequence"} . "\n";
                                print $outfile $remo->{"strand"} . "\n";
                                print $outfile $remo->{"five_prime_pos"} . "\n";
                                print $outfile $remo->{"three_prime_pos"} . "\n";
                                print $outfile $remo->{"conservation"} . "\n";
                                print $outfile $remo->{"belief_score"} . "\n";
                                print $outfile $remo->{"repeat_ratio"} . "\n";
                            }
                            print $outfile "--InnerEnd\n";
                        }
                        print $outfile "--RemoEnd\n";
                        print $outfile "--PairEnd\n";
                                        
                        print "\n$s1_gene_accession -> $s2_gene_accession ($alignmax)\n";
                        
                        foreach my $remo_set (@remo_sets)
                        {
                            print "\n\t" . $remo_set->{"remo_set"}->[0]->{"belief_score"} . " - " . $remo_set->{"remo_set"}->[1]->{"belief_score"};
                        }
                        print "-.-\n";
                    }
                    #print "\nGen:" . substr($gene_1_sequence->[0]->seq, 0, 10);
                    #print "\n" . Dumper($species_1_sequence->[0]->seq);
                    #exit;
                    print $logfile "\n";
                } # if(  $s2_gene_accession ne "none" )
            } #foreach my $s2_gene_accession ( @{$rbhs{$s1_gene_accession}} )
        } #if(defined($rbhs{$s1_gene_accession}))
        # } #if(defined($rbhs{$s1_gene_accession}) && $rbhs{$s1_gene_accession} ne "none")
    } #if($begin)
}
