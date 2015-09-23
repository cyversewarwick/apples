#!/usr/bin/perl

use strict;

# BASEDIR contains path to ORANGES base directory
# use lib "../../";
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

# use Devel::PrettyTrace;

#Grannysmith
use Runtime;
use Sequences::Database::Relative_Location;
use Sequences::Database::Sequence::Ensembl;
use Sequences::Database::Sequence::Genbank;
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

#<=== SET PARAMETERS ===>#
my $species_1 = "vitis vinifera";
my $species_2 = "arabidopsis thaliana";
my $species_3 = "populus trichocarpa";
my $species_4 = "medicago truncatula";

# my $species_1 = "grape";
# my $species_2 = "arabidopsis";
# my $species_3 = "poplar";
# my $species_4 = "medicago_truncatula";

my $db_1 = "vitis_vinifera_core_20_73_3";
my $db_2 = "arabidopsis_thaliana_core_20_73_10";
my $db_3 = "populus_trichocarpa_core_20_73_20";
my $db_4 = "medicago_truncatula_core_20_73_1";

# my $db_1 = "vitis_vinifera_core_20_81_3";
# my $db_2 = "arabidopsis_thaliana_core_20_81_10";
# my $db_3 = "populus_trichocarpa_core_20_81_20";
# my $db_4 = "medicago_truncatula_core_20_81_1";

# my $rbh_1_2 = "/home/grannysmith/data/RBH/dummy_VvAt_rbh.txt";
# my $rbh_1_3 = "/home/grannysmith/data/RBH/dummy_VvPt_rbh.txt";
# my $rbh_1_4 = "/home/grannysmith/data/RBH/dummy_VvMt_rbh.txt";



# my $species_1 = "plantV";
# my $species_2 = "plantM";
# my $species_3 = "plantA";
# my $species_4 = "plantP";

# my $db_1 = "plantV_core";
# my $db_2 = "plantM_core";
# my $db_3 = "plantA_core";
# my $db_4 = "plantP_core";

my $rbh_1_2 = "../output/rbhSearchForked_result_plantV_plantA.txt";
my $rbh_1_3 = "../output/rbhSearchForked_result_plantV_plantP.txt";
my $rbh_1_4 = "../output/rbhSearchForked_result_plantV_plantM.txt";

my $sequence_length = 2000;
my $window_size = 50;

my $pseudo_orthologs = 0;

my $outfile_fn = "../conservation_result_bundled_plantV_plantM_plantA_plantP.txt";
open my $outfile, ">$outfile_fn";

my %rbhs;

#<=== LOAD RBHS ===>#
my %rbh_hash_1_2 = ();
my %rbh_hash_1_3 = ();
my %rbh_hash_1_4 = ();

my @rbh_files = ($rbh_1_2, $rbh_1_3, $rbh_1_4);
for(my $i=0;$i<3;$i++)
{
    open my $rbhs_data, "<$rbh_files[$i]", or die "\nError: Could not open rbh file";
    $_ = <$rbhs_data>;
    while(<$rbhs_data>)
    {
        chomp;
        my @split = split(/\t/);
        if($i==0)
        {
            $rbh_hash_1_2{$split[3]} = $split[2];
        }
        elsif($i==1)
        {
            $rbh_hash_1_3{$split[3]} = $split[2];
        }
        elsif($i==2)
        {
            $rbh_hash_1_4{$split[3]} = $split[2];
        }
    }
    close $rbhs_data;
}

#<== SET PSEUDO ORTHOLOGS ==>#
if($pseudo_orthologs)
{
    for(my $i=0;$i<3;$i++)
    {
        my @useful_genes = ();
        my @useful_rbhs = ();
        
        if($i==0)
        {
            foreach my $rbh (keys %rbh_hash_1_2)
            {
                if($rbh_hash_1_2{$rbh} ne "none")
                {
                    push(@useful_genes, $rbh);
                    push(@useful_rbhs, $rbh_hash_1_2{$rbh});
                }
            }
            
            %rbh_hash_1_2 = ();
            
            @useful_genes = shuffle(@useful_genes);
            @useful_rbhs = shuffle(@useful_rbhs);
            
            for(my $i=0;$i<scalar(@useful_genes);$i++)
            {
                $rbh_hash_1_2{$useful_genes[$i]} = $useful_rbhs[$i];
            }
        }
        elsif($i==1)
        {
            foreach my $rbh (keys %rbh_hash_1_3)
            {
                if($rbh_hash_1_3{$rbh} ne "none")
                {
                    push(@useful_genes, $rbh);
                    push(@useful_rbhs, $rbh_hash_1_3{$rbh});
                }
            }
            
            %rbh_hash_1_3 = ();
            
            @useful_genes = shuffle(@useful_genes);
            @useful_rbhs = shuffle(@useful_rbhs);
            
            for(my $i=0;$i<scalar(@useful_genes);$i++)
            {
                $rbh_hash_1_3{$useful_genes[$i]} = $useful_rbhs[$i];
            }
        }
        elsif($i==2)
        {
            foreach my $rbh (keys %rbh_hash_1_4)
            {
                if($rbh_hash_1_4{$rbh} ne "none")
                {
                    push(@useful_genes, $rbh);
                    push(@useful_rbhs, $rbh_hash_1_4{$rbh});
                }
            }
            
            %rbh_hash_1_4 = ();
            
            @useful_genes = shuffle(@useful_genes);
            @useful_rbhs = shuffle(@useful_rbhs);
            
            for(my $i=0;$i<scalar(@useful_genes);$i++)
            {
                $rbh_hash_1_4{$useful_genes[$i]} = $useful_rbhs[$i];
            }
        }
    }
}

#<=== GET GENES ===>#
my $local_db = get_sequence_database("ensembl_local");
# my $local_db = get_sequence_database("ensemblgenomes");
# my $local_db = get_sequence_database("ensembl"); # Sequences::Database::Sequence::Ensembl->new


my @species_1_genes = @{$local_db->get_all_accessions($species_1)};
my @species_2_genes = @{$local_db->get_all_accessions($species_2)};
my @species_3_genes = @{$local_db->get_all_accessions($species_3)};
my @species_4_genes = @{$local_db->get_all_accessions($species_4)};

#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes
#If you want to start on a particular gene, set begin to 0
my $start_id = "VIT_07s0005g04410";
my $begin = 0;
foreach my $s1_gene_accession (@species_1_genes)
{
    print "\nGene: $s1_gene_accession";
    if($s1_gene_accession eq $start_id)
    {
        $begin = 1;
    }
    
    if($begin)
    {
        my @curgene = ();
        #Now do the conservation for this gene
        for(my $i=0;$i<3;$i++)
        {
            my $species_two;
            my $db_two;
            my $s2_gene_accession = "undef";
            if($i==0)
            {
                $species_two = $species_2;
                $db_two = $db_2;
                if(defined($rbh_hash_1_2{$s1_gene_accession}) && $rbh_hash_1_2{$s1_gene_accession} ne "none")
                {
                    $s2_gene_accession = $rbh_hash_1_2{$s1_gene_accession};
                }
            }
            elsif($i==1)
            {
                $species_two = $species_3;
                $db_two = $db_3;
                if(defined($rbh_hash_1_3{$s1_gene_accession}) && $rbh_hash_1_3{$s1_gene_accession} ne "none")
                {
                    $s2_gene_accession = $rbh_hash_1_3{$s1_gene_accession};
                }
            }
            elsif($i==2)
            {
                $species_two = $species_4;
                $db_two = $db_4;
                if(defined($rbh_hash_1_4{$s1_gene_accession}) && $rbh_hash_1_4{$s1_gene_accession} ne "none")
                {
                    $s2_gene_accession = $rbh_hash_1_4{$s1_gene_accession};
                }
            }
            
            print "\n$species_two ortholog: $s2_gene_accession";
            
            if($s2_gene_accession ne "undef")
            {
                #Start making sequences
                print "\nConservation script 1/2 call to get_gene_sequence_by_accession, (species 1: $species_1, gene_accession: $s1_gene_accession)\n";
                my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species_1, $s1_gene_accession);
                print "\nConservation script 2/2 call to get_gene_sequence_by_accession, (species 2: $species_two, gene_accession: $s2_gene_accession)\n";
                my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_two, $s2_gene_accession);
                
                #Positive strand error fixing (some genes start at the beginning of chromosomes)
                my $gene_1_start = $sequence_length;
                my $gene_2_start = $sequence_length;
                my $g1_strand = $gene_1_sequence->[0]->{"strand"};
                my $g2_strand = $gene_2_sequence->[0]->{"strand"};
                
                $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
                $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
                
                #Upstream sequence length defaults
                my $species_1_length = $sequence_length;
                my $species_2_length = $sequence_length;
                
                #If we can't physically take that length of sequence, then we need to cut it down
                #Don't need to bother doing this for negative strand really
                if($gene_1_start < $species_1_length)
                {
                    $species_1_length = $gene_1_start;
                }
                if($gene_2_start < $species_2_length)
                {
                    $species_2_length = $gene_2_start;
                }
                
                #Take both of the upstream sequences
                print "\nConservation script 293, call get_sequence_by_location\n";
                my $species_1_sequence = $local_db->get_sequence_by_location
                (
                Sequences::Database::Relative_Location->new
                (
                identifier => $s1_gene_accession,
                offset => - $species_1_length,
                'length' => $species_1_length,
                stop_at_neighbours => 1,
                'anchor' => "5' end"
                )
                );
                
                print "\nConservation script 306, call get_sequence_by_location\n";
                my $species_2_sequence = $local_db->get_sequence_by_location
                (
                Sequences::Database::Relative_Location->new
                (
                identifier => $s2_gene_accession,
                offset => - $species_2_length,
                'length' => $species_2_length,
                stop_at_neighbours => 1,
                'anchor' => "5' end"
                )
                );
                print "\nConservation script 318\n";
                if(length($species_1_sequence->[0]->seq) >= $window_size && length($species_2_sequence->[0]->seq) >= $window_size)
                {
                    #If we're ready to conduct the conservation analysis, then...
                    #Create the sequences
                    my $species_1_final = Datatypes::Sequence::Local->new_from_string
                    (
                    $species_1_sequence->[0]->seq, $species_1,
                    );
                    my $species_2_final = Datatypes::Sequence::Local->new_from_string
                    (
                    $species_2_sequence->[0]->seq, $species_two,
                    );
                    print "\nConservation script 331, create seaweed job.\n";
                    #Set up the job
                    my $job = Jobs::Subtasks::Seaweed_Job->new
                    (
                    sequence_1 => $species_1_final,
                    sequence_2 => $species_2_final,
                    windowsize => $window_size,
                    );
                    print "\nConservation script 339, run seaweed job.\n";
                    #And run the job
                    my $result = $job->run;
                    print "\nConservation script 342, finished seaweed job.\n";
                    my $alignmax = 0;
                    foreach my $cur_val (values %{$result->{"plot"}})
                    {
                        if($cur_val > $alignmax)
                        {
                            $alignmax = $cur_val;
                        }
                    }
                    
                    #Convert alignment results to APPLES format
                    my $alignment_results = convert_gs_seaweed_result_to_ap($result->plot);
                    
                    #<==Bundling==>#
                    my $rdb1 = Ensembl_Database_Parameters->new("dbname" => $db_1,
                    "alias" => $species_1,
                    "location" => "local");
                    
                    my $rdb2 = Ensembl_Database_Parameters->new("dbname" => $db_two,
                    "alias" => $species_two,
                    "location" => "local");
                    
                    #Get sequence information for the genomic intervals
                    my ($g1_region, $g1_coordsys) = $local_db->get_gene_information_by_accession($species_1, $s1_gene_accession);
                    my ($g2_region, $g2_coordsys) = $local_db->get_gene_information_by_accession($species_two, $s2_gene_accession);
                    
                    my $g1_5prime = $species_1_sequence->[0]->{"five_prime_pos"};
                    my $g1_3prime = $species_1_sequence->[0]->{"three_prime_pos"};
                    my $g1_strand = $species_1_sequence->[0]->{"strand"};
                    
                    my $g2_5prime = $species_2_sequence->[0]->{"five_prime_pos"};
                    my $g2_3prime = $species_2_sequence->[0]->{"three_prime_pos"};
                    my $g2_strand = $species_2_sequence->[0]->{"strand"};
                    
                    #Set genomic interval data
                    my $s1_ginterval = Genomic_Interval->new("genome_db" => $rdb1,
                    "region" => $g1_region,
                    "five_prime_pos" => $g1_5prime,
                    "three_prime_pos" => $g1_3prime,
                    "strand" => $g1_strand,
                    "working_sequence" => "ref_sequence",
                    "coord_sys_name" => $g1_coordsys);
                    
                    #Just in case
                    $s1_ginterval->{"gi_sequence"} = $species_1_sequence->[0]->seq;
                    $s1_ginterval->{"gi_sequence_repeatmasked"} = $species_1_sequence->[0]->{"masked_sequence"};
                                        
                    #Set genomic interval data
                    my $s2_ginterval = Genomic_Interval->new("genome_db" => $rdb2,
                    "region" => $g2_region,
                    "five_prime_pos" => $g2_5prime,
                    "three_prime_pos" => $g2_3prime,
                    "strand" => $g2_strand,
                    "working_sequence" => "ref_sequence",
                    "coord_sys_name" => $g2_coordsys);
                    
                    #Just in case
                    $s2_ginterval->{"gi_sequence"} = $species_2_sequence->[0]->seq;
                    $s2_ginterval->{"gi_sequence_repeatmasked"} = $species_2_sequence->[0]->{"masked_sequence"};
                    
                    $curgene[$i]->{"gi"} = $s2_ginterval;
                    $curgene[$i]->{"alignment"} = $alignment_results;
                    $curgene[$i]->{"gene_2_start"} = $gene_2_start;
                    $curgene[$i]->{"gene_2_acc"} = $s2_gene_accession;
                    $curgene[$i]->{"gene_2_strand"} = $g2_strand;
                    $curgene[$i]->{"gene_2_alignmax"} = $alignmax;
                    
                    
                    $curgene[3]->{"gi"} = $s1_ginterval;
                    $curgene[3]->{"gene_1_start"} = $gene_1_start;
                    $curgene[3]->{"gene_1_strand"} = $g1_strand;
                }
            }
        }

        print "\nConservation script 417\n";
        
        my @arr = ();
        my @remo_sets;
        if(defined($curgene[0]->{"gi"}))
        {
            push(@arr, $curgene[0]->{"gi"});
        }
        if(defined($curgene[1]->{"gi"}))
        {
            push(@arr, $curgene[1]->{"gi"});
        }
        if(defined($curgene[2]->{"gi"}))
        {
            push(@arr, $curgene[2]->{"gi"});
        }
        
        if(scalar(@arr) > 0)
        {
            #Make an interval set for comparison
            my $interval_set = Genomic_Interval_Set->new("genomic_interval_set" => \@arr);
            
            #Make an evolutionary tree (will not be used yet)
            my $evtree = Evolutionary_Tree->new("root" => "nasonia vitripennis",
            "sub_trees" => []);
            
            #Make the parameters
            my $ptm = Partial_Threshold_Matrix->new("evolutionary_tree" => $evtree);
            
            my $wpap = Seaweed_Algorithm_Parameters->new("stepwidth" => 1,
            "windowlength" => 50,
            "cutoff_for_uninteresting_alignments" => 57);
            
            my $sbp = Star_Bundler_Parameters->new("overlap_tolerance" => 20,
            "belief_value" => 0.05,
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
            
            if(defined($curgene[0]->{"gene_2_acc"}))
            {
                push(@{$alignment_results}, $curgene[0]->{"alignment"});
            }
            if(defined($curgene[1]->{"gene_2_acc"}))
            {
                push(@{$alignment_results}, $curgene[1]->{"alignment"});
            }
            if(defined($curgene[2]->{"gene_2_acc"}))
            {
                push(@{$alignment_results}, $curgene[2]->{"alignment"});
            }
            
            print "\nScalar: " . scalar(@{$alignment_results});
            print "\nScalar: " . scalar(@arr);

            print "\nConservation script 482\n";
            
            my $bundler = Star_Bundler->new();

            print "\nConservation script 486\n";

            @remo_sets = $bundler->truncated_bundle($curgene[3]->{"gi"}, $interval_set, $parameters, $alignment_results);
            
            print $outfile "--PairStart\n";
            print $outfile $s1_gene_accession . "\n";
            print $outfile $curgene[3]->{"gene_1_start"} . "\n";
            print $outfile $curgene[3]->{"gene_1_strand"} . "\n";
            
            if(exists($curgene[0]->{"gene_2_strand"}))
            {
                print $outfile $curgene[0]->{"gene_2_acc"} . "\n";
                print $outfile $curgene[0]->{"gene_2_start"} . "\n";
                print $outfile $curgene[0]->{"gene_2_strand"} . "\n";
                print $outfile $curgene[0]->{"gene_2_alignmax"} . "\n";
            }
            else
            {
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
            }
            if(exists($curgene[1]->{"gene_2_strand"}))
            {
                print $outfile $curgene[1]->{"gene_2_acc"} . "\n";
                print $outfile $curgene[1]->{"gene_2_start"} . "\n";
                print $outfile $curgene[1]->{"gene_2_strand"} . "\n";
                print $outfile $curgene[1]->{"gene_2_alignmax"} . "\n";
            }
            else
            {
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
            }
            if(defined($curgene[2]->{"gene_2_strand"}))
            {
                print $outfile $curgene[2]->{"gene_2_acc"} . "\n";
                print $outfile $curgene[2]->{"gene_2_start"} . "\n";
                print $outfile $curgene[2]->{"gene_2_strand"} . "\n";
                print $outfile $curgene[2]->{"gene_2_alignmax"} . "\n";
            }
            else
            {
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
                print $outfile "none" . "\n";
            }
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
