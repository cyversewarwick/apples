#!/usr/bin/perl

# derived form Nathaniel's conservation_multiple.pl
# my @species = ("nasonia vitripennis", "apis mellifera", "atta cephalotes", "bombyx mori", "acyrthosiphon pisum", "danaus plexippus", "heliconius melpomene", "tribolium castaneum");
# my @dbs = ("nasonia_vitripennis_core_18_71_1", "apis_mellifera_core_18_71_1", "atta_cephalotes_core_18_71_1", "bombyx_mori_core_18_71_1", "acyrthosiphon_pisum_core_18_71_1", "danaus_plexippus_core_18_71_1", "heliconius_melpomene_core_18_71_1", "tribolium_castaneum_core_18_71_3");


use strict;

# BASEDIR contains path to ORANGES base directory
use lib "../../";

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

my @species = ("nasonia vitripennis", "apis mellifera", "atta cephalotes", "bombyx mori", "acyrthosiphon pisum", "danaus plexippus", "heliconius melpomene", "tribolium castaneum");
my @dbs = ("nasonia_vitripennis_core_18_71_1", "apis_mellifera_core_18_71_1", "atta_cephalotes_core_18_71_1", "bombyx_mori_core_18_71_1", "acyrthosiphon_pisum_core_18_71_1", "danaus_plexippus_core_18_71_1", "heliconius_melpomene_core_18_71_1", "tribolium_castaneum_core_18_71_3");
#Rbh file for 1-2, 1-3, 1-4, etc
my @rbhs = ("/home/grannysmith/data/new_rbh/nasonia_apis_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_atta_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_bombyx_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_acrythosiphon_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_danaus_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_heliconiusmelpomene_rbh_new.txt", "/home/grannysmith/data/new_rbh/nasonia_tribolium_rbh_new.txt");

my $sequence_length = 2000;
my $window_size = 50;

my $pseudo_orthologs = 1;

my $outfile_fn = "testing_pseudo.txt";
open my $outfile, ">$outfile_fn";

#<=== LOAD RBHS ===>#
my @rbh_results = ();

for(my $i=0;$i<scalar(@species)-1;$i++)
{
    open my $rbhs_data, "<$rbhs[$i]", or die "\nError: could not open rbh file '$rbhs[$i]'";
    $_ = <$rbhs_data>;
    while(<$rbhs_data>)
    {
        chomp;
        my @split = split(/\t/);
        
        #if($split[2] ne "none")
        #{
            push(@{$rbh_results[$i]->{$split[0]}}, $split[1]);
        #}
    }
}

#<== SET PSEUDO ORTHOLOGS ==>#
if($pseudo_orthologs)
{
    for(my $i=0;$i<scalar(@species);$i++)
    {
        my @useful_genes = ();
        my @useful_rbhs = ();
        my %rbh_numbers;
        
        foreach my $rbh (keys %{$rbh_results[$i]})
        {
            push(@useful_genes, $rbh);
            my $no_rbhs = 0;
            foreach my $ortholog (@{$rbh_results[$i]->{$rbh}})
            {
                push(@useful_rbhs, $ortholog);
                $no_rbhs++;
            }
            $rbh_numbers{$rbh} = $no_rbhs;
        }
        
        $rbh_results[$i] = ();
        
        @useful_genes = shuffle(@useful_genes);
        @useful_rbhs = shuffle(@useful_rbhs);
        
        my $j = 0;
        for(my $k=0;$k<scalar(@useful_genes);$k++)
        {
            while($rbh_numbers{$useful_genes[$k]})
            {
                push(@{$rbh_results[$i]->{$useful_genes[$k]}}, $useful_rbhs[$j]);
                $j++;
                $rbh_numbers{$useful_genes[$k]}--;
            }
        }
    }
}

#<=== GET GENES ===>#
my $local_db = get_sequence_database("ensembl_local");

my @genes = ();
for(my $i=0;$i<scalar(@species);$i++)
{
    $genes[$i] = $local_db->get_all_accessions($species[$i]);
}

#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes
my $start_id = "NV10220";
my $begin = 1;
foreach my $s1_gene_accession (@{$genes[0]})
{
    #<-- If we're ready to begin, start with a gene -->#
    print "\nGene: $s1_gene_accession";
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
        
        #<-- Firstly, set up S1 gene -->#
        my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species[0], $s1_gene_accession);
        
        #<-- Positive strand error fixing -->#
        my $gene_1_start = $sequence_length;
        my $g1_strand = $gene_1_sequence->[0]->{"strand"};
        $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
        
        #<-- Upstream sequence length default -->#
        my $species_1_length = $sequence_length;
        
        #If we can't physically take that length of sequence, then we need to cut it down
        #Don't need to bother doing this for negative strand really
        if($gene_1_start < $species_1_length)
        {
            $species_1_length = $gene_1_start;
        }
        
        #Take the upstream sequence
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
        
        #If the length is above the window size
        if(length($species_1_sequence->[0]->seq) >= $window_size)
        {
            print "\n$species[0] sequence ok!";
            #Remove IUPAC codes
            my $sequence_one = $species_1_sequence->[0]->seq;
            $sequence_one =~ s/[^(A|T|C|G)]/N/g;
            
            #Create the new sequences
            my $species_1_final = Datatypes::Sequence::Local->new_from_string
            (
                $sequence_one, $species[0],
            );
            
            #Set up the stuff we need for the bundler run with gene 1
            my $rdb1 = Ensembl_Database_Parameters->new("dbname" => $dbs[0],
                "alias" => $species[0],
                "location" => "local");
            
            #Get sequence information
            my ($g1_region, $g1_coordsys) = $local_db->get_gene_information_by_accession($species[0], $s1_gene_accession);
            my $g1_5prime = $species_1_sequence->[0]->{"five_prime_pos"};
            my $g1_3prime = $species_1_sequence->[0]->{"three_prime_pos"};
            my $g1_strand = $species_1_sequence->[0]->{"strand"};
            
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
            
            #Set s1 gene comparison information
            $curgene_comparisons[0]->{"species"} = $species[0];
            $curgene_comparisons[0]->{"db"} = $dbs[0];
            $curgene_comparisons[0]->{"gi"} = $s1_ginterval;
            $curgene_comparisons[0]->{"gene_1_start"} = $gene_1_start;
            $curgene_comparisons[0]->{"gene_1_strand"} = $g1_strand;
            
            #Now that species 1 is done, we can do the same for all of the other species
            #Now do the conservation for this gene
            #For each species
            for(my $i=1;$i<scalar(@species);$i++)
            {
                #Which species is it?
                my $species_two = $species[$i];
                my $db_two = $dbs[$i];
                
                print "\nSpecies: " . $species[$i];
                
                #Get orthologues for this gene
                my @s2_genes = ();
                if(defined($rbh_results[$i-1]->{$s1_gene_accession}))
                {
                    @s2_genes = @{$rbh_results[$i-1]->{$s1_gene_accession}};
                }
                
                #Do we have any orthologues for this gene?
                if(scalar(@s2_genes) > 0)
                {
                    #Go through all s2 genes
                    foreach my $s2_gene_accession (@s2_genes)
                    {
                        print "\n\t$species_two ortholog: $s2_gene_accession";
                        
                        #Start making sequences
                        my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_two, $s2_gene_accession);
                        
                        #Positive strand error fixing (some genes start at the beginning of chromosomes)
                        my $gene_2_start = $sequence_length;
                        my $g2_strand = $gene_2_sequence->[0]->{"strand"};
                        $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
                        
                        #Upstream sequence length defaults
                        my $species_2_length = $sequence_length;
                        
                        #If we can't physically take that length of sequence, then we need to cut it down
                        #Don't need to bother doing this for negative strand really
                        if($gene_2_start < $species_2_length)
                        {
                            $species_2_length = $gene_2_start;
                        }
                        
                        #Take the upstream sequences
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
                        
                        if(length($species_2_sequence->[0]->seq) >= $window_size)
                        {
                            #Remove IUPAC codes
                            my $sequence_two = $species_2_sequence->[0]->seq;
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
                            my $result = $job->run;
                            
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
                            
                            my $rdb2 = Ensembl_Database_Parameters->new("dbname" => $db_two,
                            "alias" => $species_two,
                            "location" => "local");
                            
                            #Get sequence information for the genomic intervals
                            my ($g2_region, $g2_coordsys) = $local_db->get_gene_information_by_accession($species_two, $s2_gene_accession);
                            
                            my $g2_5prime = $species_2_sequence->[0]->{"five_prime_pos"};
                            my $g2_3prime = $species_2_sequence->[0]->{"three_prime_pos"};
                            my $g2_strand = $species_2_sequence->[0]->{"strand"};
                            
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
                            
                            $no_comparisons++;
                            
                            #Set s2 gene
                            $curgene_comparisons[$no_comparisons]->{"species"} = $species_two;
                            $curgene_comparisons[$no_comparisons]->{"db"} = $db_two;
                            $curgene_comparisons[$no_comparisons]->{"gi"} = $s2_ginterval;
                            $curgene_comparisons[$no_comparisons]->{"alignment"} = $alignment_results;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_start"} = $gene_2_start;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_acc"} = $s2_gene_accession;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_strand"} = $g2_strand;
                            $curgene_comparisons[$no_comparisons]->{"gene_2_alignmax"} = $alignmax;
                        }
                    }
                }
            }
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
            my $evtree = Evolutionary_Tree->new("root" => "nasonia vitripennis",
            "sub_trees" => []);
            
            #Make the parameters
            my $ptm = Partial_Threshold_Matrix->new("evolutionary_tree" => $evtree);
            
            my $wpap = Seaweed_Algorithm_Parameters->new("stepwidth" => 1,
            "windowlength" => $window_size,
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

            for(my $i=1;$i<scalar(@curgene_comparisons);$i++)
            {
                my $cur_g = $curgene_comparisons[$i];

                push(@{$alignment_results}, $cur_g->{"alignment"});
            }

            my $bundler = Star_Bundler->new();
            my @remo_sets = $bundler->truncated_bundle($curgene_comparisons[0]->{"gi"}, $interval_set, $parameters, $alignment_results);
            
            print $outfile "--PairStart\n";
            print $outfile $s1_gene_accession . "\n";
            print $outfile $curgene_comparisons[0]->{"gene_1_start"} . "\n";
            print $outfile $curgene_comparisons[0]->{"gene_1_strand"} . "\n";
            
            for(my $i=1;$i<scalar(@curgene_comparisons);$i++)
            {
                my $cur_g = $curgene_comparisons[$i];
                
                print $outfile $cur_g->{"species"} . "\n";
                print $outfile $cur_g->{"gene_2_acc"} . "\n";
                print $outfile $cur_g->{"gene_2_start"} . "\n";
                print $outfile $cur_g->{"gene_2_strand"} . "\n";
                print $outfile $cur_g->{"gene_2_alignmax"} . "\n";
            }
            
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
