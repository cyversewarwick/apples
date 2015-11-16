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
#Two species that you're going to be comparing
my $species_1 = "vitis vinifera";#"oryza sativa";
my $species_2 = "arabidopsis thaliana";#"arabidopsis thaliana";

my $db_1 = "vitis_vinifera_core_29_82_3";
my $db_2 = "arabidopsis_thaliana_core_29_82_10";

#How much upstream sequence to take
my $sequence_length = 2000;
#Window size for seaweeds algorithm
my $window_size = 60;

my $pseudo_orthologs = 1; # 1=TRUE

my $outfile_fn = "../output/conservation_result_wB_two_species_plantV_plantA_short.txt";
# my $outfile_fn = "../output/conservation_result_two_species_plantV_plantA_long.txt";
open my $outfile, ">$outfile_fn";

#<=== LOAD RBHS ===>#
my %rbhs = ();

my $rbh_file = "../output/rbhSearchForked_result_plantV_plantA.txt"; # short version
# my $rbh_file = "../output/rbhSearchForked_result_Vitis_vinifera_Arabidopsis_thaliana.txt"; # long version
open my $rbhs_data, "<$rbh_file", or die "\nError: Could not open rbh file";
$_ = <$rbhs_data>;
while(<$rbhs_data>)
{
    chomp;
    my @split = split(/\t/);
    # push(@{$rbhs{$split[0]}}, $split[1]);
    $rbhs{$split[3]} = $split[2];
}
close $rbhs_data;
# print Dumper (%rbhs);
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

#<=== GET GENES ===>#
my $local_db = get_sequence_database("ensembl_local");

my @species_1_genes = @{$local_db->get_all_accessions($species_1)};
my @species_2_genes = @{$local_db->get_all_accessions($species_2)};

#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes

my $start_id = "AT3G01850";
my $begin = 1; # 0 to begin with the gene specified in $start_id, 1 otherwise.

# my $total = 0;
# my $count = 0;

foreach my $s1_gene_accession (@species_1_genes)
{
    # $total++;
    
    # print $total . "\n";

    if($s1_gene_accession eq $start_id)
    {
        $begin = 1;
    }
    
    if($begin)
    {
        # print $rbhs{$s1_gene_accession} . "\n";
        # foreach my $s2_gene_accession (@{$rbhs{$s1_gene_accession}})
        if(defined($rbhs{$s1_gene_accession}) && $rbhs{$s1_gene_accession} ne "none")
        {
            #Get the RBH in species 2
            my $s2_gene_accession = $rbhs{$s1_gene_accession};
            
            #Now we have to check if we can take the sequence we want to take
            #So, get the sequences
            my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species_1, $s1_gene_accession);
            my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_2, $s2_gene_accession);
            
            #        print Dumper($gene_2_sequence);
            #Where do they start on the chromosome?
            #Positive strand error fixing (some genes start at the beginning of chromosomes)
            my $gene_1_start = $sequence_length;
            my $gene_2_start = $sequence_length;
            my $g1_real_start;
            my $g2_real_start;
            my $g1_strand = $gene_1_sequence->[0]->{"strand"};
            my $g2_strand = $gene_2_sequence->[0]->{"strand"};
            
            $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
            $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
            
            $g1_real_start = $gene_1_sequence->[0]->{"five_prime_pos"};
            $g2_real_start = $gene_2_sequence->[0]->{"five_prime_pos"};
            
            #Upstream sequence length defaults
            my $species_1_length = $sequence_length;
            my $species_2_length = $sequence_length;
            
            #If we can't physically take that length of sequence, then we need to cut it down
            if($gene_1_start < $species_1_length)
            {
                $species_1_length = $gene_1_start;
            }
            if($gene_2_start < $species_2_length)
            {
                $species_2_length = $gene_2_start;
            }
            
            #Take both of the upstream sequences
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
            
            #Both sequences must be above or equal to the window size
            if(length($species_1_sequence->[0]->seq) >= $window_size && length($species_2_sequence->[0]->seq) >= $window_size)
            {
                #Remove IUPAC codes
                # my $sequence_one = $species_1_sequence->[0]->seq;
                # $sequence_one =~ s/[^(A|T|C|G)]/N/g;
                # my $sequence_two = $species_2_sequence->[0]->seq;
                # $sequence_two =~ s/[^(A|T|C|G)]/N/g;

                #If we're ready to conduct the conservation analysis, then...
                #Create the sequences
                my $species_1_final = Datatypes::Sequence::Local->new_from_string
                (
                    $species_1_sequence->[0]->seq, $species_1,
                );
                my $species_2_final = Datatypes::Sequence::Local->new_from_string
                (
                    $species_2_sequence->[0]->seq, $species_2,
                );
                
                
                #Set up the job
                my $job = Jobs::Subtasks::Seaweed_Job->new
                (
                sequence_1 => $species_1_final,
                sequence_2 => $species_2_final,
                windowsize => $window_size,
                );

                my @db_handles;
                # my $ensembl_registry = 'Bio::EnsEMBL::Registry';
                my $ensembl_registry = $local_db->{"registry"};
                # print Dumper($local_db);
                
                my @db_adaptors = @{ $ensembl_registry->get_all_DBAdaptors() };
                foreach my $db_adaptor (@db_adaptors) {
                    # print "db_adaptor: " . Dumper($db_adaptor);
                    my $db_connection = $db_adaptor->dbc();
                    # push @db_handles, $db_connection->db_handle();
                    $job->add_handle($db_connection->db_handle());

                    # push @db_handles, $db_connection->db_handle();
                    # print "script: $db_connection->db_handle()". Dumper($db_connection->db_handle()) . "\n";
                }
                # print "db_handles:" . Dumper (@db_handles);

                
                #And run the job
                my $result = $job->run;

                print "\nConservation script two species, removing temporary files.\n";
                unlink glob "'./tempfiles/*'";
                print "\nConservation script two species, seaweed tempfiles removed.\n";
                
                #What is the max alignment value for this run of seaweeds?
                my $alignmax = 0;
                foreach my $cur_val (values %{$result->{"plot"}})
                {
                    if($cur_val > $alignmax)
                    {
                        $alignmax = $cur_val;
                    }
                }
                
                my $alignment_results;

                #Convert alignment results to APPLES format
                push(@{$alignment_results}, convert_gs_seaweed_result_to_ap($result->plot));
                
                #<==Bundling==>#
                my $rdb1 = Ensembl_Database_Parameters->new("dbname" => $db_1,
                "alias" => $species_1,
                "location" => "local");
                
                my $rdb2 = Ensembl_Database_Parameters->new("dbname" => $db_2,
                "alias" => $species_2,
                "location" => "local");
                #Get sequence information for the genomic intervals
                my ($g1_region, $g1_coordsys) = $local_db->get_gene_information_by_accession($species_1, $s1_gene_accession);
                my ($g2_region, $g2_coordsys) = $local_db->get_gene_information_by_accession($species_2, $s2_gene_accession);
                
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
                
                #Make an interval set for comparison
                my @arr = ($s2_ginterval);
                my $interval_set = Genomic_Interval_Set->new("genomic_interval_set" => \@arr);
                
                #Make an evolutionary tree (will not be used yet)
                my $evtree1 = Evolutionary_Tree->new("root" => "oryza sativa",
                "sub_trees" => []);
                my $evtree2 = Evolutionary_Tree->new("root" => "musa acuminata",
                "sub_trees" => []);
                my $evtree3 = Evolutionary_Tree->new("root" => "",
                "sub_trees" => [$evtree1, $evtree2]);
                
                #Make the parameters
                my $ptm = Partial_Threshold_Matrix->new("evolutionary_tree" => $evtree3);
                
                my $wpap = Seaweed_Algorithm_Parameters->new("stepwidth" => 1,
                "windowlength" => 60,
                "cutoff_for_uninteresting_alignments" => 57);
                my $sbp = Star_Bundler_Parameters->new("overlap_tolerance" => 20,
                "belief_value" => 0.05,
                "partial_threshold_matrix" => $ptm); # LB - windowlength was set to 50 here, changed it to 60
                
                my $sp = Sequence_Parameters->new("region" => "upstream",
                "min_length_to_return" => 50,
                "max_length_to_search" => 10000,
                );
                
                my $parameters = ReMo_Set_Phylogenetic_Constructor_Parameters->new("window_pair_algorithm_parameters" => $wpap,
                                                                                    "star_bundler_parameters" => $sbp,
                                                                                    "sequence_parameters" => $sp,
                                                                                    "sequence_databases_to_use_for_homologs" => []);
                
                my $bundler = Star_Bundler->new();
                my @remo_sets = $bundler->truncated_bundle($s1_ginterval, $interval_set, $parameters, $alignment_results);
                
                #Print the appropriate output
                print $outfile "--PairStart\n";
                print $outfile $s1_gene_accession . "\n";
                print $outfile $s2_gene_accession . "\n";
                print $outfile $g1_strand . "\n";
                print $outfile $g2_strand . "\n";
                print $outfile $g1_real_start . "\n";
                print $outfile $g2_real_start . "\n";
                print $outfile $species_1_sequence->[0]->seq . "\n";
                print $outfile $species_2_sequence->[0]->seq . "\n";
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
                                
                print "\n$s1_gene_accession -> $s2_gene_accession ($alignmax)";
                
                foreach my $remo_set (@remo_sets)
                {
                    print "\n\t" . $remo_set->{"remo_set"}->[0]->{"belief_score"} . " - " . $remo_set->{"remo_set"}->[1]->{"belief_score"};
                }
                
            }
            #print "\nGen:" . substr($gene_1_sequence->[0]->seq, 0, 10);
            #print "\n" . Dumper($species_1_sequence->[0]->seq);
            #exit;
        }
    }
}