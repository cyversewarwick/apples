#!/usr/bin/perl

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

use Runtime;

use JSON;
use Data::Dumper;

use Sequences::Database::Relative_Location;
use Sequences::Database::Sequence::Ensembl;
use Sequences::Database::Sequence::Genbank;

use Datatypes::Sequence::Local;
use Jobs::Subtasks::Seaweed_Job;

use Serialization::Serializable;

use List::Util qw(shuffle);

#<=== SET PARAMETERS ===>#
#Two species that you're going to be comparing
my $species_1 = "vitis vinifera";#"oryza sativa";
my $species_2 = "arabidopsis thaliana";#"arabidopsis thaliana";

#How much upstream sequence to take
my $sequence_length = 2000;
#Window size for seaweeds algorithm
my $window_size = 60;

my $pseudo_orthologs = 0; # 1=TRUE

my $outfile_fn = "../output/conservation_result_two_species_plantV_plantA_long.txt";
open my $outfile, ">$outfile_fn";

#<=== LOAD RBHS ===>#
my %rbhs = ();

# my $rbh_file = "../output/rbhSearchForked_result_plantV_plantA.txt"; # short version
my $rbh_file = "../output/rbhSearchForked_result_Vitis_vinifera_Arabidopsis_thaliana.txt"; # long version
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
print Dumper (%rbhs);
#<== SET PSEUDO ORTHOLOGS ==>#
if($pseudo_orthologs)
{
    my @useful_genes = ();
    my @useful_rbhs = ();
    my %rbh_numbers;
    
    foreach my $rbh (keys %rbhs)
    {
        push(@useful_genes, $rbh);
        my $no_rbhs = 0;
        foreach my $ortholog (@{$rbhs{$rbh}})
        {
            push(@useful_rbhs, $ortholog);
            $no_rbhs++;
        }
        $rbh_numbers{$rbh} = $no_rbhs;
    }
    
    %rbhs = ();
    
    @useful_genes = shuffle(@useful_genes);
    @useful_rbhs = shuffle(@useful_rbhs);
    
    my $j = 0;
    for(my $i=0;$i<scalar(@useful_genes);$i++)
    {
        while($rbh_numbers{$useful_genes[$i]})
        {
            push(@{$rbhs{$useful_genes[$i]}}, $useful_rbhs[$j]);
            $j++;
            $rbh_numbers{$useful_genes[$i]}--;
        }
    }
}

#<=== GET GENES ===>#
my $local_db = get_sequence_database("ensembl_local");

my @species_1_genes = @{$local_db->get_all_accessions($species_1)};
my @species_2_genes = @{$local_db->get_all_accessions($species_2)};

#<=== BEGIN CONSERVATION SEARCH ===>#
#Go through all S1 genes

my $start_id = "AT3G01850";
my $begin = 1; # 0 to begin with the gene specified in $start_id, 1 otherwise.

my $total = 0;
my $count = 0;

foreach my $s1_gene_accession (@species_1_genes)
{
    $total++;
    
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
            $count++;
            #Get the RBH in species 2
            my $s2_gene_accession = $rbhs{$s1_gene_accession};
            
            #Now we have to check if we can take the sequence we want to take
            #So, get the sequences
            my $gene_1_sequence = $local_db->get_gene_sequence_by_accession($species_1, $s1_gene_accession);

            my $gene_2_sequence = $local_db->get_gene_sequence_by_accession($species_2, $s2_gene_accession);
            
            #        print Dumper($gene_2_sequence);
            #Where do they start on the chromosome?
            my $gene_1_start = $gene_1_sequence->[0]->{"five_prime_pos"};
            my $gene_2_start = $gene_2_sequence->[0]->{"five_prime_pos"};
            
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
                my $sequence_one = $species_1_sequence->[0]->seq;
                $sequence_one =~ s/[^(A|T|C|G)]/N/g;
                my $sequence_two = $species_2_sequence->[0]->seq;
                $sequence_two =~ s/[^(A|T|C|G)]/N/g;

                #If we're ready to conduct the conservation analysis, then...
                #Create the sequences
                my $species_1_final = Datatypes::Sequence::Local->new_from_string
                (
                $sequence_one, $species_1,
                );
                my $species_2_final = Datatypes::Sequence::Local->new_from_string
                (
                $sequence_two, $species_2,
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
                
                #Print the appropriate output
                print $outfile "--PairStart\n";
                print $outfile $s1_gene_accession . "\n";
                print $outfile $s2_gene_accession . "\n";
                print $outfile $species_1_sequence->[0]->seq . "\n";
                print $outfile $species_2_sequence->[0]->seq . "\n";
                print $outfile $alignmax . "\n";
                print $outfile "--PairEnd\n";
                
                print "\n$s1_gene_accession -> $s2_gene_accession ($alignmax)";
                
            }
            #print "\nGen:" . substr($gene_1_sequence->[0]->seq, 0, 10);
            #print "\n" . Dumper($species_1_sequence->[0]->seq);
            #exit;
        }
    }
}

print "\nTotal total: $total vs $count";
exit;
