#!/usr/bin/perl
### RBH script for iPlant. Based on LB's RBH script LB_rbh_search_grape_sorghum.pl ###
### Bo Gao ###

use strict;

# BASEDIR contains path to APPLES/ORANGES base directory
use lib "../webseaweeds/";
BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
	1;
}

use Runtime;
use Data::Dumper;
use Orthology::Search::RBH;
use Links::Output::SVG;
#use Bio::EnsEMBL::Registry;

#<== Parameters ==>#

# Quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: rbhSearch.pl species1 species2\n";
    exit;
}

# Whether there is a '.fa' in the species name does not affect the script.
# Use the following to force .fa file
# if ($ARGV[0] !~ /.fa$/ || $ARGV[1] !~ /.fa$/) {
#     print "\nPlease provide .fa files as arguments.\n";
#     print "Usage: rbhSearch.pl species1.fa species2.fa\n";
#     exit;
# }
if ($ARGV[0] =~ /.fa$/) {
    $ARGV[0] = substr($ARGV[0], 0, -3);
}
if ($ARGV[1] =~ /.fa$/) {
    $ARGV[1] = substr($ARGV[1], 0, -3);
}

# my $species_1 = "vitis vinifera blastp";
# my $species_2 = "medicago truncatula blastp"; #"arabidopsis thaliana blastp";
# my $species_1 = "plantB blastp";
# my $species_2 = "plantA blastp"; #"arabidopsis thaliana blastp";

# Added 'blastp' for default.
# Error if no 'blastp': Name is not unique: plantV.fa (matches : plantV.fa_blastp plantV.fa_blastx)
my $species_1 = $ARGV[0] . ' blastp'; 
my $species_2 = $ARGV[1] . ' blastp';

print "RBH from $ARGV[0] to $ARGV[1].\n";

#my $in_fn = "/home/grannysmith/data/nasonia_bombyx_rbh.txt";
my $out_fn = '../outputs/rbhSearch_result_' . $ARGV[0] . '_' . $ARGV[1] . '.txt';
open my $outfile, ">$out_fn";
my $total_genes = 0;
my $rbh_succ = 0;
my %already_done = ();

### from HERE is optional, if half-way through a batch and need to restart ###

#<== Get already filled-in ones ==>#
#my $in_fn = "/home/grannysmith/data/grape_poplar_rbh.txt";
#open my $infile, "<$in_fn", or die "\nError: $in_fn is missing";
#$_ = <$infile>;
#
#while(<$infile>)
#{
#    chomp;
#    my @split = split("\t");
#    
#     $already_done{$split[0]}->{"gid"} = $split[3];
#     $already_done{$split[0]}->{"rbhid"} = $split[1];
#    $already_done{$split[0]}->{"rbhgid"} = $split[2];
#}
#close $infile;
### to HERE ###

#<== Get genes ==>#
my $db_1 = get_sequence_database($species_1);
my $db_2 = get_sequence_database($species_2);
my @species_1_genes = ();

my @tempgenes = keys %{$db_1->{"sindex"}->{"_DB"}};
#print Dumper (@tempgenes);exit;# debugging for OS statement below
foreach my $s1_accession (@tempgenes)
{
    # print "db_1:\n$db_1\n";
    # print "tempgenes:\n@tempgenes\n";
    # print @tempgenes;
    # if($s1_accession =~ m/VIT/)
    # if($s1_accession =~ m/OS/)
    if($s1_accession !~ /^__/){
        push(@species_1_genes, $s1_accession);
    } else {
        print "filtered accession code: $s1_accession\n";
    }
}

#<== Perform RBH ==>
my $searcher = Orthology::Search::RBH->new(
$species_1,
$species_2
);

foreach my $s1_accession (@species_1_genes)
{
    #Info about genes
    my $source_gid = $db_1->get_sequence($s1_accession)->[0]->desc;
    if($source_gid =~ m/gene\:([^\s]+)/) {
        $source_gid = $1;
    } else {
        $source_gid = "unknown";
    }
    $total_genes++;

    print "\n($total_genes) Calculating RBH for $s1_accession ($source_gid)";

    if(defined($already_done{$s1_accession}))
    {
        my $target_gid = $already_done{$s1_accession}->{"rbhgid"};
        my $target_accession = $already_done{$s1_accession}->{"rbhid"};
        my $source_gid = $already_done{$s1_accession}->{"gid"};
        print " --Loaded! ($target_accession)";
        
        print $outfile "\n$s1_accession\t$target_accession\t$target_gid\t$source_gid";
    }
    else
    {
        my $target_accession = "none";
        my $target_gid = "none";
        
        #Get result
        my $result = $searcher->search($s1_accession, 1);
        my $links = $result->get_links();
        
        foreach my $l (@$links)
        {
            #If we have found an RBH
            if(defined $l->data('RBH'))
            {
                #If it's from this accession
                if($l->source->accession eq $s1_accession)
                {
                    if($target_accession eq "none")
                    {
                        $target_accession = $l->target->accession;
                        $rbh_succ++;
                    }
                    else
                    {
                        die "\nTwo RBHs found for one gene";
                    }
                }
            }
        }
        
        if($target_accession ne "none")
        {
            $target_gid = $db_2->get_sequence($target_accession)->[0]->desc;
            if($target_gid =~ m/gene\:([^\s]+)/) {
                $target_gid = $1;
            } else {
                $target_gid = "unknown";
            }
            print " --Found! ($target_gid)";
        }
        else
        {
            print " --Not found!";
        }
        
        print $outfile "\n$s1_accession\t$target_accession\t$target_gid\t$source_gid";
    }

}
print "\nTotal genes analysed: $total_genes";
print "\nSuccessful RBH hits: $rbh_succ";
print "\n";
exit;
