#!/usr/bin/perl
### RBH script for iPlant. Based on LB's RBH script LB_rbh_search_grape_sorghum.pl ###
### Bo Gao ###

use strict;
use Parallel::ForkManager;
use IPC::Shareable;

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
my $out_fn = '../outputs/rbhSearchForked_result_' . $ARGV[0] . '_' . $ARGV[1] . '.txt';
open my $outfile, ">$out_fn";
my $out_fn2 = $out_fn . '2';
open my $outfile2, ">$out_fn2";
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

## Fork parameters and initilalisation
my $max_procs = 6;
my $pm = Parallel::ForkManager->new($max_procs, '../tempdir/');
my $rbh_succ = 0;
my $rbh_succ_per_proc = 0;

## Shared variable
my %shareable_options = (
     create    => 'yes',
     exclusive => 0,
     mode      => 0644,
     destroy   => 'yes',
 );
my $total_found_shared = 0;
tie $total_found_shared, 'IPC::Shareable', 'data', {%shareable_options};


# my $outputTotal = '';
my @outputArray = ("");

$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
        
        print "\n>>>>Child no.$ident(PID: $pid) found $exit_code.";
        $rbh_succ = $rbh_succ + $exit_code;

        if(defined($data_structure_reference)){
            # $outputTotal = $outputTotal . ${$data_structure_reference};
            $outputArray[$ident] = ${$data_structure_reference};
        }

    }
);

$pm->run_on_start( sub {
    my ($pid, $ident) = @_;
    print "\n** Child process no.$ident started, pid: $pid.";
});

my $childCount = 0;
my $genes_count = scalar @species_1_genes;

foreach my $s1_accession (@species_1_genes)
{

    # if($childCount % 100 == 0) {
    #     print $outfile @outputArray;
    #     @outputArray = ("");
    # }

    if(scalar @outputArray % 100 == 0) {
        print $outfile @outputArray;
        @outputArray = ("");

        # print "\n--progress:$childCount/$genes_count (" . int($childCount*100/$genes_count) . "%)";
    }

    #Info about genes
    my $source_gid = $db_1->get_sequence($s1_accession)->[0]->desc;
    if($source_gid =~ m/gene\:([^\s]+)/) {
        # my @source_gid_temp = $source_gid =~ m/gene\:([^\s]+)/;
        $source_gid = $1;
        # $source_gid = $source_gid_temp[0];
        # print "\nSource gene: " . $source_gid . ", " . $source_gid_temp[0];
    } else {
        $source_gid = "unknown";
    }
    $total_genes++;

    print "\n($total_genes) Calculating RBH for $s1_accession ($source_gid)";


    ## Forking starts here

    $childCount++;

    my $pid = $pm->start($childCount) and next;

    $rbh_succ_per_proc = 0;

    my $outStatement = "";

# my $searcher_proc = Orthology::Search::RBH->new(
# $species_1,
# $species_2
# );

# my $db_1 = get_sequence_database($species_1);
# my $db_2 = get_sequence_database($species_2);

    # if(defined($already_done{$s1_accession}))
    # {
    #     my $target_gid = $already_done{$s1_accession}->{"rbhgid"};
    #     my $target_accession = $already_done{$s1_accession}->{"rbhid"};
    #     my $source_gid = $already_done{$s1_accession}->{"gid"};
    #     print " --Loaded! ($target_accession)";
        
    #     print $outfile "\n$s1_accession\t$target_accession\t$target_gid\t$source_gid";
    # }
    # else
    # {
        my $from_accession = $s1_accession;
        my $from_gid = $source_gid;

        my $target_accession = "none";
        my $target_gid = "none";
        
        #Get result
        
        my $searcher_proc = Orthology::Search::RBH->new(
            $species_1,
            $species_2
        );
        
        my $result = $searcher_proc->search($s1_accession, 1);
        # my $result = $searcher->search($s1_accession, 1);

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
                        # $rbh_succ++;
                        $rbh_succ_per_proc++;

                        (tied $total_found_shared)->shlock;
                        $total_found_shared++;
                        (tied $total_found_shared)->shunlock;
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
        
        # $outStatement = "\n$s1_accession\t$target_accession\t$target_gid\t$source_gid";
        $outStatement = "\n$from_accession\t$target_accession\t$target_gid\t$from_gid";
        
        # print $outfile "\n$s1_accession\t$target_accession\t$target_gid\t$source_gid";
        print $outfile2 "\n$from_accession\t$target_accession\t$target_gid\t$from_gid";

    # }

    $pm->finish($rbh_succ_per_proc, \$outStatement);

}
$pm->wait_all_children;

print $outfile @outputArray;

print "\nTotal genes analysed: $total_genes";
print "\nSuccessful RBH hits: $rbh_succ, $total_found_shared";
print "\n";

# print @outputArray;
exit;
