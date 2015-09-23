#!/usr/bin/perl

my $out_fn = './random.txt';
open my $outfile, ">$out_fn";
my  @outputArray = ("");
$outputArray[0] = "\n0";
$outputArray[10] = "\n10";
$outputArray[5] = "\n5";

print $outfile @outputArray;
print "Array size: " . scalar @outputArray . "\n";

@outputArray = ("");
$outputArray[100] = "\n100";

print $outfile @outputArray;
print "Array size: " . scalar @outputArray . ", $outputArray[10]\n";



use strict;

use Getopt::Std;

use File::Glob ':glob';
use File::Spec qw(catfile);
use File::Basename;
use File::Temp qw(tempfile);

my %opts;

getopts('u:p:d:x', \%opts);

my $mysql_user = $opts{u} || 'root';
my $mysql_pass = $opts{p} ||'geheim';

my $d = $opts{d};
my $execute = $opts{x};

my @files = <*>;

print @files;

exit;