### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

# add any non-APPLES-dependency that may be useful in user scripts
use Cwd;
use Data::Dumper;
use diagnostics;
use feature qw (state);
use File::Temp qw (tempdir);
use File::Spec;
use MooseX::Declare;
use SOAP::Lite;
use Storable qw (nstore retrieve);

1; # needed in files that do have a 'use MooseX::Declare' statement
