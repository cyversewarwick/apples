### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

# $Id: APPLES_library.pm 809 2010-12-19 21:48:55Z sascha $

our $APPLES_SVN_REVISION = q$Rev: 809 $;

## should include all APPLES components, update this file when you create new components

use APPLES_Datatypes;
use constant {FALSE => 0,
	      TRUE	=> 1};	
use Alignment_Set;
use Alignment_Utilities;
use BiFa_Raw_Result;
use BiFa_Server_Interface;
use BiSi;
use BiSi_Set;
use Cache;
use Character_Matrix;
use Character_Matrix_Maker;
use ChIP_Seq_Fragment_Set;
use Conservation_Profiles;
use Conservation_Profiles_Maker;
use Evolutionary_Tree;
use Exception;
use General_Utilities;
use Generic_Pattern_Matching_Model;
use Generic_Sequence_Pattern;
use Genome_DB_Utilities;
use Genomic_Interval;
use Genomic_Interval_Set;
use Genomic_Interval_Set_Maker;
use Genomic_Data_Miner;
use Generic_Graph;
use Graph_Clustering;
use Graphics;
use Hypergeometric_Test_Result;
use Job_Handler;
use MEME_wrapper;
use Memory_Cache;
use Molecule;
use Motif_Finder;
use Motif_Finding_Result;
use Motif_Finding_Result_Maker;
use Ortholog_Mapping;
use Ortholog_Mapping_Maker;
use Parameters;
use Parameters_Maker;
use Partial_Threshold_Matrix;
use Partial_Threshold_Matrix_Maker;
use ReCo;
use ReCo_Maker;
use Reg_Loc;
use Reg_Loc_Maker;
use Release_Notes;
use ReMo;
use ReMo_Maker;
use ReMo_Set;
use ReMo_Set_Maker;
use ReRe;
use ReRe_Maker;
use ReRe_Set;
use ReRe_Set_Maker;
use Running_Perseverance;
use Running_Perseverance_Maker;
use Scheduler;
use Star_Bundler;
use Sequence_Model_Learner;
use Statistics;
use Tree_Maker;
use Version_Number;
use WM_Groups;
use Window_Pair_Alignment;
use WM_Scores;
use WM_Utilities;
use XML_Writer;

1; # needed in files that do have a 'use MooseX::Declare' statement
