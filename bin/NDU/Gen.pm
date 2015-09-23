package NDU::Gen;
#Version 1.0.0

use lib "../../";

BEGIN {
	use Configuration::AppleSeeds;
	Configuration::AppleSeeds::load_APPLES();
	1;
}

use Runtime;
use strict;
use Exporter;
use Data::Dumper;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
use ReMo_Data;

our $VERSION = 1.00;
our @all_functions = qw(&randomise_sequence &convert_gs_seaweed_result_to_ap);

our @ISA = qw(Exporter);
#Things to export, things that are okay to export, shorthand tags
our @EXPORT      = qw();
our @EXPORT_OK   = @all_functions;
our %EXPORT_TAGS = ( ALL => [@all_functions],);

sub randomise_sequence
{
    my ($sequence) = @_;
    
    my @seq = split("", $sequence);
    
    my $new_seq = "";
    for(my $i=0;$i<length($sequence);$i++)
    {
        my $char_pos = int(rand(scalar(@seq)));
        my $rand_char = splice(@seq, $char_pos, 1);
        $new_seq .= $rand_char;
    }
    return $new_seq;
}

sub convert_gs_seaweed_result_to_ap
{
    my ($gs_plot) = @_;
    
    my @ap_result = ();
    foreach my $window (keys %{$gs_plot})
    {
        my $wp_result = WindowPair->new;
        
        my @window_pair_split = split(/_/, $window);
        
        $wp_result->{"offset1"} = $window_pair_split[0];
        $wp_result->{"offset2"} = $window_pair_split[1];
        
        $wp_result->{"score"} = $gs_plot->{$window};
        
        push(@ap_result, $wp_result);
    }
    
    return \@ap_result;
}