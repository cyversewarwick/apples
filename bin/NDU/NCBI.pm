package NDU::NCBI;
#Version 1.0.0

use strict;
use Exporter;
use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);

our $VERSION = 1.00;
our @all_functions = qw(&send_web_blastx_query &retrieve_web_query &batch_retrieve_web_queries &batch_submit_web_queries &parse_top_e_value_from_blastx);

our @ISA = qw(Exporter);
#Things to export, things that are okay to export, shorthand tags
our @EXPORT      = qw();
our @EXPORT_OK   = @all_functions;
our %EXPORT_TAGS = ( ALL => [@all_functions],);

sub send_web_blastx_query
{
    my ($query) = @_;
    
    #User agent for sending HTTP requests
    my $ua = LWP::UserAgent->new;
    my $program = "blastx";
    
    #Which database?
    my $database = "nr";
    
    #uri_escape query
    my $encoded_query = uri_escape($query);
    
    #Build the request
    my $args = "CMD=Put&PROGRAM=$program&DATABASE=$database&QUERY=" . $encoded_query;
    
    #Make the request
    my $req = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
    $req->content_type('application/x-www-form-urlencoded');
    $req->content($args);
    
    # get the response
    my $response = $ua->request($req);
    
    # parse out the request id
    $response->content =~ /^    RID = (.*$)/m;
    my $rid=$1;
    
    return $rid;
}

sub retrieve_web_query
{
    my ($rid) = @_;
 
    #User agent for sending HTTP requests
    my $ua = LWP::UserAgent->new;
    
    my $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
    my $response = $ua->request($req);
    
    if ($response->content =~ /\s+Status=WAITING/m)
    {
        return (1, "");
    }
    
    if ($response->content =~ /\s+Status=FAILED/m)
    {
        return (2, "");
    }
    
    if ($response->content =~ /\s+Status=UNKNOWN/m)
    {
        return (3, "");
    }
    
    if ($response->content =~ /\s+Status=READY/m)
    {
        if ($response->content =~ /\s+ThereAreHits=yes/m)
        {
            $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=$rid";
            $response = $ua->request($req);

            return (0, $response->content);
        }
        else
        {
            return (4, "");
        }
    }

    return (5, "");
}

sub batch_submit_web_queries
{
    my ($queries_ref) = @_;
    
    my %query_output = ();

    my @queries = @{$queries_ref};
    
    foreach my $query (@queries)
    {
        my $id = send_web_blastx_query($query);
        my @query_split = split("\n", $query);
        print "\nSubmitting $query_split[0]";

        $query_output{$query_split[0]} = $id;
        
        sleep 3;
    }
    
    return \%query_output;
}


sub batch_retrieve_web_queries
{
    my ($queries_ref) = @_;
    
    my %query_output = ();
    
    my @queries = @{$queries_ref};
    
    foreach my $query (@queries)
    {
        print "\nRetrieving $query...";
        my ($code, $info) = retrieve_web_query($query);
        
        $query_output{$query}->{"status"} = $code;
        
        if($code == 0)
        {
            my ($code2, $top_eval) = parse_top_e_value_from_blastx($info);
            if($code2 == 1)
            {
                print "\nError: $query unknown blastx error";
            }
            $query_output{$query}->{"top_e_val"} = $top_eval;
        }
        sleep 3;
    }
    
    return \%query_output;
}

sub parse_top_e_value_from_blastx
{
    my ($blast_output) = @_;
    
    my @lines = split(/\n/, $blast_output);
    
    for (my $i=0;$i<scalar(@lines);$i++)
    {
        if(substr($lines[$i], 0, 13) eq "Sequences pro")
        {
            my $top_val;
            $top_val = substr($lines[$i+2], 75);
            return (0, $top_val);
        }
    }
    
    return (1, "NULL");
    
}
1;