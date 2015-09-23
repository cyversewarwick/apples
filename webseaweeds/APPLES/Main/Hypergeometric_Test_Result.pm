### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Hypergeometric_Test_Result class ###
use MooseX::Declare;

class Hypergeometric_Test_Result {

    use Generic_Sequence_Pattern;
    use General_Utilities;
	
    my $GU = General_Utilities->new();

    use Moose;
    has 'cluster_sizes' => (is => 'ro', isa => 'HashRef', required => 1);
      # keys are cluster-IDs, content is the size
    has 'cluster_hits' => (is => 'ro', isa => 'HashRef', required => 1);
      # keys are cluster-IDs, content is an array of records containing
      # three entries: P_VALUE, PSSM_ID, PSSM_NAME
    has 'weight_matrices' => (is => 'ro', isa => 'HashRef', required => 1);
      # Generic_Weight_Matrix objects for all PSSM-IDs given in cluster_hits
    
    method render_html(Str $output_directory, $html_filename) {
      # will create html-file and subdirectory for motif logos

	my $full_filename = $output_directory.'/'.$html_filename;
	my $logo_directory = $output_directory.'/logos/';
	open HTML, ">$full_filename" or die 'System error!';
	print HTML "<html>\n";
	print HTML "<body>\n";
	mkdir($logo_directory) or die 'System error!';
	my @one_set_of_keys = keys %{$self->cluster_sizes};
	my @second_set_of_keys = keys %{$self->cluster_hits};
	my @all_keys = (@one_set_of_keys,@second_set_of_keys);
	my @clusters = $GU->remove_duplicates_from_list(\@all_keys);
	foreach my $cluster_ID (@clusters) {
	    my $cluster_size = ${$self->cluster_sizes}{$cluster_ID};
	    print HTML "Cluster: ".$cluster_ID.", size of cluster: ".$cluster_size."\n";
	    my @cluster_hits = @{${$self->cluster_hits}{$cluster_ID}};
	    foreach my $cluster_hit (@cluster_hits) {
		my $pssm_id = $cluster_hit->{PSSM_ID};
		my $p_value = $cluster_hit->{P_VALUE};
		my $pssm_name = $cluster_hit->{PSSM_NAME};
		my $gwm = ${$self->weight_matrices}{$pssm_id};
		print HTML "<p>\n";
		print HTML $p_value."\t".$pssm_name."\t".$pssm_id;
		my @icm = $gwm->generate_logo($logo_directory);
		print HTML "<img src= \"logos/$pssm_id.png\" width=\"144\" height=\"50\" />\n";
		print HTML "</p>";
	    }
	    print HTML "<hr>";
	}
	print HTML "</body>";
	print HTML "</html>";
	close HTML;
	system("rm $logo_directory/*.txt");
    } # render_html #
    
} # Hypergeometric_Test_Result #
