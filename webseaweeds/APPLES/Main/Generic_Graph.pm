### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Generic_Graph Class ###
use MooseX::Declare;

class Generic_Graph{
	use General_Utilities;
	use APPLES_Datatypes qw (Boolean);
	
	has 'node_labels' => (is => 'ro', isa => 'ArrayRef', required => 0, trigger => \&check_unique_labels);
	has 'node_labels_hash' => (is => 'rw', isa => 'HashRef');
	has 'edge_weights' => (is => 'rw', isa => 'ArrayRef[ArrayRef]', trigger => \&validate_2d_array);
	has 'node_connections' => (is => 'rw', isa => 'ArrayRef[ArrayRef[Boolean]]');
	
	my $GU = General_Utilities->new();
	
	# Method to check all labels are unique and to create hash containing key => value pairs: label => @$self->node_labels index
	sub check_unique_labels(){
		
		my ($self, $arg) = @_;
		
		my $labels = $self->node_labels;
		my $occurrence = 0;

		# Check for two occurrences fo the same label in array
		# This is an inefficient method but should suffice for this problem
		foreach my $x ( @$labels ){
			foreach my $y ( @$labels ){
				if ( $x eq $y ){
					$occurrence++;
				}
			}	
			if ( $occurrence > 1 ){
				die "Node labels provided are not unique\n";
			}
			$occurrence = 0;
		}
		# Create hash where each label has its own key and the value is the index of that
		# label in the node_labels array. This will be used to decrease the time taken to 
		# look up distances between weight matrices. 
		my %label_hash;
		my $index = 0;
		foreach( @{$self->node_labels} ){
			$label_hash{$_} = $index;
			$index++;
		}
		
		$self->node_labels_hash(\%label_hash);
		
	} # check_unique_labels #
	
	sub validate_2d_array(){
		
		my ($self, $arg) = @_;
		
		if( scalar(@{$self->node_labels}) != scalar(@{$self->edge_weights}) ){
			die "There should be as many rows in edge_weights matrix as there are node_labels\n";
		}
		
	} # validate_2d_array #
	
	method retrieve_weight(Str $labela, Str $labelb){
			
		my $index_a = ${$self->node_labels_hash}{$labela};
		my $index_b = ${$self->node_labels_hash}{$labelb};
		my $weight;
		
		if ( ${$self->node_connections}[$index_a][$index_b] ){
			$weight = ${$self->edge_weights}[$index_a][$index_b];
		}
		else{
			$weight = undef;
		}
		
		return $weight;
	} # retrieve_weight #
	
	method print_edge_weights_table(){
			
		foreach my $line ( @{$self->edge_weights}){
				
			
			foreach( @{$line}){
				# To fit everything nicely within one tab, truncate any number to 4 sig. figs.
				my $toprint = substr($_, 0, 4);
				print $toprint."\t";
			}
			print "\n";
		}
	} # print_edge_weights_table #

} # Generic_Graph #

