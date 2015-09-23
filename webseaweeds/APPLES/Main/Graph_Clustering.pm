### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Graph_Clustering Class ###
use MooseX::Declare;

class Graph_Clustering{
	
	use General_Utilities;
	use Generic_Graph;
	use Data::Dumper;
	
	has 'clusters' => (is => 'rw', isa => 'ArrayRef[ArrayRef]');
	
	my $GU = General_Utilities->new();
	
	method hierarchical_clustering( Generic_Graph $graph, Num $threshold ){
		
		$GU->user_info(3,"Performing hierarchical clustering\n");
		my @clusters;
		
		# Separate node labels into array of arrays, where each label has its own array
		foreach my $label ( @{$graph->node_labels} ){
			
			my @cluster;
			push(@cluster, $label);
			push(@clusters, \@cluster);
		}
		
		# Compute pairwise distances and if < threshold then merge nodes
		my $iter = 1;
		
		while ( 1 ){
			
			my $num_clusters = scalar(@clusters);
			my %cluster_pairs;
			for ( my $i=0; $i<$num_clusters; $i++){
				
				my %inter_cluster_distances;
				
				for ( my $j=$i+1; $j<$num_clusters; $j++){
					
					my $av_dist = $self->get_distance_between_two_clusters($graph, $clusters[$i], $clusters[$j]);
					#$GU->user_info(3,$i.":".$j."\t".$av_dist);
					$inter_cluster_distances{$j} = $av_dist;
				}
			
				$cluster_pairs{$i} = \%inter_cluster_distances;
			}
			# Merge clusters. Merge the two clusters that have shortest distance between one another
			my @keys = keys %cluster_pairs;
			@keys = sort{ $a <=> $b } @keys;
			
			my $min_key1;
			my $min_key2;
			my $min_dist;
			
			while ( my($key, $val) = each %cluster_pairs){
					
				#$GU->user_info(3,$key ."=>". $val."\n");
				while ( my($k, $v) = each %$val){
					
					if( !defined($min_dist) or $v < $min_dist){
						$min_dist = $v;
						$min_key1 = $key;
						$min_key2 = $k;
					}
				}
				
			}
			
			#$GU->user_info(3, "$iter: $min_dist\t$min_key1\t$min_key2\n");
			$GU->user_info(3,"$iter: $min_dist\t$min_key1\t$min_key2\n");
			if( $min_dist < $threshold){
				
				foreach ( @{$clusters[$min_key2]} ){
				
					push( @{$clusters[$min_key1]}, $_);
				}
				# Remove one of the cluster that we just merged
				splice( @clusters, $min_key2, 1 );
			}
			else{
				last;
			}
			
			
			$iter++;
			
		}
		
		return \@clusters;
		
	} # hierarchical_clustering #
	
	method get_distance_between_two_clusters(Generic_Graph $graph, ArrayRef $clusterA, ArrayRef $clusterB){
		
		my $av_distance = 0;
		my $count = 0;
		
		foreach my $labelA ( @$clusterA ){
				
			foreach my $labelB( @$clusterB ){
				
				my $d = $graph->retrieve_weight($labelA, $labelB);
				$av_distance += $d;
				$count++;
			}
		}
		
		return $av_distance/$count;
	} # get_distance_between_two_clusters #

} # Graph_Clustering #
