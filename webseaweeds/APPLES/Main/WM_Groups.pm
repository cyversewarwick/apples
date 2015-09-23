### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### WM_Groups Class ###
use MooseX::Declare;

class WM_Groups{
	
	use General_Utilities;
	use Generic_Graph;
	use Data::Dumper;
	use APPLES_Datatypes qw (PositiveNum);
	use File::Path qw (rmtree);
	use constant {FALSE => 0,TRUE => 1};
	
	my $WMU = WM_Utilities->new();
	my $GU = General_Utilities->new();
	
	has 'weight_matrices' => ( is => 'ro', isa => 'ArrayRef[Generic_Weight_Matrix]', required => 1, trigger => \&create_graph);
	has 'weight_matrices_hash' => ( is => 'rw', isa => 'HashRef');
	has 'graph' => (is => 'rw', isa => 'Generic_Graph'); # will be created via trigger associated when weight_matrices is initialised
	has 'threshold' => (is => 'rw', isa => PositiveNum);
	has 'groups' => (is => 'rw', isa => 'ArrayRef[ArrayRef]');
	
	sub create_graph(){
		
		my ($self, $arg) = @_;	
		
		my @wm_labels;
		# Get list of wm id's
		my $index = 0;
		my %hash;
		foreach( @{$self->weight_matrices} ){
			
			push(@wm_labels, $_->wm_identifier);
			
			$hash{$_->wm_identifier} = $_;
			
		}
		
		$self->weight_matrices_hash( \%hash );
		
		my $graph = Generic_Graph->new( node_labels => \@wm_labels);
		
		$self->{graph} = $graph;
		
		#my $matrix = $self->private_get_distance_matrix();
		
		#$self->graph->edge_weights($matrix);
		#$self->graph->print_edge_weights_table();
		
		#my $connections = $self->private_set_node_connections();
		
		#$self->graph->node_connections($connections);
		
	} # create_graph #
	
	method output_group_as_pssm_set_for_bifa_server( Str $set_name, Str $outputdir ){
		
		# Create dir if it does not exist
		
		if ( -d "$outputdir" ){
			
		}
		else{
			mkdir ($outputdir);
		}

		# print all weight matrices in group to this directory
		
		foreach my $wm ( @{$self->weight_matrices} ){
				
			$wm->print_to_custom_pssm_format($outputdir);
		}
		
		# Now print custom pssm set file
		open FH, ">$outputdir/$set_name.pssm_set" or die "Cannot open a file to print custom pssm set info to\n";
		my $printme;
		foreach my $wm ( @{$self->weight_matrices} ){
			$printme .= $wm->wm_identifier." ";
		}
		chomp($printme);
		print FH $printme."\n";
		close FH;
	} # output_group_as_pssm_set_for_bifa_server #
	
	method render_groups(Str $outputdir, Str $tvp){
		if ( -d "$outputdir" ){
			
		}
		else{
			mkdir ($outputdir);
		}
		
		# Make logos subfolder
		
		if ( -d "$outputdir/logos" ){
			
		}
		else{
			mkdir ($outputdir."/logos");
		}
		
#		open FH, "info" or die "Cannot open info file\n";
#
#        my %info = ();
#
#        my $id1;
#        my $id2;
#    
#        my %loc_info;
#    
#        while (<FH>) {
#
#            chomp;
#            my @split = split( /,/, $_ );
#
#            if( $split[0] eq 'LOC' ) {
#                $id1 = $split[1];
#                $id2 = $split[2];
#                $info{$id1."-".$id2}{'length'} = $split[3];
#                $info{$id1."-".$id2}{'reverse'} = $split[4];
#            }
#    
#            elsif( $split[0] eq 'INFO'){
#                $info{$id1."-".$id2}{'front1'} = $split[1];
#                $info{$id1."-".$id2}{'back1'} = $split[2];
#                $info{$id1."-".$id2}{'front2'} = $split[3];
#                $info{$id1."-".$id2}{'back2'} = $split[4];
#            }
#        }
#
#        close FH;
		
		my $html_filepath = $outputdir."/groups.html";
		open HTML, ">$html_filepath" or die;
		print HTML "<html>\n";
		print HTML "<head><title>PSSM Clusters</title></head>\n";
		print HTML "<body>\n";
		print HTML "<center><h1>Clusters with threshold: </h1>\n";		
		
		my $group_num = 1;
		foreach my $group (@{$self->groups()}) {
			my $max_length;
			
			foreach my $id (@$group) {
				
			}
			
			print HTML "<table  border=\"1\" bordercolor=\"#FFCC00\" cellpadding=\"3\">\n";
			print HTML "<tr><th colspan=\"2\">Group: ".$group_num."</th></tr>\n";
			$group_num++;
			
			foreach my $id (@$group) {
			$GU->user_info(3,"Rendering logo for ID: $id \n");
			    print HTML "<tr><td>";
				my $gwm = ${$self->weight_matrices_hash}{$id};
				
				my $info;
				if ($gwm->isa('BiFa_Server_Weight_Matrix')) {
					$info = $WMU->get_bifa_pssm_info($id, $tvp);
					print HTML $id."\t".$$info{name}."\t";
				} else {
					print HTML $id."\t(no name)\t";
				}
				
			    print HTML "</td><td>";
				
				my $logos_filepath = $outputdir."/logos/";
				my $logo_fp = $gwm->generate_logo($logos_filepath);
				my $logo_filename = $id.".png";

				print HTML "<img src= \"logos/$logo_filename\" width=\"144\" height=\"50\" />";
				print HTML "</td></tr>\n";
	
			}
			
			print HTML "</table>\n";
            print HTML "<br>\n";
		}
		
		
		print HTML "</center>\n";
		print HTML "</body>";
		print HTML "</html>";
		close HTML;
				
	} # render_groups #
	
	
	method cluster_weight_matrices(Int $sorting){
			
		$self->private_cluster_this($sorting);
	} # cluster_weight_matrices #
	
	method private_cluster_this(Int $sorting){
		
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		my $wms_fp = $WMU->print_multiple_pssms_to_single_file($self->weight_matrices, $tempdir);
		my $num_pssms_to_compare = scalar(@{$self->weight_matrices});
		my $threshold = $self->threshold();
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_matrix_distance_c_exe = $APPLES_C_binaries."/pssm_distance_functionsREV789";
		
		$GU->user_info(3, "$path_to_matrix_distance_c_exe -f C -i $wms_fp -n $num_pssms_to_compare -o $tempdir -d h -t $threshold\n");
		system("$path_to_matrix_distance_c_exe", "-f", "C", "-i", "$wms_fp", 
		              "-n", "$num_pssms_to_compare", "-o", "$tempdir", 
		              "-d", "h", "-t", "$threshold", "-s", "$sorting") == 0 || die "System error!";
		
		open FH, "$tempdir/final_clusters" or die "Cannot open pairwise distance matrix file\n";
		my $cluster_number;
		my $number_of_members;
		my @self_weight_matrices_index;
		my @groups;
		
		while(<FH>) {
			
			chomp;
			my @split = split('\t', $_);
			
			if( $split[0] eq 'CL' ){
				$cluster_number = $split[1];
			}
			elsif( $split[0] eq 'NM'){
				$number_of_members = $split[1];
			}
			elsif( $split[0] eq 'MB' ){
				my @wm_indexs = split(',', $split[1]);
				my @group;
				
				foreach my $i ( @wm_indexs ) {
					push(@group, ${$self->weight_matrices}[$i]->wm_identifier );
				}
				push(@groups, \@group);
					
			}
			
		}
		#rmtree($tempdir);
		$self->groups( \@groups )
	} # private_cluster_this #
	
	method private_set_node_connections(){
			
		my @matrix;
		
		for( my $i = 0; $i < scalar(@{$self->weight_matrices}); $i++){
			
			for( my $j = 0; $j < scalar(@{$self->weight_matrices}); $j++){
				
				if( $self->weight_matrices->[$i]->wm_identifier eq $self->weight_matrices->[$j]->wm_identifier ){
					$matrix[$i][$j] = 0;
				}
				else{
					$matrix[$i][$j] = 1;
				}
			}
		}
		
		return \@matrix;
	} # private_set_node_connections #
	
	method private_get_distance_matrix(){
			
		my @matrix;
		my $addme=0;
		
		my $tempdir = $GU->get_temp_random_directory(FALSE);
		my $wms_fp = $WMU->print_multiple_pssms_to_single_file($self->weight_matrices, $tempdir);
		my $num_pssms_to_compare = scalar(@{$self->weight_matrices});
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_matrix_distance_c_exe = $APPLES_C_binaries."/produce_pssm_distance_matrixREV493";
		
		#my $path_to_matrix_distance_c_exe = "/cluster/richardhickman/jan4_apples_working_copy/Main/C/produce_pssm_distance_matrix";
		print "$path_to_matrix_distance_c_exe $wms_fp $num_pssms_to_compare $tempdir\n";
		system("$path_to_matrix_distance_c_exe", "$wms_fp", "$num_pssms_to_compare", "$tempdir") == 0 || die "System error!";
		
		open FH, $tempdir."/pairwise_distance_matrix" or die "Cannot open pairwise distance matrix file\n";
		
		while(<FH>){
			chomp;
			
			my @split = split('\t', $_);
			
			push(@matrix, \@split);
			
		}
		rmtree($tempdir);
		return \@matrix;
	} # private_get_distance_matrix #
	
	method compute_node_node_distance(Generic_Weight_Matrix $wm, ArrayRef[Generic_Weight_Matrix] $GWMA){
		
		#die "Method not yet implemented\n"; # (RH)

	    die "This method cannot work - needs checking."; # Sascha: $tempdir is not used in the system command below - cannot be right

		my $tempdir = $GU->get_temp_random_directory(FALSE);
		
		# Print pssms to temporary directory
		my $wm_fp = $WMU->print_pssm_to_file_for_empirical_scoring($wm, $tempdir);
		my $wm_array_fp = $WMU->print_multiple_pssms_to_single_file($GWMA, $tempdir);
		my $num_pssms_to_compare = scalar(@{$GWMA});
		
		my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
		my %APPLES_config = $APPLES_conf->getall();
		my $APPLES_C_binaries = $APPLES_config{path_to_APPLES_C_binaries};
		my $path_to_matrix_distance_c_exe = $APPLES_C_binaries."/matrix_distanceREV493";
		
		# Get distance between pssms
		#my $path_to_matrix_distance_c_exe = "/cluster/richardhickman/jan4_apples_working_copy/Main/C/matrix_distance";
		print "$path_to_matrix_distance_c_exe\t". "$wm_fp\t". "$wm_array_fp\t"."$num_pssms_to_compare\n";
		
		system("$path_to_matrix_distance_c_exe", "$wm_fp", "$wm_array_fp", "$num_pssms_to_compare") == 0 || die "System error!";
		my $output = `$path_to_matrix_distance_c_exe $wm_fp $wm_array_fp $num_pssms_to_compare`;
		my @split_output = split('X', $output);
		my @shortest_distances;
		
		foreach ( @split_output){
				
			my @distances = split(',', $_);
			
			my @ranked_distances = sort { $a <=> $b } @distances;
			
			my $shortest_distance = $ranked_distances[0];
			push(@shortest_distances, $shortest_distance);
		}
		
    	        rmtree($tempdir);
		return \@shortest_distances;
	} # compute_node_node_distance # 
	
} # WM_Groups #
