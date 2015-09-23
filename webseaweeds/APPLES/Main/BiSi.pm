### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### BiSi Class: describes a Binding Site (in a single species) ###
use MooseX::Declare;
use Serializable;

class BiSi extends Serializable {
	use Parameters;
	use Genomic_Interval_Set;
	use General_Utilities;
	use Evidence;
	use Data::Dumper;
	
	has 'genomic_interval_set'		=>	(is => 'rw', 
							 isa => 'Genomic_Interval_Set', 
							 required => 1); 				 # must have a Genomic_Interval
	has 'tissue'				=> 	(is => 'rw', isa => 'Str');
	has 'single_site_pvalue'	=>	(is => 'rw', isa => 'Num');
	has 'comment'				=>	(is => 'rw', isa => 'Str');
	has 'bifa_assignment'		=>	(is	=> 'rw', isa => 'ArrayRef[Str]');
	has 'bisi_score'			=>	(is => 'rw', isa => 'Num');      
	has 'bisi_evidence'			=>	(is => 'rw', isa => 'ArrayRef[Evidence]');   
	has 'bisi_orthologs' 		=> 	(is => 'rw', isa => 'ArrayRef[BiSi]');
	has 'label'                 =>  (is => 'rw', isa => 'Str');
	
	my $GU = General_Utilities->new();

	method bisi_bifa_assigner(Str $identifier) {			# assigns a binding site to potential binding factors

		die 'method not implemented yet!';
		
                my $BiFa_Assignment;
		return $BiFa_Assignment;
	} # bisi_bifa_assigner #
		
	method get_bisi_evidence {			
	
		die 'method not implemented yet!';
	
		my $BiSi_Evidence;
	} # get_bisi_evidence #
		
	method get_bisi_tissue {			
	
		die 'method not implemented yet!';
	
		my $Tissue;
	} # get_bisi_tissue #

	method render_text {
	    if (defined $self->single_site_pvalue) {
		$GU->user_info(1,"p-value: ".$self->single_site_pvalue."\n");
	    }
	    my @sequences = $self->genomic_interval_set->return_unmasked_sequences;
	    my $count = 1;
	    foreach my $sequence (@sequences) {
		$GU->user_info(1,$count.': '.$sequence."\n");
		$count++
	    }
	} # render_text #
	
#	method to_array() {
#		my @gi_array;
#        
#        foreach my $gi (@{$self->genomic_interval_set->genomic_interval_set}) {
#
#           my $bisi_hash = {'genomic_interval'=>$gi->to_hash(), 'single_site_pvalue' => $self->single_site_pvalue, 'label' => $self->label};
#           push(@gi_array, $bisi_hash); 
#        }
#        
#        return @gi_array;
#	}
#	
#	override to_hash() {
#        my @gi_array;
#        
#        foreach my $gi (@{$self->genomic_interval_set->genomic_interval_set}) {
#
#           my $bisi_hash = {'genomic_interval'=>$gi->to_hash(), {'single_site_pvalue' => $self->single_site_pvalue}, {'label' => $self->label}};
#           push(@gi_array, $bisi_hash); 
#        }
#        
#        my $hash = {'bisi' => [{'single_site_pvalue' => $self->single_site_pvalue}]};
#        
#        return $hash;  
#	}

} # BiSi #



