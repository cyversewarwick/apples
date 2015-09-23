### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReCo_Maker class ###

use MooseX::Declare;


class ReCo_Maker {
    use APPLES_Datatypes qw (Boolean);
    use Statistics::R;
    use Parameters;
    use ReCo;
    use Generic_Pattern_Matching_Model;
    use constant {FALSE => 0,TRUE => 1};
    use File::Temp qw (tempdir);
    use File::Path qw (rmtree);
	
    method reco_maker_binomial(Genomic_Interval_Set $gi_set, ArrayRef[Generic_Pattern_Matching_Model] $pattern_models, ReCo_Parameters $reco_parameters) {
	# returns a ReCo, uses binomial overrepresentation test

	# score patterns on sequences, get binding molecules for each pattern, and take maximum score for each molecule
	my @emptyarray = ();
	my $input_ReCo = ReCo->new(molecules => \@emptyarray);
	my %molecules_seen_before;
	foreach my $pattern_model (@{$pattern_models}) {
	    my $result = $pattern_model->binomial_overrepresentation_test($gi_set,$reco_parameters->max_number_of_sites);
	    my $p_value = $result->{PVALUE};
	    my $sampling_likelihood = $self->private_convert_p_value_into_sampling_likelihood($p_value,
											      $reco_parameters->minimum_node_probability,
											      $reco_parameters->conversion_parameters);
	    my @binding_molecules = $pattern_model->pattern->list_binding_molecules();
	    foreach my $molecule (@binding_molecules) {
		my $unique_ID = $molecule->unique_ID();
		if (!$molecules_seen_before{$unique_ID}) {
		    $input_ReCo->add_molecule($molecule);
		}
		$molecules_seen_before{$unique_ID} = TRUE;
		$input_ReCo->lower_bound_on_likelihood($unique_ID,$sampling_likelihood);
	    }
	}
	# get nuclear localisation data
	my $nuclear_localisation_ReCo = ReCo->new(molecules => \@emptyarray);
        # do sampling
	my $output_ReCo = $self->reco_maker($input_ReCo, $nuclear_localisation_ReCo, $reco_parameters);
	return $output_ReCo;
    } # reco_maker_binomial #

    method reco_maker(ReCo $input_ReCo, ReCo $nuclear_localisation_ReCo, ReCo_Parameters $reco_parameters) {

	my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
        my %APPLES_config = $APPLES_conf->getall();
	my $path_to_APPLES_R_programs = $APPLES_config{path_to_APPLES_R_programs};

	# Create directory to write results files to.
	my $GU = General_Utilities->new();
        my $tempdir = $GU->get_temp_random_directory(FALSE);
        
        if ($^O eq 'MSWin32') {
          # Temporary directory created in windows uses backslashes, which confuses R
          $tempdir =~ s/\\/\//g;}
        

        # do sampling using Perl to R interface
        my $R = Statistics::R->new();
        $R->startR();
        # Pass input parameters to R
	$R->send(qq /path_to_code                          <- \"/.$path_to_APPLES_R_programs.qq /\"/ );
        $R->send(qq /p_p_network                           <- \"/.$reco_parameters->p_p_network.qq /\"/); 
           # human_intact or no_network
        $R->send(qq /minimum_node_probability              <- /.$reco_parameters->minimum_node_probability); 
        $R->send(qq /number_of_sampling_iterations         <- /.$reco_parameters->number_of_sampling_iterations); 
        #$R->send(qq /max_number_of_sites                   <- /.$reco_parameters->max_number_of_sites); 
           # not needed in R?
        #$R->send(qq /conversion_parameters                 <- /.$reco_parameters->conversion_parameters); 
           # not needed in R?
        $R->send(qq /scale_number_of_experiments           <- /.$reco_parameters->scale_number_of_experiments); 
        $R->send(qq /scale_number_of_experimental_methods  <- /.$reco_parameters->scale_number_of_experimental_methods); 
        $R->send(qq /weight_for_two_logistic_functions     <- /.$reco_parameters->weight_for_two_logistic_functions); 
	$R->send(qq /results_tempdir			   <- \"/.$tempdir.qq /\"/);

	$R->send(qq /print(results_tempdir)/);
	my $R_return = $R->read;
	print "\nReturned from R $R_return\n\n";


	# Apply node weights by writing out the $input_ReCo into R.
	$R->send(qq (nodeProbs.v <- c() ));
	$R->send(qq (nodeIds.v <- c() ));
	my $molecules = $input_ReCo->molecules;
	foreach my $molecule (@$molecules){
 	    my $molecule_ID = $molecule->unique_ID();
	    my $likelihood = $molecule->likelihood();
	    $R->send(qq /nodeId <- \"/.$molecule_ID.qq /\"/ );
	    $R->send(qq /nodeIds.v <- c(nodeIds.v, nodeId)/ );
	    $R->send(qq /nodeProb <- /.$likelihood );
	    $R->send(qq /nodeProbs.v <- c(nodeProbs.v, nodeProb)/ );
	}


	# With all the parameters passed to R, load the script with 
	# the functions and commands to do the sampling
	print(qq /source("/.$path_to_APPLES_R_programs.qq/\/ReCo_Sampler.R")\n/);
        $R->send(qq /source("/.$path_to_APPLES_R_programs.qq/\/ReCo_Sampler.R")/);

	# Return the new ReCo from R (ReCo_scores.v) CAN'T GET THIS TO WORK.
	#print "\nPrinting ReCo.\n\n";
	#$R->send(qq /print(ReCo_scores.v)/);
	my $R_source_return = $R->read;
	print "Returned from R $R_source_return\n\n";

        $R->stopR() ; 
	
	# Read the ReCo from file
	my @emptyarray = ();
	my $reco = ReCo->new(molecules => \@emptyarray);
	
	# May want to change this to a temp dir
	my $fname = "$tempdir/ReCo.tsv";
	open(INFILE, $fname)
		or die "Can't open ReCo result $fname\n";
	while (my $line = <INFILE>) {
		if ( $line =~ /(\w+)\t(\d+\.?\d?+)/ ){
			print "Protein  $1\n";
			my $protein = Protein->new(unique_ID => $1, likelihood => $2);
			$reco->add_molecule($protein);
		}
	}
	rmtree($tempdir);
	return $reco;
    } # reco_maker #

    method private_convert_p_value_into_sampling_likelihood(Num $p_value, Num $minimum_likelihood, P_Value_Sampling_Likelihood_Conversion_Parameters $conversion_parameters) {
	# uses logistic function to turn a p-value into a likelihood used in sampling

	my $likelihood = 2*(1/(1+$p_value)-0.5); # $pvalue is abbreviation of exp(log($p_value))	
	if ($likelihood<$minimum_likelihood) {
	    $likelihood = $minimum_likelihood;
	}
	return $likelihood;
    } # private_convert_p_value_into_sampling_likelihood #
    
} # ReCo_Maker #
