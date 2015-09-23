### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### ReRe_Maker class ###
use MooseX::Declare;

class ReRe_Maker {
	use Genome_DB_Utilities;
	use Parameters;
	use ReRe;
	use ReMo_Set_Maker;
	use Data::Dumper;
	use General_Utilities;
	use Exception;
	use constant {FALSE => 0,
			TRUE	=> 1};
	my $GU = General_Utilities->new();
	
	method rere_maker(ReRe_Constructor_Parameters $parameters, Reg_Loc $core_reg_loc) {
		my @all_remo_sets = ();
		foreach my $parameter_set (@{$parameters->remo_set_constructor_parameters}) {
			my @remo_sets;
			if ($parameter_set->isa('ReMo_Set_Core_Promoter_Constructor_Parameters')) {
				my $remo_set = 	ReMo_Set_Maker->new()->make_remo_sets_from_core_promoter($parameter_set, $core_reg_loc);
				push (@remo_sets, $remo_set);
			}
			elsif ($parameter_set->isa('ReMo_Set_Phylogenetic_Constructor_Parameters')) {
				@remo_sets = ReMo_Set_Maker->new()->make_phylogenetic_remo_sets($parameter_set, $core_reg_loc);
			}
			else {
				$GU->user_info( 1, "invalid parameter set\n" );
			}
			@all_remo_sets = (@all_remo_sets,@remo_sets);
		}
		my $rere = ReRe->new(remo_sets => \@all_remo_sets, core_reg_loc => $core_reg_loc); 
		
		return $rere;
	} # rere_maker #

	method rere_maker_through_job_handler (ReRe_Constructor_Parameters $parameters, Reg_Loc $core_reg_loc) {
	  # will not do any exception handling

	  my $object = ReRe_Maker->new();
	  my $function = 'rere_maker';
	  my @parameters = ($parameters,$core_reg_loc);
	  my $high_memory = TRUE;
	  my @result = $GU->standard_call_to_job_handler_without_exception_handling($object,$function,\@parameters,$high_memory,FALSE);
	  return $result[0];
        } # rere_maker_through_job_handler #

} # ReRe_Maker #
