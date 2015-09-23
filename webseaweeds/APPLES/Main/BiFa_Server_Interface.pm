use SOAP::Lite;

use MooseX::Declare;

class BiFa_Server_Interface{
use Config::General;	

#	Constants to use as the algorithm parameter.  
#	After "use BiFa_Server_Interface;" in the client program
#	they can be referred to using: BiFa_Server_Interface::ALG_BAYESIAN 

use constant ALG_OTT => 0;
use constant ALG_BAYESIAN => 1;

has 'sIf' => (is => 'rw', isa => 'SOAP::Lite');
has 'errorStr' => (is => 'rw');
has 'tfVersion' => (is => 'rw');

method PssmFreqs($pssmName,$dataVersion)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};
	my $var1 = SOAP::Data->type('string') ->name('pssmName') ->value($pssmName);
	my $var = SOAP::Data->type('complex') ->value($var1) ;
	
	my $res = $self -> sIf ->PssmFreqs($var);
	my $result = $res -> result();
	
	if ($self -> IsFault($res)) {return undef;}

	my @output;

	#	Find out if an array has been returned 
	if (ref($result) eq "stringArray")
		{
			#	Convert from strings to a two dimensional araay
			my $freqs = $result -> {'string'};
			foreach my $freq (@$freqs)
			{
				my @z = split(',',$freq);
				push (@output,\@z);			
			}
			return @output;
		}
	else
		{ return ();}
}

method PssmCounts($pssmName,$dataVersion)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};
	my $var1 = SOAP::Data->type('string') ->name('pssmName') ->value($pssmName);
	my $var = SOAP::Data->type('complex') ->value($var1) ;
	
	my $res = $self -> sIf ->PssmCounts($var);
	my $result = $res -> result();
	
	if ($self -> IsFault($res)) {return 0;}

	my @output;

	#	Find out if an array has been returned 
	if (ref($result) eq "stringArray")
		{
			#	Convert from strings to a two dimensional araay
			my $freqs = $result -> {'string'};
			foreach my $freq (@$freqs)
			{
				my @z = split(',',$freq);
				push (@output,\@z);			
			}
			return @output;
		}
	else
		{ return ();}
}


#	Initialises the soap interface.  For efficient usage a BiFa_Server_Interface
#	object should be created, initialised and reused.
method init ($username,$password){
	
	# If run on the wsbc server or IBM cluster use the server address, otherwise use local
        # address, assuming that an ssh tunnel has been setup (by putty?) to the wsbc server
        # (both port 4338 and 4339 on the XServe cluster are accessible for all compute nodes of
        # the IBM cluster, note that without using SSH your password will be sent from the IBM cluster
        # to the wsbc server without encryption).
        my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
        my %APPLES_config = $APPLES_conf->getall();
        my $server = $APPLES_config{BiFa_Server_Name};
        my $port = $APPLES_config{BiFa_Server_Port};

	my $Uri = $server.':'.$port;
#	if ($^O eq "MSWin32")
#		{$Uri = 'localhost:'.portNo;}

	$self -> sIf(SOAP::Lite 
			->uri('http://'.$Uri)
			->proxy('http://'.$username.':'.$password.'@'.$Uri));

	#	Now get the cookie, using the 'returningString' method
	my $var = SOAP::Data->type('string') ->name('key')->value('connection_info');
	my $var2 = SOAP::Data->type('complex') ->value($var);
	my $result = $self -> sIf ->returningString($var2);

	if ($self -> IsFault($result)) {return 0;}

	#	Extract the cookie and the transfac version from the returned message and store 
	#	it for subsequent connections
	my $value = $result -> result();
	my $cookie = "";
	if (substr($value,0,12) eq "established:") {
		$cookie = substr($value,12);
		my $p = index($cookie,":");
		if ($p != -1) {
			$self ->tfVersion(substr($cookie,$p+1));
			$cookie = substr($cookie,0,$p); 
			}}
	else 
		{return 0;}

	#	And change the soap interface data to use the cookie rather than the password	
	$self -> sIf ->proxy('http://'.$username.':'.$cookie.'@'.$Uri);
	$self -> sIf ->proxy ->timeout(500);
	return 1;
}  



method IsFault($result)
{
	if ($result -> fault)
	{
		$self ->errorStr (join ',',$result->faultcode,$result->faultdetail);
		return 1;
	}
	return 0;
}
#	If the previous operation failed, this returns a string
#	describing the error
method ErrorString()
{
	return $self -> errorStr;
}

#	returns a hash value with members {major} and {minor}
method ServerVersion()
{
	#  Have to pass in an explicit variable to the underlying method
	#	rather than no variable at all, otherwise the xml syntax within 
	#	soap request is incorrect.
	#   I suspect that this is because the perl code does not have access to the wdsl
	#	to create the correct syntax
	my $var2 = SOAP::Data->type('none');
	my $result = $self -> sIf ->serverVersion($var2);
	
	if ($self -> IsFault($result)) {return undef;}
	
	return $result -> result();
}
#	returns a hash value with members {major} and {minor}

method CustomPssmVersion()
{
	my $var2 = SOAP::Data->type('none');
	my $result = $self -> sIf ->customPssmVersion($var2);
	
	if ($self -> IsFault($result)) {return undef;}
	
	return $result -> result();
}


#	returns the maximum chain/Maximum number of sequences.
method MaxChainMaxNumSequences()
{
	my $result = $self -> sIf ->maxChainMaxNumSequences(SOAP::Data->type('none'));
	if ($self -> IsFault($result)) {return 0;}
	return $result -> result();
}

#	Returns a string indicating the current transfac version
method TransfacVersion()
{
	my $result = $self -> sIf ->transfacVersion( SOAP::Data->type('none'));
	
	if ($self -> IsFault($result)) {return undef;}
	
	return $result -> result();
}

#	Returns an array of strings, one for each PssmSetName
method PssmSetNames()
{
	my $result = $self -> sIf ->PssmSetNames(SOAP::Data->type('none'));

	if ($self -> IsFault($result)) {return undef;}
	
	# The need for the @{} syntax took hours to sort out, it extracts the underlying array
	if (defined $result -> result() -> {'string'})
		{ return @{$result -> result() -> {'string'}}; }
	else
		#	So have to explictly return an empty set of hits
		{ return ();}
}

method Pssms($useConsensusSequences, $matrixSpecies, $matrixNameMatch, $dataVersion,$pssmSets)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};
	my $var1 = SOAP::Data->type('boolean') ->name('useConsensusSequences')->value($useConsensusSequences);
	my $var2 = SOAP::Data->type('string') ->name('matrixSpecies')->value($matrixSpecies);
	my $var3 = SOAP::Data->type('string') ->name('matrixNameMatch')->value($matrixNameMatch);
	#	And put them into the array holding the soap data to be output
	my $var = SOAP::Data->type('complex') ->value($var1,$var2,$var3);

	if (defined $pssmSets)
	{
		my $i = @$pssmSets;
		if ($i > 0)
			{ my $var4 = SOAP::Data->type('tns:stringArray') ->name('pssmSets') ->value($pssmSets);
			$var ->value($var ->value,$var4); }
	}	
	my $result = $self -> sIf ->Pssms($var);
	
	if ($self -> IsFault($result)) {return undef;}

	if (defined $result -> result() -> {'string'})
		{ return @{$result -> result() -> {'string'}}; }
	else
		#	So have to explictly return an empty set of hits
		{ return ();}
}



#	returns an array containg the score for a the match for a specific pssm at every
#	location in a sequence
method scorePssmOnSequence	($pssm_name,$sequence)
{
	my $var1 = SOAP::Data->type('string') ->name('pssm_name')->value($pssm_name);
	my $var2 = SOAP::Data->type('string') ->name('sequence')->value($sequence);
	my $var3 = SOAP::Data->type('float') ->name('threshold')->value(0.0);
	my $var = SOAP::Data->type('complex') ->value($var1,$var2,$var3) ;
	
	my $result = $self -> sIf ->scorePssmOnSequence($var);
	
	if ($self -> IsFault($result)) {return undef;}
	
	for my $r (@{$result -> result() -> {'float'}})
		{
			if ($r == -1 )	{ $r = undef;}
		}
	
	return @{$result -> result() -> {'float'}};
}

#	As above, but multiple pssms or multiple sequences can be analysed at the same time
#	(but not both) The algorithm parameter is either BiFa_Server_Interface::ALG_OTT or 
#	BiFa_Server_Interface::ALG_BAYSIAN.
#	BiFa_Server_Interface::ALG_BAYSIAN is the algorithm supported by scorePssmOnSequence and is the
#	only algorithm currently supported in this method.
method scorePssmsOnSequences	($pssmNames,$sequences,$algorithm,$dataVersion)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};

	my $var1 = SOAP::Data->type('tns:stringArray') ->name('pssmNames') ->value($pssmNames);
	my $var2 = SOAP::Data->type('tns:stringArray') ->name('sequences') ->value($sequences);
	my $var3 = SOAP::Data->type('integer') ->name('algorithm') ->value($algorithm);

	my $var = SOAP::Data->type('complex') ->value($var1,$var2,$var3) ;
	
	my $result = $self -> sIf ->scorePssmsOnSequences($var);
	
	if ($self -> IsFault($result)) {return undef;}
	
	my $res = $result -> result();

	my @output;

	#	Find out if we have an array, or an array of arrays 
	if (ref($res) eq "stringArray")
		{
			my $s = $res -> {'string'};
			if (ref($s) eq 'ARRAY')
			{ 
				for my $i (0.. @$s-1)
				{
					my @z = split(',',$$s[$i]);
					for my $r (@z)
					{
						if ($r == -1 ) { $r = undef;}
					}
					$$s[$i]="";
					push (@output,\@z);
				}			
			}
			else
			{	
				my @z = split(',',$s);
				for my $r (@z)
				{
					if ($r == -1 ) {
					  $r = undef;
					}
				}
				push (@output,\@z);			
			}
			return @output;
		}
	else
		#	No arrays, so have to explictly return an empty set of hits
		#	which indicates a fault
		{ return ();}
		
}

#	Get more detailed info related to a specific PSSM, as indicated by its 
#	M00001, R00001 or custom pssm identifier
method PssmInfo	($pssmName,$dataVersion)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};
	my $var1 = SOAP::Data->type('string') ->name('pssmName') ->value($pssmName);

	my $var = SOAP::Data->type('complex') ->value($var1) ;
	
	my $result = $self -> sIf ->PssmInfo($var);
	
	if ($self -> IsFault($result)) {return undef;}

	#	The factors are returned as a list, and Perl and SOAP between them
	#	put them in an array if there is nore than one of them, or
	#	return them as a scalar if there is only one,  making it 
	#	difficult to write client code that will cope with both.
	#	This code spots when there is a scalar and converts it into an array
	#	with one entry	
	my $factors = $result -> result() ->{'factors'}->{'string'};
	if (defined($factors)) {
		if (ref($factors) ne "ARRAY"){
			$result -> result() ->{'factors'}->{'string'} = [$factors];
		}
	}
	else {
		$result -> result() ->{'factors'}->{'string'} = [];
	}
	
	return $result -> result();
}


#	Creates an svg picture showing likely binding regions on a stretch of DNA
method bifa	($sequence, $threshold, $title, $algorithm, $showLabels,  $phyloSequences,$useConsensusSequences, $matrixSpecies,$phyloThreshold, $matrixNameMatch, $useOldAlgorithm, $useCumulativeLikelihoods, $pssmSets)
{

	#	Create the soap parameters.  types and names as per wdsl	
	my $var1= SOAP::Data->type('string') ->name('sequence')->value($sequence);
	my $var2 = SOAP::Data->type('float') ->name('threshold')->value($threshold);
	my $var3 = SOAP::Data->type('string') ->name('title')->value($title);
	my $var4 = SOAP::Data->type('int') ->name('algorithm')->value($algorithm);
	my $var5 = SOAP::Data->type('boolean') ->name('showLabels')->value($showLabels);
	my $var7 = SOAP::Data->type('boolean') ->name('useConsensusSequences')->value($useConsensusSequences);
	my $var8 = SOAP::Data->type('string') ->name('matrixSpecies')->value($matrixSpecies);
	my $var9 = SOAP::Data->type('float') ->name('phyloThreshold')->value($phyloThreshold);
	my $var10 = SOAP::Data->type('string') ->name('matrixNameMatch')->value($matrixNameMatch);
	my $var11 = SOAP::Data->type('boolean') ->name('useOldAlgorithm')->value($useOldAlgorithm);
	my $var12 = SOAP::Data->type('boolean') ->name('useCumulativeLikelihoods')->value($useCumulativeLikelihoods);
	
	#	And put them into the array holding the soap data to be output
	my $var = SOAP::Data->type('complex') ->value($var1,$var2,$var3,$var4,$var5,$var7,$var8,$var9,$var10,$var11,$var12);

	#Only add the phylo sequences that are required.  The Phylosequences are an array, and while most parameters
	# are passed 'by value' Perl automatically decides to pass by arrays by reference.  Inside the subroutine some
	#  array operations are unaffacted by the fact that this is a reference.
	#  Getting the size of the array is not, and the array has to be explicitly deferenced
	#  In the following line '@$XX' means get the array referenced by XX, and by setting it to a scalar we get the length
	#  So, the following line sets i to the length of the array phylosequences.   But that was obvious (not)
	my $i = @$phyloSequences;
	if ($i > 0)
		#  Only add the sequences to the soap call if there are some

		{ 
		my $var6 = SOAP::Data->type('tns:stringArray') ->name('phyloSequences') ->value($phyloSequences);
		$var ->value($var ->value,$var6); }

	if (defined $pssmSets)
	{
		my $i = @$pssmSets;
		if ($i > 0)
			{ my $var13 = SOAP::Data->type('tns:stringArray') ->name('pssmSets') ->value($pssmSets);
			$var ->value($var ->value,$var13); }
	}	
	my $result = $self -> sIf ->BiFaAnalyser($var);
	
	if ($self -> IsFault($result)) 
	{
		return undef;
	}

	return $result -> result();
}

#	As 'bifa' except that the result is in the form of an array of hits.
method bifaHits	($sequence, $threshold, $algorithm, $phyloSequences,$useConsensusSequences, $matrixSpecies,$phyloThreshold, $matrixNameMatch, $useCumulativeLikelihoods,$dataVersion,$pssmSets)
{
	if ($dataVersion ne $self->tfVersion)
		{ return 0;};

	#	Create the soap parameters.  types and names as per wdsl	
	my $var1= SOAP::Data->type('string') ->name('sequence')->value($sequence);
	my $var2 = SOAP::Data->type('float') ->name('threshold')->value($threshold);
	my $var3 = SOAP::Data->type('int') ->name('algorithm')->value($algorithm);
	my $var5 = SOAP::Data->type('boolean') ->name('useConsensusSequences')->value($useConsensusSequences);
	my $var6 = SOAP::Data->type('string') ->name('matrixSpecies')->value($matrixSpecies);
	my $var7 = SOAP::Data->type('float') ->name('phyloThreshold')->value($phyloThreshold);
	my $var8 = SOAP::Data->type('string') ->name('matrixNameMatch')->value($matrixNameMatch);
	my $var9 = SOAP::Data->type('boolean') ->name('useCumulativeLikelihoods')->value($useCumulativeLikelihoods);
	
	#	And put them into the array holding the soap data to be output
	my $var = SOAP::Data->type('complex') ->value($var1,$var2,$var3,$var5,$var6,$var7,$var8,$var9);

	#  Only add the phylo sequences that are required.  The Phylosequences are an array, and while most parameters
	#  are passed 'by value' Perl automatically decides to pass by arrays by reference.  Inside the subroutine some
	#  array operations are unaffacted by the fact that this is a reference.
	#  Getting the size of the array is not, and the array has to be explicitly deferenced
	#  In the following line '@$XX' means get the array referenced by XX, and by setting it to a scalar we get the length
	#  So, the following line sets i to the length of the array phylosequences.   But that was obvious (not)
	my $i = @$phyloSequences;
	if ($i > 0)
		#  Only add the sequences to the soap call if there are some
		{ my $var4 = SOAP::Data->type('tns:stringArray') ->name('phyloSequences') ->value($phyloSequences);
		$var ->value($var ->value,$var4); }

	if (defined $pssmSets)
	{
		my $i = @$pssmSets;
		if ($i > 0)
			{ my $var11 = SOAP::Data->type('tns:stringArray') ->name('pssmSets') ->value($pssmSets);
			$var ->value($var ->value,$var11); }
	}	
	
	my $result = $self -> sIf ->BiFaHits($var);
	
	if ($self -> IsFault($result)) {return undef;}

	#	Only defined if there are a number of hits 
	if (defined $result -> result() -> {'BiFaHit'})
	{
		return @{$result -> result() -> {'BiFaHit'}};
	}	
	else
	{
		#	So have to explictly return an empty set of hits
		return ();
	}
}
}

