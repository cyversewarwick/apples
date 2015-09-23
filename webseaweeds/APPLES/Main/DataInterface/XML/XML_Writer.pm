## (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### XML_Writer ###
# Class for outputting various data structures to XML for re-use in webservices.
#

BEGIN {
  our $global_print_level = "1";
  push(@INC,'../');
}


use MooseX::Declare;

class XML_Writer {
	use Conservation_Profiles;
	use Genomic_Interval;
	
	use Data::Dumper;
	use List::Util qw(min max);

	use constant { FALSE => 0, TRUE => 1 };
	use APPLES_Datatypes qw (Boolean);
	use XML::Simple;
#	use XML::Validator::Schema;
#    	use XML::SAX::ParserFactory;

	my $GU = General_Utilities->new();

	has 'XML' => (is => 'rw', isa => 'HashRef');

	# write conservation profiles to XML file.
	# XML format for genomic intervals is as follows:
	# <dna> <!-- (not written by this function, we just write the fragment) -->
	# <!-- each dna element can contain multiple chromosome fragments -->
	# <chromosomefragment>
	#  <!-- start of the fragment, int, absolute coordinates -->
	#  <start>...</start>
	#  <!-- end of the fragment, int, absolute coordinates -->
	#  <end>...</end> <!-- absolute coordinates -->
	#  <!-- strandedness. this was used to compute start and end
	#       coordinates and is for information only -->
	#  <strand>[positive|negative]</strand>
	#
	# </chromosomefragment>
	#<!-- ... -->
	#</dna> <!-- (not written by this function, we just write the fragment) -->
	# one parameter: $gi of type genomic interval
	method write_conservation_profiles( Conservation_Profiles $profiles ) {		
		my @gene_ids = ();
		my $gene_regions = {};
		$self->{XML} = {};
		foreach my $gi (@{$profiles->genomic_intervals}) {
			push @{$self->XML->{dna}[0]->{chromosomefragment}}, 
				{
				    'id' => $gi->label,
				    'region' => $gi->region,
					'strand' => $gi->strand,
			        'five_prime_pos' => $gi->five_prime_pos,
			        'three_prime_pos' => $gi->three_prime_pos,
			        'sequence' => [{ 'data' => $gi->gi_sequence }],
			        'masked_sequence' => [{ 'data' => $gi->get_repeatmasked_sequence() }],
					'source' => $gi->genome_db,
				};


			my $label = $gi->label;
			push @gene_ids, $label;
			$gene_regions->{$label} = {
				'strand' => $gi->strand,
		        'five_prime_pos' => $gi->five_prime_pos,
		        'three_prime_pos' => $gi->three_prime_pos,
			};
		}

		my %gene_profiles = ();
		
		while ( my ($key, $value) = each(%{$profiles->profile_pairs})) {
			# APPLES hashes things such that $one < $two
			my ($one, $two) = split /\_/, $key;
			
			my $id1 = $gene_ids[$one];
			my $id2 = $gene_ids[$two];
			
			my $profile1 = $value->profile1;
			my @profile1_xy = ();
			
			my $step1 = $value->step_width1;
			my $ofs = 0;
			my $winlen = $value->window_length;

			foreach (@{$profile1}) {
				push @profile1_xy, { 'pos' => $ofs, 'score' => $_ };
				$ofs+= $step1;
			}
			
			my $profile2 = $value->profile2;
			my @profile2_xy = ();
						
			my $step2 = $value->step_width2;

			$ofs = 0;
			foreach (@{$profile2}) {
				push @profile2_xy, { 'pos' => $ofs, 'score' => $_ };
				$ofs+= $step2;
			}

			push @{$gene_profiles{'conservationprofile'}},
			    {
						'gene' => $gene_ids[$one],
						'ortholog' => $gene_ids[$two],
						'window_length' => $winlen,
						'score' => \@profile1_xy,
				};

			push @{$gene_profiles{'conservationprofile'}},
				{
						'gene' => $gene_ids[$two],
						'ortholog' => $gene_ids[$one],
						'window_length' => $winlen,
						'score' => \@profile2_xy,
				};
		}
		
		$self->XML->{conservationprofiles}[0] = \%gene_profiles;
	}    # write_conservation_profiles

=pod
    
=head2 write_XML( Str $filename, Str $schema)

Write the content of $self->XML to the specified file. Then validate using specified schema (not implemented).

=cut

	method write_XML(Str $filename, Str $schema) {
		
	my $xs = new XML::Simple(NoAttr => 1);
    	open FILE, ">$filename"  or die $!;
    	print FILE $xs->XMLout($self->XML, RootName => 'applesdata', XMLDecl    => "<?xml version='1.0'?>");
	close FILE;
#		qx/gzip $filename/;
        
#        $validator = XML::Validator::Schema->new(file => $schema);
#
#        $parser = XML::SAX::ParserFactory->parser(Handler => $validator);
#
#        eval { $parser->parse_uri('$filename') };
#
#        if ($@) {
#        	#unlink($filename);
#            die "File failed validation: $@";
#        } else {
#        	print "Validated!\n";
#        }
		
	}

}    # XML_Writer
