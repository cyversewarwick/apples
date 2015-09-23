### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Release_Notes_Maker Class ###
use MooseX::Declare;

class Release_Notes_Maker {
	use APPLES_Datatypes qw (WebService);
	use Release_Notes;
	use Version_Number;
	
	method make_release_notes () {
		# hard-code release notes in this class
		
		my @release_notes;
		
		my $release_note = Release_Notes_Triplet->new(
		    'version_number' => Version_Number->new(
			'major' => 1,
			'minor' => 0,
			'patch' => 0),
		    'release_note' => 'First release'
		    );
		push (@release_notes, $release_note);
		
		$release_note = Release_Notes_Triplet->new(
		    'version_number' => Version_Number->new(
			'major' => 1,
			'minor' => 0,
			'patch' => 1),
		    'web_service' => 'meme',
		    'release_note' => 'update to meme',
		    );
		push (@release_notes, $release_note);

		# template code to use for additional release note:
		
		#$release_note = Release_Notes_Triplet->new(
		#												'version_number' => Version_Number->new(
		#													'major' => ,
		#													'minor' => ,
		#													'patch' => ),
		#												'web_service' => '',
		#												'release_note' => '',
		#												);
		#push (@release_notes, $release_note);

		my $all_release_notes = Release_Notes->new('release_notes' => \@release_notes);
				
		return $all_release_notes;
	} # make_release_notes #

} # Release_Notes_Maker #
