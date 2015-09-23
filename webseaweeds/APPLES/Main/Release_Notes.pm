### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Release_Notes Class ###
# Release notes are only to be used for the web interface
use MooseX::Declare;

class Release_Notes {
	
	use Data::Dumper;
	use APPLES_Datatypes qw (WebService);
	
	has 'release_notes' => (is => 'rw', isa => 'ArrayRef[Release_Notes_Triplet]');
		
	method print_webservice_release_notes (WebService $webservice) {
		# gives all version numbers and associated release notes for a specified webservice
		print "*** Release Notes for ".$webservice." ***\n";
		foreach my $release_note (@{$self->release_notes}) {
			unless (!$release_note->web_service) {
				if ($release_note->web_service eq $webservice) {
					print "Version: ".$release_note->version_number->major.".".$release_note->version_number->minor.".".$release_note->version_number->patch."\tRelease note: ".$release_note->release_note."\n";
				}
			}
		}
	} # print_webservice_release_notes #
	
	method print_all_release_notes {
		# prints all release notes to date
		print "*** Release Notes ***\n";
		foreach my $release_note (@{$self->release_notes}) {
			print "release note: ".$release_note->release_note."\n";
		}

	} # print_all_release_notes #
	
} # Release_Notes #

class Release_Notes_Triplet {
	use APPLES_Datatypes qw (WebService);

	has 'version_number' => (is => 'rw', isa => 'Version_Number', required => 1);
	has 'release_note' => (is => 'rw', isa => 'Str', required => 1);
	has 'web_service' => (is => 'rw', isa => WebService);

} # Release_Notes_Triplet #
