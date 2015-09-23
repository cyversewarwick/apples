#!/usr/bin/perl

package Datatypes::Moose_Types;

=head1 Moose Type Validation

In this package, we pre-declare datatypes when we want to be able
to

1) validate data using Moose type constraints

2) Get default values for these types (i.e. values we can use to initalize
   required members of this type)

3) a String-based suggestion mechanism, so we can suggest values via
   the Validation Web Service

This package is likely to grow quite large, so a FIXME is finding
a way to split it up.

=cut

use strict;

our @ISA = qw(Exporter);
require Exporter;
our @EXPORT_OK = qw(default_value suggest_values validate_value
  %datatypes_type_defaults %datatypes_type_proposals can_suggest
  %datatypes_input_class type_isa_class type_is_subtype_of
  type_is_arrayref type_is_hashref contained_type);
our @EXPORT = qw(default_value suggest_values validate_value can_suggest
  type_isa_class type_is_subtype_of type_is_arrayref type_is_hashref
  contained_type
);

use Configuration::AppleSeeds;
use Runtime;

use Moose;
use Moose::Util::TypeConstraints;

## This hash stores valid default initializers for each type constraint
our %datatypes_type_defaults = ();

## this has stores type proposals
## such proposals can either be an array of valid values
## or a ref to a sub taking a single parameter, returning a value
our %datatypes_type_proposals = ();

## If we want to use a special input class for reading a certain
## datatype, we can put this in here. See e.g. the DatasetName subtype
our %datatypes_input_class = ();

=head1 Datatype Declarations

Here we include all parameter types
Add your own parameters in the same way, note that they need to
be in _this_ package (Moose_Types)

Also have a look at ORANGES/javascript/oranges_validation.js,
there, you can define faster javascript validators.

=cut

=head2 All search Databases in Runtime

=cut

subtype 'Search_Database_ID' => as Str => where {
	my $val = $_;
	my $db  = Runtime::get_search_database($val);
	return defined($db) && UNIVERSAL::isa( $db, "Sequences::Database::Search" );
} => message {
	"Please choose a valid search database identifier.";
};

$datatypes_type_defaults{Search_Database_ID} =
  sub { return get_search_database_ids()->[0] };
$datatypes_type_proposals{Search_Database_ID} = sub {
	my $prefix = shift;
	return filter( get_search_database_ids(), $prefix );
};
$datatypes_input_class{Search_Database_ID} = 'Forms::Input::MultipleChoice';

=head2 All sequence Databases in Runtime

=cut

subtype 'Sequence_Database_ID' => as Str => where {
	my $val = $_;
	return 1 if (lc($val) eq "none");
	my $db  = Runtime::get_sequence_database($val);
	return defined($db)
	  && UNIVERSAL::isa( $db, "Sequences::Database::Sequence" );
} => message {
	"Please choose a valid sequence database identifier.";
};

$datatypes_type_defaults{Sequence_Database_ID} =
  sub { return get_sequence_database_ids()->[0] };

$datatypes_type_proposals{Sequence_Database_ID} = sub {
	my $prefix = shift;
	return filter( get_sequence_database_ids(), $prefix );
};

$datatypes_input_class{Sequence_Database_ID} = 'Forms::Input::MultipleChoice';

=head2 MIME Types for ResultDB

=cut

subtype 'MIMEType' => as Str => where {
	my $val = $_;
	return defined($val)
	  && (
		scalar( grep { $_ eq $val } @{ $datatypes_type_proposals{MIMEType} } )
	  );
} => message {
	"Please choose a supported MIME type.";
};

$datatypes_type_defaults{MIMEType} = sub { return "text/plain" };

$datatypes_type_proposals{MIMEType} = [
	"application/javascript",   "application/json",
	"application/pdf",          "application/zip",
	"application/x-gzip",       "application/octet-stream",
	"application/octet-stream", "image/svg+xml",
	"image/gif",                "image/jpeg",
	"image/png",                "image/tiff",
	"text/plain",               "text/css",
	"text/csv",                 "text/html",
	"text/xml",
];

$datatypes_input_class{MIMEType} = 'Forms::Input::MultipleChoice';

=head2 BLAST Search tool types

=cut

subtype 'BLASTTool' => as Str => where {
	my $val = $_;
	return defined($val)
	  && (
		scalar( grep { $_ eq $val } @{ $datatypes_type_proposals{BLASTTool} } )
	  );
} => message {
	"Please choose a supported BLAST-like tool type.";
};

$datatypes_type_defaults{BLASTTool} = sub { return "blastn" };

$datatypes_type_proposals{BLASTTool} = [
	"blastn",
	"tblastx",
];

$datatypes_input_class{BLASTTool} = 'Forms::Input::MultipleChoice';

=head2 Anchor types for sequence locations

=cut

subtype 'SequenceAnchor' => as Str => where {
	my $val = $_;
	return defined($val)
	  && (
		scalar( grep { $_ eq $val } @{ $datatypes_type_proposals{SequenceAnchor} } )
	  );
} => message {
	"Please choose a supported sequence anchor type.";
};

$datatypes_type_defaults{SequenceAnchor} = sub { return "5' end" };

$datatypes_type_proposals{SequenceAnchor} = [
	"5' end", "3' end"
];

$datatypes_input_class{SequenceAnchor} = 'Forms::Input::MultipleChoice';

=head2 Anchor types for sequence locations

=cut

subtype 'Sequence_Search' => as 'Str';

$datatypes_type_defaults{Sequence_Search} = sub { return "" };

$datatypes_input_class{Sequence_Search} = 'Forms::Input::Sequence_Search';

=head2 Empirical scoring score types

=cut

subtype 'PSSM_Scoringfunction' => as Str => where {
	my $val = $_;
	return defined($val)
	  && (
		scalar( grep { $_ eq $val } @{ $datatypes_type_proposals{PSSM_Scoringfunction} } )
	  );
} => message {
	"Please choose a supported pssm scoring function.";
};

$datatypes_type_defaults{PSSM_Scoringfunction} = sub { return "bio" };

$datatypes_type_proposals{PSSM_Scoringfunction} = [
	"bio",
	"mult none",
	"mult sqrt",
	"mult linear",
	"mult bifa",
];

$datatypes_input_class{PSSM_Scoringfunction} = 'Forms::Input::MultipleChoice';

=head2 Motif overrepresentation test types

=cut

subtype 'Overrep_Test' => as Str => where {
	my $val = $_;
	return defined($val)
	  && (
		scalar( grep { $_ eq $val } @{ $datatypes_type_proposals{Overrep_Test} } )
	  );
} => message {
	"Please choose a overrepresentation filter.";
};

$datatypes_type_defaults{Overrep_Test} = sub { return "none" };

$datatypes_type_proposals{Overrep_Test} = [
	"none",
	"binomial",
	"binomial_R",
	"poissonbinomial_R",
	"qvalue",
	"binomial_and_qvalue",
];

$datatypes_input_class{Overrep_Test} = 'Forms::Input::MultipleChoice';

=head2 ResultDB dataset names

We do not suggest names here, since we don't want to go through AJAX calls to
do this. Actually, we could, but this would take more time than it's worth.

Mainly, this type serves as a placeholder so that ORANGES forms can figure
out when to use a custom control for input.

=cut

subtype 'ResultDataset' => as 'Str';

$datatypes_input_class{ResultDataset} = 'Forms::Input::ResultDataset';

=head2 Default values for built in types

=cut

$datatypes_type_defaults{Int} = sub { return 0 };
$datatypes_type_defaults{Num} = sub { return 0 };
$datatypes_type_defaults{Str} = sub { return "" };

###############################################################################
###############################################################################
###############################################################################
###############################################################################

=head1 Implementation Part

=cut

=head2 Suggest a suitable value for a datatype, based on a string prefix.

 Parameters:
 $suggest_type : a type
 $prefix : part of the value string to base the suggestion on

 Returns:
 an ARRAYREF of valid values for initialising the given type

=cut

sub suggest_values {
	my $suggest_type = shift;
	my $prefix = shift || "";
	our %datatypes_type_proposals;

	if ( can_suggest($suggest_type) ) {
		my $suggestions = $datatypes_type_proposals{$suggest_type};

		if ( ref($suggestions) eq 'ARRAY' ) {
			return filter( $suggestions, $prefix );
		} elsif ( ref($suggestions) eq 'CODE' ) {
			return $suggestions->($prefix);
		} else {
			confess "Invalid suggestion info.";
		}

	} else {
		return undef;
	}
}

=head2 Suggest a suitable value for a datatype, based on a string prefix.

 Parameters:
 $suggest_type : a type

 Returns:
 True if we can suggest values for this type

=cut

sub can_suggest {
	my $suggest_type = shift;
	return
	  scalar( grep { $_ eq $suggest_type } ( keys %datatypes_type_proposals ) )
	  > 0;
}

=head2 Get default value for a datatype, based on a string prefix.

 Parameters:
 $suggest_type : a type

 Returns:
 a valid value for initialising the given type

=cut

sub default_value {
	my $suggest_type = shift;
	if ( grep { $_ eq $suggest_type } ( keys %datatypes_type_defaults ) ) {
		return $datatypes_type_defaults{$suggest_type};
	}
	if ( type_isa_class( $suggest_type, "Serialization::Serializable" ) ) {
		return sub { return $suggest_type->new; }
	}
	return undef;
}

=head2 Validate a value against a type

 Parameters:
 $suggest_type : a type
 $value : the value to test

 Returns:
 undef if the value was accepted, a message string otherwise

=cut

sub validate_value {
	my $suggest_type = shift;
	my $value        = shift;
	my $constraint   = find_type_constraint($suggest_type);
	if ( !defined($constraint) ) {
		confess("Unknown type.");
	}
	return $constraint->validate($value);
}

=head2 Filter to only include strings with a given prefix

 Parameters:
 $in_array: an ARRAYREF of strings
 $prefix : the prefix

 Returns:
 an ARRAYREF to all elements in $in_array that start with $prefix
=cut

sub filter {
	my $in_array = shift;
	my $prefix   = shift;

	if ( $prefix eq "" ) {
		return $in_array;
	}

	my @out = ();
	foreach my $v (@$in_array) {
		if ( $v =~ m/$prefix/i ) {
			push @out, $v;
		}
	}
	return \@out;
}

=head2 Check if a given type is subtype of Num

 Parameters:
 $type : a type identifier
 $subtype_of : a type to check whether $type is a subtype of

 Returns:
 1 if $type is subtype of $subtype_of
 0 otherwise

=cut

sub type_is_subtype_of {
	my $type       = shift;
	my $subtype_of = shift;

	if ( $type eq $subtype_of ) {
		return 1;
	}

	my $moosetype =
	  Moose::Util::TypeConstraints::find_or_create_isa_type_constraint($type);
	if ( !defined($moosetype) ) {
		return 0;
	}
	return $moosetype->is_subtype_of($subtype_of);
}

=head2 Check if a given type checks for a class

 Parameters:
 $type : a type identifier
 $class : the class to check for

 Returns:
 1 if $type checks for objects of type $class
 0 otherwise

=cut

sub type_isa_class {
	my $type  = shift;
	my $class = shift;
	my $moosetype =
	  Moose::Util::TypeConstraints::find_or_create_isa_type_constraint($type);

	if ( defined($moosetype)
		&& UNIVERSAL::isa( $moosetype, "Moose::Meta::TypeConstraint::Class" ) )
	{
		my $val = 0;
		eval {
			eval "require " . $moosetype->class;
			my $test = bless {}, $moosetype->class;
			$val = UNIVERSAL::isa( $test, $class );
		};
		return $val;
	}
	return 0;
}

=head2 Check if a given type is an arrayref

 Parameters:
 $type : a type identifier

 Returns:
 1 if $type checks for an arrayref
 0 otherwise

=cut

sub type_is_arrayref {
	my $type = shift;
	return $type =~ m/^ArrayRef/;
}

=head2 Get the type name contained in an array or hashref

 (This does not check if the type is actually an array or
  hashref -- so for simple types, it will return the type
  name itself)

 Parameters:
 $type : a type identifier

 Returns:
 the contained typename (or "")

=cut

sub contained_type {
	my $type = shift;
	if ( $type =~ m/^ArrayRef(\[.*\])?/ ) {
		$type = $1;
	}
	if ( $type =~ m/^HashRef(\[.*\])?/ ) {
		$type = $1;
	}
	$type =~ s/^\[(.*)\]$/$1/;
	return $type;
}

=head2 Check if a given type is a hash ref

...or serializable object...

 Parameters:
 $type : a type identifier
 $class : the class to check for

 Returns:
 1 if $type checks for an hash ref
 0 otherwise

=cut

sub type_is_hashref {
	my $type = shift;
	return $type =~ m/^HashRef/
	  || type_isa_class( $type, "Serialization::Serializable" );
}

1;
