### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Character_Matrix_Maker Class ###

use MooseX::Declare;

class Bundler::Dependencies::Character_Matrix_Maker {
    use Bundler::Dependencies::Character_Matrix;
    use Bundler::Dependencies::APPLES_Datatypes qw(NonNegativeInt);

    method make_character_matrix(Int $width, Int $height) {
	my @character_matrix;
	for (my $j=0;$j<$height;$j++) {
	    my @array = $self->private_make_array_of_spaces($width);
	    push(@character_matrix,\@array);
	}
	my $result = Bundler::Dependencies::Character_Matrix->new(character_matrix => \@character_matrix, width => $width, height => $height);
	return $result;
    } # make_character_matrix #

    method private_make_array_of_spaces(Int $number_of_spaces) {
	my @result;
	for (my $i=0;$i<$number_of_spaces;$i++) {
	    push(@result,' ');
	}
	return @result;
    } # private_make_array_of_spaces #

} # Character_Matrix_Maker #
