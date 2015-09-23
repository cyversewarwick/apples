### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Character_Matrix Class ###
# used for printing evolutionary trees, may be useful for generating other complex text output

use MooseX::Declare;

class Bundler::Dependencies::Character_Matrix {
    use Bundler::Dependencies::APPLES_Datatypes qw(NonNegativeInt);
    use Bundler::Dependencies::General_Utilities;

    has 'character_matrix' => (is => 'rw', isa => 'ArrayRef', required => 1);
    has 'width' => (is => 'rw', isa => 'NonNegativeInt', required => 1);
    has 'height' => (is => 'rw', isa => 'NonNegativeInt', required => 1);

    my $GU = Bundler::Dependencies::General_Utilities->new();

    method insert_line(Str $string, Int $line_index) {
	# writes a string into the centre of a line
	# line number count starts at zero

	my $length = length($string);
	if ($line_index<0) {
	    die 'invalid line index.';
	}
	if ($length > $self->width) {
	    die 'string is too long to insert into character matrix.';
	}
	if ($line_index>=$self->height) {
	    die 'this line does not exist in character matrix.';
	}
	my $remaining_spaces = $self->width-$length;
	my $offset = int($remaining_spaces/2);
	my @line = @{${$self->character_matrix}[$line_index]};
	for (my $i=0;$i<$length;$i++) {
	    $line[$offset+$i] = substr($string,$i,1);
	}
	${$self->character_matrix}[$line_index] = \@line;
    } # insert_line #

    method insert_matrix(Bundler::Dependencies::Character_Matrix $character_matrix, Int $upper_left_column, Int $upper_left_row) {
	# inserts another character matrix into $self 
	# $upper_left_column is the index of the column where the upper left corner of $character_matrix
	# will be placed, $upper_left_row accordingly

	my $last_column = $upper_left_column+($character_matrix->width-1);
	if ($last_column>=$self->width) {
	    die 'matrix too wide to fit into this one.';
	}
	my $last_line = $upper_left_row+($character_matrix->height-1);
	if ($last_line>=$self->height) {
	    die 'matrix too high to fit into this one.';
	}
	if (($upper_left_column<0)||($upper_left_row<0)) {
	    die 'invalid index.';
	}
	for (my $i=0;$i<$character_matrix->width;$i++) {
	    for (my $j=0;$j<$character_matrix->height;$j++) {
		my $column = $upper_left_column+$i;
		my $row = $upper_left_row+$j;
		${${$self->character_matrix}[$row]}[$column] = ${${$character_matrix->character_matrix}[$j]}[$i];
	    }
	}
    } # insert_matrix #

    method render_text() {
	for (my $i=0;$i<$self->height;$i++) {
	    my @line = @{${$self->character_matrix}[$i]};
	    for (my $j=0;$j<$self->width;$j++) {
		$GU->user_info(1,$line[$j]);
	    }
	    $GU->user_info(1,"\n");
	}
    } # render_text #

} # Character_Matrix #
