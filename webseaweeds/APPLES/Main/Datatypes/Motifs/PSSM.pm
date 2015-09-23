#!/usr/bin/perl

use MooseX::Declare;

=head1 PSSM Motif class

Stores a PSSM and related information.

=cut

class Datatypes::Motifs::PSSM extends Links::Node {

	use Runtime;
	use Configuration::AppleSeeds;
	use File::Temp qw(tempfile);

	require Datatypes::Concepts::Transcription_Factor;

	## we have a name
	has 'name' => (
		is => "rw",
		isa => "Str",
		default => sub {return "BLANK"},
	);

	## e.g. Transfac or BiFa assign an accession
	has 'accession' => (
		is => "rw",
		isa => "Str",
		default => sub {return "BLANK"},
	);

	## PSSM matrix :
	has 'matrix' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
		default => sub {return []},
	);

	## background character distribution
	has 'background' => (
		is => 'rw',
		isa => 'ArrayRef[Num]',
		## A / C / G / T
		default => sub {return [ 0.3, 0.2, 0.2, 0.3 ]},
	);

	## PSSM number of locations
	has 'length' => (
		is => 'rw',
		isa => 'Int',
		default => sub {return 0;},
	);

	## PSSM shift offset
	has 'offset' => (
		is => 'rw',
		isa => 'Int',
		default => sub {return 0;},
	);

	## Metainformation from BiFa/Transfac
	has 'metainfo' => (
		is => 'rw',
		isa => 'Any',
		default => sub {return  {};},
	);

	## Pseudocount. See Boris' discussion in the APPLES documentation repo
	has 'pseudocount' => (
		is => 'rw',
		isa => 'Num',
		default => sub {return  0.25;},
	);

	## We might also know a pathway (index of some sort into Transpath)
	has 'pathway' => (
		is => 'rw',
		isa => 'Str',
		default => sub {return  "Unknown";},
	);

	## we can also have a list of associated transcription factors
	has 'factors' => (
		is => 'rw',
		isa => 'ArrayRef[Ref]',
		default => sub {return  [];},
	);

	## number of sites this motif came from
	has 'nsites' => (
		is => 'rw',
		isa => 'Int',
		default => sub {return  1;},
	);

=head2 Create a dot label

=cut
	method dot_style () {
		return
		    "label=\"PSSM: "
		  . $self->name
		  . "\" shape=\"trapezium\"";
	}

=head2 Set Probability for a character and location.

Allowed Characters are acgt.

 Parameters:
 $pos : the position
 $char : the character
 $prob : the probability

=cut
	method set_probability (Int $pos, Str $char, Num $prob) {
		$self->matrix->[($self->offset + $pos) * 4 + int ($char)] = $prob;
		$self->length ( int ( ( scalar @{$self->matrix} ) / 4 ) );
	}

=head2 Set Probability for a character and location.

Allowed Characters are acgt.

 Parameters:
 $pos : the position
 $char : the character

 Returns:
 the probability of getting character $char at $pos

=cut
	method get_probability (Int $pos, Str $char) {
		my $chr = $self->_map_char($char);
		return $self->matrix->[($self->offset + $pos) * 4 + $chr] || $self->background->[ $chr ] || 0.0;
	}


=head2 Add a related transcription factor

=cut
	method add_factor (Ref $factor) {
		if (!UNIVERSAL::isa ($factor, 'Datatypes::Concepts::Transcription_Factor')) {
			die "Can only add Datatypes::Concepts::Transcription_Factor objects.";
		}
		push @{$self->factors}, $factor;
	}

=head2 Map character to index

 Parameters:
 $char : a characters

 Returns:
 The matrix column index

=cut

	method _map_char (Str $char) {
		my %chrs = (
			'a' => 0,
			'c' => 1,
			'g' => 2,
			't' => 3,
		);
		my $rv = $chrs{ lc ($char) };
		if (!defined ($rv)) {
			die "Unsupported PSSM Character: $char" ;
		}
		return $rv;
	}

=head2 Helper to produce a dump of the matrix in a legible format

=cut

	method dump_probabilities () {
		my $rv = $self->name . " (" . $self->length . ") : \n";
		my %chrs = (
			'a' => 0,
			'c' => 1,
			'g' => 2,
			't' => 3,
		);
		$rv.= "$_\t" foreach (sort {$a cmp $b} (keys %chrs));
		$rv.= "\n";
		for (my $j = 0; $j < $self->length(); ++$j) {
			$rv.= ($self->get_probability($j, $_) . "\t") foreach  (sort {$a cmp $b} (keys %chrs));
			$rv.= "\n";
		}
		return $rv;
	}

=head2 dump in MEME format

=cut

	method dump_meme () {
		my %chrs = (
			'A' => 0,
			'C' => 1,
			'G' => 2,
			'T' => 3,
		);
		my $rv = <<END;
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from
END
		foreach my $chr (sort {$a cmp $b} (keys %chrs)) {
			$rv.= $chr .  " " . $self->background->[ $chrs{$chr} ] . " " ;
		}
		my $n = $self->name;
		my $l = $self->length;
		my $ns = $self->nsites;
		$rv.= <<END;

MOTIF $n
letter-probability matrix: alength= 4 w= $l nsites= $ns E= 0
END
		for (my $j = 0; $j < $self->length(); ++$j) {
			$rv.= " ";
			$rv.= ($self->get_probability($j, $_) . "\t") foreach  (sort {$a cmp $b} (keys %chrs));
			$rv.= "\n";
		}
		$rv.= "\n";
		return $rv;
	}
=head2 Get a PNG image for the motif

 Returns:
 binary data for a PNG file.

=cut

	method get_png (Int $char_w = 100, Int $h = 500) {
		my $cacheid = {
			what => "MOTIF_PNG_SIZE_${char_w}_$h",
			acc => $self->accession,
		};

		my $buf = cache->cache_get_hash ($cacheid);
		if (defined ($buf)) {
			return $buf;
		}

		## This stuff works when you have installed bioconductor, and their seqLogo package
		## check http://www.bioconductor.org/packages/2.3/bioc/html/seqLogo.html
		my ($fh, $filename) = tempfile();
		close ($fh);

		my $w = $char_w * $self->length + 20;

## uncomment to have the output from R available directly,
## by default this will go to a temp file which gets deleted afterwards
#		$filename = $self->accession . "_.png";

		my $info = $self->name . " : " . $self->pathway;
		$info =~ s/\"/\"\"/;

		my $matrix = join ",", @{$self->matrix};

		## this seems to leak memory like crazy,
		## which is why we run it in a separate job
		eval_R  ( <<END
library (seqLogo)
v = c( $matrix )
png(filename = "$filename", width = $w, height = $h,
	units = "px", pointsize = 12, bg = "white",
	res = NA)

mySeqLogo = seqLogo::seqLogo
bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") |
        sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] = NULL
norm = function(x) scale(x, center=FALSE, scale=colSums(x))

grid.newpage()
pwm = norm ( matrix (v, nrow=4, byrow=FALSE) )

pushViewport(viewport(x=0.5, y=0.5, width=1.0, height=1.0 ) )
mySeqLogo(pwm)
grid.text(sprintf("$info"), x=0.5, y=0.97,
          hjust=0.5, vjust=1, gp=gpar(fontsize=15, col="darkblue"))
popViewport()
dev.off()
END
		);
		$buf = read_file ( "$filename" );

		cache->cache_put_hash ($cacheid, $buf);

		## also put under name.
		$cacheid = {
			what => "MOTIF_PNG_SIZE_${char_w}_$h",
			acc => $self->name,
		};

		cache->cache_put_hash ($cacheid, $buf);

		return $buf;

		## Alternatives

		## use PS intermediate stage
		#$R->run(qq`postscript("$filename" , horizontal=FALSE , width=500 , height=500 , pointsize=1)`);
		#my $convert = make_executable (trim (`which convert`) || 'convert');
		#system $convert, "$filename", "$filename.png";

		# this seems to cause an error.

		# png(filename = "$filename", width = $w, height = $h,
		# 	units = "px", pointsize = 12, bg = "white",
		# 	res = NA, type = c("cairo", "cairo-png", "Xlib", "quartz"))
	}

=head2 Unique ID override

=cut

	method unique_id () {
		return blessed($self) . " : $self->{accession}";
	}

}
