# This class contains modified code taken from:
# Lenhard B., Wasserman W.W. (2002) TFBS: Computational framework for transcription factor binding site analysis. Bioinformatics 18:1135-113
#

use MooseX::Declare;

class Graphics {
	
	use General_Utilities;
	use Data::Dumper;
	use GD;
	use Statistics::R; 
	use APPLES_Datatypes qw (PositiveNum);
		
	my $GU = General_Utilities->new();

	method create_polar_plot_for_pvalues(ArrayRef[ArrayRef] $p_values, ArrayRef $pssm_names, Str $output_file_name, Str $plot_title, PositiveNum $minus_log_threshold) {

	#------------------------
	# Check input - no need to continue if something is wrong
	# 1. The first dimension in p_values should equal the number of entries in pssm_names.
	# 2. The maximum number of p-value arrays that can be plotted in one figure is 8 (colour limitation).
	# 3. The number of p-values must be the same for each pssm.
	#------------------------ 
		
	my $n_exp = @{$p_values};
	my $n_pssm = @{$pssm_names}; 
	if ($n_exp != $n_pssm){	
		die "The number of PSSM-names is not consistent with size of the number of p-value arrays.\n";
	} elsif ($n_exp > 8){
		die "Cannot create plot: more pssms than colours available for plotting (8).\n"		
	}
	
	my @n_p_values;
	foreach my $exp (@{$p_values}){
		my $n_p_values = @{$exp};
		push(@n_p_values,$n_p_values);
	}
	@n_p_values = sort(@n_p_values);
	if ($n_p_values[0] != $n_p_values[$n_exp-1]){
		die "The p-value arrays have different lengths.\n";
	}
										
#	------------------------
#	 Get path information from APPLES.dat 
#	------------------------ 
	my $APPLES_conf = new Config::General($ENV{'APPLES_DAT'});
	my %APPLES_config = $APPLES_conf->getall();
	my $path_to_heRMisc_R_package = $APPLES_config{path_to_heRMisc_R_package};
	my $path_to_APPLES_R_programs = $APPLES_config{path_to_APPLES_R_programs};
	

#	------------------------
#	 Start the Perl-R-interface 
#	------------------------
	my $R = Statistics::R->new() ;
	$R->startR ;

#	------------------------
#	 Send the parameters/variables
#	------------------------
	$R -> send(qq /filename <- \"/.$output_file_name.qq /\"/ );
	$R -> send(qq /path_to_heRMisc_R_package <- \"/.$path_to_heRMisc_R_package.qq /\"/ );
	$R -> send(qq /plot_title <- \"/.$plot_title.qq/\"/);
	$R -> send(qq /threshold <- /.$minus_log_threshold);
	
	$R -> send(qq /legendnames = c()/);
	foreach my $pssm_name (@{$pssm_names}){
		$R -> send(qq /legendname <- \"/.$pssm_name.qq/\"/);
		$R -> send(qq /legendnames <- c(legendnames,legendname)/);				}
		
#	for debugging (i.e. see what is sent to R)
#	$R->send(qq /print(parametername)/);
#	my $R_return = $R->read;
#	$GU->user_info(3,"\nReturned from R $R_return\n\n");		
	
	$R -> send(qq/ncol <- /.$n_p_values[0]); # number of p-values in each array
	$R -> send(qq/nrow <- /.$n_exp); # number of p-value arrays
	$R -> send(qq /r_mat <- matrix(0,nrow,ncol)/);
	my $i = 0;	
	foreach my $exp (@{$p_values}){		
		$i = $i+1;
		$R->send(qq /i <- /.$i);
		my $j = 0;
		foreach my $p_value (@{$exp}){
			$j = $j+1;
			if (($p_value > 0)&&($p_value <=1)) {
				$p_value = -log($p_value)/log(10);
			} elsif ($p_value == 0) {
				$GU->user_info(3,"\nP value of 0 found, replacing negative logged p_value with 20\n");
				$p_value = 20;
			} else {
				die "One of your p-values in the input to the polar plot function is not between 0 and 1";
			}
			$R->send(qq /j <- /.$j);
			$R->send(qq /r_mat[i,j] <- /.$p_value);	
			
		} 
	}

#	------------------------
#	 Call the script that does the actual plotting
#	------------------------
	$R->send(qq /source("/.$path_to_APPLES_R_programs.qq/\/Create_Polar_Plot.R")/);

#	------------------------
#	 Stop the Perl-R-interface 
#	------------------------
	$R->stopR() ;

	} # create_polar_plot_for_pvalues #

	
	method generate_PSSM_logo(ArrayRef $info_content_matrix, Str $pssm_id, Str $outdir){
		
		my @icm_matrix = @{$info_content_matrix};
		
		my @transposed_icm_matrix;
		
		for (my $i=0; $i < 4; $i++){
			
			for (my $j=0; $j < scalar(@icm_matrix); $j++){
				
				$transposed_icm_matrix[$i][$j] = $icm_matrix[$j][$i];
			}
		}
		
		my $logo_filepath = $outdir."/$pssm_id.png";
		
		draw_logo(\@transposed_icm_matrix, -file=> $logo_filepath, 
						-full_scale =>2.25,
						-xsize=>500,
						-ysize =>250,  
						-x_title=>"position", 
						-y_title=>"bits");
		
		return $logo_filepath;
		
		sub draw_logo {
			
			#no strict;
			my $matrixref = shift;
			my %args = (-xsize      => 600,
			-full_scale => 2.25,
			-graph_title=> "",
			-x_title    => "",
			-y_title    => "",
			@_);
			# Other parameters that can be specified:
			#       -ysize -line_width -margin
			# do not have a fixed default value 
			#   - they are calculated from xsize if not specified
			
			# draw postscript logo if asked for
			#if ($args{'-ps'} || $args{'-pdf'}){
			# return _draw_ps_logo($self, @_);  
			#}
			
			#require GD;
			
			
			my ($xsize,$FULL_SCALE, $x_title, $y_title)  = @args{qw(-xsize -full_scale -x_title y_title)} ;
			
			my $PER_PIXEL_LINE = 300;
			
			# calculate other parameters if not specified
			
			my $line_width = ($args{-line_width} or int ($xsize/$PER_PIXEL_LINE) or 1);
			my $ysize      = ($args{-ysize} or $xsize/1.6); 
			# remark (the line above): 1.6 is a standard screen x:y ratio
			my $margin     = ($args{-margin} or $ysize*0.15);
			
			my $image = GD::Image->new($xsize, $ysize);
			my $white = $image->colorAllocate(255,255,255);
			my $black = $image->colorAllocate(0,0,0);
			#my $motif_size = $self->pdl_matrix->getdim(0);
			my $motif_size = scalar(@{$$matrixref[0]});
			my $font = ((&GD::gdTinyFont(), &GD::gdSmallFont(), &GD::gdMediumBoldFont(), 
			&GD::gdLargeFont(), &GD::gdGiantFont())[int(($ysize-50)/100)]
			or &GD::gdGiantFont());
			my $title_font = ((&GD::gdSmallFont(), &GD::gdMediumBoldFont(), 
			&GD::gdLargeFont(), &GD::gdGiantFont())[int(($ysize-50)/100)]
			or &GD::gdGiantFont());
			
			
			# WRITE LABELS AND TITLE
			
			# graph title   #&GD::Font::MediumBold
			$image->string($title_font,
			$xsize/2-length($args{-graph_title})* $title_font->width() /2,
			$margin/2 - $title_font->height()/2,
			$args{-graph_title}, $black);
			
			# x_title
			$image->string($font,
			$xsize/2-length($args{-x_title})*$font->width()/2,
			$ysize-( $margin - $font->height()*0  - 5*$line_width)/2 
			- $font->height()/2*0,
			$args{-x_title}, 
			$black);
			# y_title
			$image->stringUp($font,
			($margin -$font->width()- 5*$line_width)/2 
			- $font->height()/2 ,
			$ysize/2+length($args{'-y_title'})*$font->width()/2,
			$args{'-y_title'}, $black);
			
			
			# DRAW AXES
			
			# vertical: (top left to bottom right)
			$image->filledRectangle($margin-$line_width, $margin-$line_width, 
			$margin-1, $ysize-$margin+$line_width, 
			$black);
			# horizontal: (ditto)
			$image->filledRectangle($margin-$line_width, $ysize-$margin+1, 
			$xsize-$margin+$line_width,$ysize-$margin+$line_width,
			$black);
			
			# DRAW VERTICAL TICKS AND LABELS
			
			# vertical axis (IC 1 and 2) 
			my $ic_1 = ($ysize - 2* $margin) / $FULL_SCALE;
			foreach my $i (1..$FULL_SCALE)  {
				$image->filledRectangle($margin-3*$line_width, 
				$ysize-$margin - $i*$ic_1, 
				$margin-1, 
				$ysize-$margin+$line_width - $i*$ic_1, 
				$black);
				$image->string($font, 
				$margin-5*$line_width - $font->width,
				$ysize - $margin - $i*$ic_1 - $font->height()/2,
				$i,
				$black);
			}
			
			# DRAW HORIZONTAL TICKS AND LABELS, AND THE LOGO ITSELF 
			
			# define function refs as hash elements
			my %draw_letter = ( A => \&draw_A,
			C => \&draw_C,
			G => \&draw_G,
			T => \&draw_T );
			
			my $horiz_step = ($xsize -2*$margin) / $motif_size;		

			foreach my $x (0..$motif_size)  {
				
			    if ($x>70) {
				my $GU = General_Utilities->new();
				$GU->user_info(1,"cut motif length down to avoid crash (".$args{-file}.")\n");
				   # Kashi, Sascha: trying to print a motif-logo for R02192 (TRANSFAC) caused segmentation
				   #                faults. The error occurred in the draw_G-function when it tries to fill
				   #                the arc of the G at one of the last positions in this consensus sequence.
                                   #                Resolution was not obvious, therefore restricted this code. No need for
                                   #                logos this long at this point in time.
				last;
			    }

				$image->filledRectangle($margin + $x*$horiz_step, 
				$ysize-$margin+1, 
				$margin + $x*$horiz_step+ $line_width, 
				$ysize-$margin+3*$line_width, 
				$black);
				last if $x==$motif_size;								
				
				# get the $i-th column of matrix
				my %ic; 
				#($ic{A}, $ic{C}, $ic{G}, $ic{T}) = list $self->pdl_matrix->slice($x);
				$ic{A} = $$matrixref[0][$x];
				$ic{C} = $$matrixref[1][$x];
				$ic{G} = $$matrixref[2][$x];
				$ic{T} = $$matrixref[3][$x];
				#$GU->user_info(3,"->$x:\n");
				#$GU->user_info(3,$$matrixref[$x][0]."\n");
				#$GU->user_info(3,$$matrixref[$x][1]."\n");
				#$GU->user_info(3,$$matrixref[$x][2]."\n");
				#$GU->user_info(3,$$matrixref[$x][3]."\n");
				# sort nucleotides by increasing information content
				my @draw_order = sort {$ic{$a}<=>$ic{$b}} qw(A C G T);
				#$GU->user_info(3,Dumper (%ic));
				# draw logo column
				my $xlettersize = $horiz_step /1.1;
				my $ybottom = $ysize - $margin;		
			
				foreach my $base (@draw_order)  {
					#$GU->user_info(3,$base."\n");
					my $ylettersize = int($ic{$base}*$ic_1 +0.5);
					#$GU->user_info(3,$ylettersize."\n");
					next if $ylettersize ==0;
					
	    			# draw letter				
					$draw_letter{$base}->($image,
							      $margin + $x*$horiz_step,
							      $ybottom - $ylettersize,
							      $xlettersize, $ylettersize, $white);
					$ybottom = $ybottom - $ylettersize-1;
						    			
				}			
				
				if ($args{'-error_bars'} and ref($args{'-error_bars'}) eq "ARRAY")  {
					my $sd_pix   = int($args{'-error_bars'}->[$x]*$ic_1);
					my $yt     = $ybottom - $sd_pix+1;
					my $yb  = $ybottom + $sd_pix-1;
					my $xpos     = $margin + ($x+0.45)*$horiz_step;
					my $half_width;
					
					if ($yb > $ysize-$margin+$line_width)  {
						$yb = $ysize-$margin+$line_width
					}
					else {
						$image->line($xpos - $xlettersize/8, $yb, 
						$xpos + $xlettersize/8, $yb, 
						$black);
					}
					
					$image->line($xpos, $yt, $xpos, $yb, $black);
					$image->line($xpos - 1 , $ybottom, $xpos+1, $ybottom, $black);
					$image->line($xpos - $xlettersize/8, $yt, 
					$xpos + $xlettersize/8, $yt, 
					$black);
					
					
				}
				
				# print position number on x axis
				$image->string($font,
				$margin + ($x+0.5)*$horiz_step - $font->width()/2,
				$ysize - $margin +5*$line_width,
				$x+1,
				$black);
			}
			
			# print $args{-file};
			if  ($args{-file}) {  
				open (PNGFILE, ">".$args{-file})
				or die;
				print PNGFILE $image->png;
				close PNGFILE;
			}
			return $image;
		}
		
		
		
		
		
		
		sub total_ic  {
			return $_[0]->pdl_matrix->sum();
		}
		
		#################################################################
		# INTERNAL METHODS
		#################################################################
		
		
		sub _check_ic_validity  {
			my ($self) = @_;
			# to do
		}
		
#		sub DESTROY  {
			# nothing
#		}
		#### ------------------------------------------------------------
		sub draw_C  {
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $blue = $im->colorAllocate(0,0,204);
			$im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$blue);
			$im->fill($x+$xsize/2, $y+$ysize/2, $blue);
			if ($ysize>12) {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				0.55*$xsize, (0.625-0.625/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize/2, $y+$ysize/2.8+1, 
				$x+$xsize*1.1, $y+(2.8*$ysize/4)-1,
				$white);
			}
			elsif ($ysize>3)  {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				(0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
				$x+$xsize*1.1, $y+$ysize/2,
				$white);
				
			}
			return 1;
		}
		
		sub draw_G  {
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $yellow = $im->colorAllocate(255,179,0);
			$im->arc($x+$xsize*0.54, $y+$ysize/2,1.08*$xsize,$ysize,0,360,$yellow);
			$im->fill($x+$xsize/2, $y+$ysize/2, $yellow);
			if ($ysize>20) {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				0.55*$xsize, (0.625-0.625/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize/2, $y+$ysize/2.8+1, 
				$x+$xsize*1.1, $y+$ysize/2-1,
				$white);
			}
			elsif($ysize>3)  {
				$im->arc($x+$xsize*0.53, $y+$ysize/2, 
				(0.75-0.75/$ysize)*$xsize, (0.725-0.725/$ysize)*$ysize,
				0,360,$white);
				$im->fill($x+$xsize/2, $y+$ysize/2, $white);
				$im->filledRectangle($x+$xsize*0.25, $y+$ysize/2, 
				$x+$xsize*1.1, $y+$ysize/2,
				$white);
				
			}
			$im->filledRectangle($x+0.85*$xsize, $y+$ysize/2,
			$x+$xsize,$y+(3*$ysize/4)-1,
			$yellow);
			$im->filledRectangle($x+0.6*$xsize, $y+$ysize/2,
			$x+$xsize,$y+(5*$ysize/8)-1,
			$yellow);
			return 1;
		}
		
		sub draw_A {
			
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $green = $im->colorAllocate(0,204,0);
			my $outPoly = GD::Polygon->new();
			$outPoly->addPt($x, $y+$ysize);
			$outPoly->addPt($x+$xsize*.32, $y);
			$outPoly->addPt($x+$xsize*.68, $y);
			$outPoly->addPt($x+$xsize, $y+$ysize);
			$outPoly->addPt($x+0.75*$xsize, $y+$ysize);
			$outPoly->addPt($x+0.635*$xsize, $y+0.75*$ysize);
			$outPoly->addPt($x+0.355*$xsize, $y+0.75*$ysize);
			$outPoly->addPt($x+0.25*$xsize, $y+$ysize);
			$im->filledPolygon($outPoly, $green);
			if ($ysize>8)  {
				my $inPoly = GD::Polygon->new();
				$inPoly->addPt($x+$xsize*.5, $y+0.2*$ysize);
				$inPoly->addPt($x+$xsize*.40, $y+0.6*$ysize-1);
				$inPoly->addPt($x+$xsize*.60, $y+0.6*$ysize-1);
				$im->filledPolygon($inPoly, $white);
			}
			return 1;
		}
		
		sub draw_T {
			
			my ($im, $x, $y, $xsize, $ysize, $white) = @_;
			my $red = $im->colorAllocate(204,0,0);
			$im->filledRectangle($x, $y, $x+$xsize, $y+0.16*$ysize, $red);
			$im->filledRectangle($x+0.35*$xsize, $y, $x+0.65*$xsize, $y+$ysize, $red);
			return 1;
		}
		
	} # generate_logo #
	
} # Graphics #
