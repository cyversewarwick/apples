#!/usr/bin/perl

=head1 Class to render a links network as a dot file

=cut

package Links::Output::Dot;

use strict;

use Runtime;

use Scalar::Util qw(blessed);

=head2 Render a dot file of a Links_Database

 Parameters:
 $ldb : the database
 $wattr : (optional) the name of the weight attribute
 			default : "weight"
 
 Returns:
 the dot code as a string

=cut

sub render {
	my $ldb = shift;
	my $wattr = shift || "weight";

	my $dot = "digraph G {\n";

	my @nodes = @{ $ldb->get_nodes };

	my %bidir_linklist = ();

	## first draw the links
	foreach my $l ( @{ $ldb->get_links } ) {	
		my $dot_label1 = "N_" . $l->source->unique_id;
		my $dot_label2 = "N_" . $l->target->unique_id;

		if ( $l->bidirectional ) {
			if ( exists( $bidir_linklist{ $dot_label1 . $dot_label2 } ) ) {
				next;
			} else {
				$bidir_linklist{ $dot_label1 . $dot_label2 } = 1;
			}
		}
		$dot_label1 =~ s/[^A-Za-z0-9]/_/g;
		$dot_label2 =~ s/[^A-Za-z0-9]/_/g;

		my $linkstyle = (
			( $l->bidirectional )
			? " penwidth=3 dir=both fontcolor=red color=red arrowhead=\"odot\" arrowtail=\"odot\" "
			: "color=grey fontcolor=grey "
		);

		my $weight = $l->data($wattr);

		if ( defined($weight) ) {
			if ( ref($weight) eq "SCALAR" ) {
				$linkstyle .= " weight=\"$weight\"";
			} elsif ( UNIVERSAL::can( $weight, 'weight' ) ) {
				my $ww = $weight->weight;
				$linkstyle .= " weight=\"$ww\"";
			}
		}

		my $label        = "";
		my $simple_attrs = 0;
		my $complex_attrs = 0;
		foreach my $attr ( $l->attributes ) {
			my $attr_val = $l->data($attr);
			if ( UNIVERSAL::isa( $attr_val, "Links::Node" ) ) {
				## to keep things simple, if we can get a weight
				## from the link attribute node, we just draw
				## it as one edge
				if ( UNIVERSAL::can( $attr_val, 'weight' ) ) {
					my $wstr = sprintf "%1.1g", $attr_val->weight;
					$label .= "$attr : $wstr\\n";
					++$simple_attrs;
				} else {
					push @nodes, $attr_val;

					my $dot_label3 = "N_" . $attr_val->unique_id;
					$dot_label3 =~ s/[^A-Za-z0-9]/_/g;
					$dot .= "$dot_label1 -> $dot_label3 [ $linkstyle ] ; \n";
					$dot .= "$dot_label3 -> $dot_label2 [ $linkstyle ] ; \n";
					++$complex_attrs;
				}
			} else {
				$label .= "$attr : $attr_val\\n";
				++$simple_attrs;
			}
		}
		$linkstyle .= " label=\"$label\"";

		$dot .= "$dot_label1 -> $dot_label2 [ $linkstyle ] ;\n"
		  if ( $simple_attrs > 0 || $complex_attrs == 0 );
	}

	foreach my $n (@nodes) {
		my $dot_label = "N_" . $n->unique_id;
		$dot_label =~ s/[^A-Za-z0-9]/_/g;

		my $label     = $n->label;
		my $nodestyle = "label=\"$label\"";

		if ( defined( $n->dot_style ) ) {
			$nodestyle = $n->dot_style;
		}

		$dot .= "$dot_label [ $nodestyle ] ;\n";
	}
	$dot .= "}\n";
	return $dot;
}

1;
