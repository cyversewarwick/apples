#!/usr/bin/perl

=head1 Class to render a links network as a JSON file for arbor.js

=cut

package Links::Output::JSON;

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

	## use cache so we don't do this many times
	my $links = {
		nodes => {},
		edges => {},
	};

	foreach my $n ( @{ $ldb->get_nodes() } ) {
		$links->{nodes}->{ $n->unique_id } = _make_node($n);
	}

	foreach my $e ( @{ $ldb->get_links() } ) {
		my $v = {
			label         => "",
			seed          => int( rand(10) ),
			bidirectional => $e->bidirectional,
		};

		my $label        = "";
		my $simple_attrs = 0;
		my $weight       = 0;

		foreach my $attr ( $e->attributes ) {
			my $attr_val = $e->data($attr);
			if ( UNIVERSAL::isa( $attr_val, "Links::Node" ) ) {

				$v->{data}->{$attr} =
				  Serialization::Serializable::to_hash($attr_val);

				## to keep things simple, if we can get a weight
				## from the link attribute node, we just draw
				## it as one edge
				if ( UNIVERSAL::can( $attr_val, 'weight' ) ) {
					my $wstr = sprintf "%1.1g", $attr_val->weight;
					$label .= "$attr : $wstr\\n";

					## geom mean of all weights
					my $w = $attr_val->weight;
					$weight = sqrt( $weight * $weight + $w * $w );

					++$simple_attrs;
				} else {
					$label .= "$attr : " . blessed($attr_val) . "\\n";
				}
			} else {
				$v->{data}->{$attr} = $attr_val;
				$label .= "$attr : $attr_val\\n";
				++$simple_attrs;
			}
		}

		$v->{label} = $label;
		$v->{'length'} = $weight == 0 ? 1 : log($weight);

		if ( $v->{label} eq "" ) {
			$v->{label} = $e->label;
		}

		if ( $e->bidirectional ) {
			$links->{edges}->{ $e->source->unique_id }
			  ->{ $e->target->unique_id } = $v;
		} else {
			$v = {
				source => { name => $e->source->unique_id, },
				target => { name => $e->target->unique_id, },
				data   => $v,
			};
			push
			  @{ $links->{nodes}->{ $e->source->unique_id }->{unidir_links} },
			  $v;
		}

	}

	my $sum = 0;
	my $min = undef;
	my $max = undef;

## normalize edge lengths
	foreach my $k1 ( keys %{ $links->{edges} } ) {
		foreach my $k2 ( keys %{ $links->{edges}->{$k1} } ) {
			my $l = $links->{edges}->{$k1}->{$k2}->{'length'};
			next unless defined($l);
			$sum += $l;
			if ( !defined($min) || $l < $min ) {
				$min = $l;
			}
			if ( !defined($max) || $l > $max ) {
				$max = $l;
			}
		}
	}

	foreach my $k1 ( keys %{ $links->{edges} } ) {
		foreach my $k2 ( keys %{ $links->{edges}->{$k1} } ) {
			my $l = $links->{edges}->{$k1}->{$k2}->{length};
			next unless defined($l);

			if ( $max == $min ) {
				$l = 1;
			} else {
				$l = ( $l - $min ) / ( $max - $min );
			}

			if ( $l < 0 ) {
				$l = 1;
			}
			$l *= 10;

			$links->{edges}->{$k1}->{$k2}->{'length'} = $l;
		}
	}
	return $links;
}

=head2 Make a JSON node for the renderer from a links node

 Parameters:
 $n : a Links::Node

 Returns:
 a HASHREF for the renderer

=cut

sub _make_node {
	my $n = shift;

	my $dotstyle = $n->dot_style;
	my $dotlabel = $n->label;
	my $dotshape = "rect";

	if ( $dotstyle =~ m/label\s*\=\s*['"]?([^'"]+)/ ) {
		$dotlabel = $1;
	}

	if ( $dotstyle =~ m/shape\s*\=\s*['"]?([^\s'"]+)/ ) {
		$dotshape = $1;
	}

	#	my $mass = $n->out_degree() * 100;
	#	if ( $mass <= 1 ) {
	#		$mass = 1;
	#	}

	return {
		type => blessed($n),

		#		mass      => $mass,
		mass         => 2,
		label        => $n->label,
		dot_style    => $n->dot_style,
		dot_shape    => $dotshape,
		dot_label    => $dotlabel,
		color        => $n->color,
		data         => Serialization::Serializable::to_hash($n),
		unidir_links => [],
	};
}

1;
