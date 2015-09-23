#!/usr/bin/perl

=head1 Class to denominate non-serializable classes

Add this class to @ISA to prevent automatic serialization of a 
class or its subclasses.

=cut

package Serialization::Non_Serializable;

use strict;
use Moose;


1;