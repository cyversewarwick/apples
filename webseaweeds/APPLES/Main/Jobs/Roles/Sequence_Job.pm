#!/usr/bin/perl

=head1 A sequence job works on all sequences which are available in a dataset

In this role, we have functionality for getting specific bits of sequence out
of the dataset.

=cut

package Jobs::Roles::Links_Job;

use feature ':5.10';

use Runtime;

use Moose::Role;
use MooseX::Declare;

use Datatypes::Moose_Types;
use Digest::SHA;
use JSON;


1;