#!/usr/bin/perl

=head1 JSON Web service client class

=cut

package Utilities::JSON_Client;

use strict;

use File::Spec qw(catfile);
use LWP::UserAgent;
use HTTP::Cookies;
use HTTP::Response;
use HTTP::Request::Common qw(POST);
use Time::HiRes qw(time);
use Data::Dumper;

use JSON;

use Configuration::AppleSeeds;

use Runtime;
use Carp;

=head2 Constructor

 Parameters:
 $class : Runtime::Remote_Job_Service
 $url   : The url where jobservice.pl is running (we will append jobservice.pl?)

 $huser : http auth user name (optional)
 $hpass : http auth password (optional)

 $retrycount : number of times to try (optional, default : 3)
 $timeout : receive timeout in s (optional, default : 180)

 Returns:

 a new JSON_Client

=cut

sub new {
	my ( $class, $url, $huser, $hpass, $retrycount, $timeout ) = @_;
	my $self = {
		url        => $url,
		huser      => $huser,
		hpass      => $hpass,
		retrycount => $retrycount || 3,
		timeout    => $timeout || 180,
	};

	my $cookies = {};

	if ( defined( get_config_key("usertempdir") )
		&& -d get_config_key("usertempdir") )
	{
		$cookies = HTTP::Cookies->new(
			file => File::Spec->catfile(
				get_config_key("usertempdir"),
				"lwpcookies.txt"
			),
			autosave => 1
		);
	}

	$self->{cookies} = $cookies;

	bless $self, $class;
	return $self;
}

=head2 Perform a JSON request to the server

This will return a HASHREF from

 Parameters:
 $self : a self object
 $urn  : the service name to access
 $parameters : the POST parameters

 Returns:
 a translated HASHREF from the JSON call

=cut

sub request {
	my $self = shift;
	my $urn  = shift;
	my $data = shift;
	
	my $error = "";

	my $uri = $self->{url};

	if ( $uri =~ m/[^\/]$/ ) {
		$uri .= "/";
	}
	$uri .= $urn;

	my $tries = $self->{retrycount};
	my $res   = undef;

	while ( $tries > 0 ) {
		eval {
			my $ua = LWP::UserAgent->new;
			## TODO this is a workaround for the fact that on www.wsbc
			## the SSL certificate does not validate
			$ua->ssl_opts ( verify_hostname => 0 ) ;
			$ua->cookie_jar($self->{cookies});
			$ua->timeout( $self->{timeout} );
			my $req = POST $uri,
			  Content_Type => 'form-data',
			  Content      => $data;

			if ( defined( $self->{huser} ) && defined( $self->{hpass} ) ) {
				$req->authorization_basic( $self->{huser}, $self->{hpass} );
			}

			my $t0 = time;
			my $response = $ua->request($req);
			my $t1 = time;
			debug ("JSON Request time: " . ($t1-$t0));

			if ( $response->is_success ) {
				my $type = $response->content_type();
				my $text = $response->decoded_content;
				if ( $text eq "") {
					die "Response is empty.";
				}

				if ( $type =~ m/json/ ) {
					eval { $res = from_json($text); };

					if ( !defined($res) || $@ ) {
						if ($@) {
							$error = $@;
						} else {
							$error = "Invalid JSON response.";
						}
					} else {
						## Here is where we return successfully (out of the eval)
						return;
					}
				} else {
					## everything non-json, we return like this:
					$res = {
						NONJSON =>
						  {
							type => $type,
							data => $text,
							size => length($text),
						  }
					};
				}

			} else {
				$error = $response->status_line;
				warn( "Retrying access to JSON Service $urn : " . $error );
			}
		};
		if ($@) {
			$error = $@;
		} elsif ( defined $res ) {
			## Here is where we return successfully
			return $res;
		}
		--$tries;
	}

	confess( "Could not access JSON Service $urn : " . $error . Dumper ($data) );
}

1;
