#!/usr/bin/perl

use strict;

my $good;
open(IN, "$ARGV[1]" );
while(<IN>){
	my ($gid) = split;	
	$good->{ $gid } = 1;
}
close(IN);

open(IN, "$ARGV[0]" );
while(<IN>){
        my @A = split;
	my $att = $A[8];
	my ($label, $id) = split( /=/, $att );
	print $_ if( exists $good->{$id});	
}
close(IN);




