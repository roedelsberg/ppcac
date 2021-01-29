#!/usr/bin/perl

use strict;

  

die "program blast1 blast

This script takes two BLASTP output files (generated -outfmt 6 option)
as input and outputs best reciprocal hits (BRH) and unidirectional
hits. 
\n" if(@ARGV != 2 );
 
#my $anno = read_anno( $ARGV[2] );

my $ab = read_file( shift @ARGV  );
my $ba =  read_file( shift @ARGV  );
foreach my $query( keys  %$ab ){
    my $t = $ab->{$query};
    if( exists $ba->{$t} ){
	if( $ba->{$t}  eq $query){
	    print  "$query\t$ab->{$query}\tBRH\n"  ;
	}
	else{
	    print "$query\t$ab->{$query}\tbest\n";
	}
    }
    else{
	print "$query\t$ab->{$query}\tbest\n";
    }
}


sub read_file{
    my ($file) = @_;
    my $best_target;
    my $best_bitscore;
    open( IN , $file );
    while(<IN>){
	my @A = split;
	my ($query, $target , $identity, $len, $evalue , $bit_score) =  ( $A[0] , $A[1] , $A[2] , $A[3] , $A[10], $A[11]) ;
	next if(  $evalue > 0.001 ); # ||  $identity < 99 ||  $len < 150 );
	if( exists $best_target->{ $query} ){
	    if( $best_bitscore->{$query} < $bit_score ){
		$best_target->{$query} = $target;
		$best_bitscore->{$query} = $bit_score;
	    }	    
	}
	else{
	    $best_target->{$query} = $target;
	    $best_bitscore->{$query} = $bit_score;
	}
	
	}
    close( IN );
    return $best_target;
}

 

sub read_anno{
	my ($file) = @_;
	my $hash;
	open( IN , $file );
	while( <IN> ){
		my @A =  split;
		my $key = shift(@A);
		$hash->{$key} = "@A";
	}
	close( IN);
	return $hash;
}
