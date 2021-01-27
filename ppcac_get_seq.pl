#!/usr/bin/perl

use strict;

my $usage = "program multi.fa n [-l predefined_list of fasta seqs] [-re (choose random elements from a non fasta list)]";
my $description = "Returns n random sequences from a fasta file or selected IDs from a list";
if( scalar(@ARGV) == 0 ){ print "$description\n\n$usage\n"; exit 0;}

my $n  = $ARGV[1];

my $max = 0;

if( $ARGV[2] eq "-re" ){
    $max = `wc -l $ARGV[0] `; chomp $max;
}
else{
$max = `grep ">" $ARGV[0] | wc -l`; chomp $max;
}

print STDERR "$max\n";

my %indices;

for( my $i = 0; $i < $n ; $i++){
	my $x = int( rand( $max ) );
	while( exists $indices{ $x } ){
		$x = int( rand( $max ) );
	}
	
	$indices{ $x } = 1;
}

if( $ARGV[2] eq "-l"){
	parseFastaFromList( $ARGV[0], $ARGV[3]);	

}
else{
	if($ARGV[2] eq "-re"){  # choose random element from a list file
			$max = `wc -l $ARGV[0]`; chomp $max;
			%indices = {};
			my @lines = `less $ARGV[0]`;
			my %already = {};
			for( my $i = 0; $i < $n ; $i++){
				my $x = int( rand( $max ) );
				while( exists $already{$x} ){
					$x = int( rand( $max ) );
					}
				print $lines[$x];
				$already{ $x } = 1;	
				}	  
				
		}
	else{	
		parseFasta( $ARGV[0], \%indices );	
	}	
}
######################################################
###################  PARSE FASTA     #################
######################################################

sub parseFasta{
	my ( $file, $indices  ) = @_;

my $sequence = "";
my $header = "";
my $positions;
my $first = 1;
my $i = 0;	
open( FASTA, "$file" );

while( my $line = <FASTA> ){
	#chomp $line;	
	
	if( ! $first ){	  				# $sequence and header are not initialized in first iteration 
		
		if(  $line =~ /\>/   ){  	# encounters a new sequence 		
			if( exists $indices->{$i} ){
				print ">".$header;
				print toLines($sequence);
			}
		
			$i++;
			$header = $line; $header =~ tr/\>//d;  # read a new header line
			$sequence = "";
			}
			
		else{ 
			$line =~ s/\s+//g;
			$sequence .= $line;   	# concatenate sequences
			}
		}
	else{                    		# read the header line
		if(  $line =~ /\>/   ){ $header = $line; $header =~ tr/\>//d;
			}
		}	
	$first = 0;
	}
	close( FASTA );
	if( exists $indices->{$i} ){
				print ">".$header;
				print toLines($sequence)."\n";	
			}
	
}



sub toLines{
        
        my ($dna) = @_;
        my $l = length( $dna );
        my $i = 0;
        my $ret = "";
        my $d = 60; 
        while( $i < $l) {
                my $line =  substr($dna, $i, 60). "\n";
                $ret .= $line;
                $i = $i + $d;
        } 
        return $ret;
        
        
}


sub parseFastaFromList{
	my ( $file, $list ) = @_;

my %queries = readHash( $list );
my $sequence = "";
my $header = "";
my $rest = "";
my $positions;
my $first = 1;
my $i = 0;	
my %already;
open( FASTA, "$file" );

while( my $line = <FASTA> ){
	#chomp $line;	
		
	if( ! $first ){	  				# $sequence and header are not initialized in first iteration 
		
		if(  $line =~ /\>/   ){  	# encounters a new sequence 		
			if( exists $queries{$header} && ! exists $already{$header} ){
				print ">".$header."\n";
				print toLines($sequence);	
				$already{$header} = 1;
			}
		
			$i++;
			( $header, $rest ) = split( /\s+/, $line );
			$header =~ tr/\>//d;  # read a new header line
			 #print STDERR "$header\n";
			$sequence = "";
			}
			
		else{ 
	 $line =~ s/\s+//g;
		$sequence .= $line;   	# concatenate sequences
			}
		}
	else{                    		# read the header line
		if(  $line =~ /\>/   ){ 
				( $header, $rest ) = split( /\s+/, $line );
				$header =~ tr/\>//d;
			#	print STDERR "$header\n";	
			}
		}	
	$first = 0;
	}
	close( FASTA );
	
	if( exists $queries{$header} && ! exists $already{$header} ){
				print ">".$header."\n";
				print toLines($sequence)."\n";	
			}
	else{
		; 
         }
}	
	
	sub readHash{
        my ( $file ) = @_ ;
        my %hash;
        open( IN, $file );
        while( my $line = <IN>){
        if( length $line > 3){
        	my ($id, $loc ) = split( /\s+/, $line );
        	$hash{ $id } = 1;
#        	print STDERR  "$id\n";	
        	}
        }
        close(IN);
        return %hash;
}
	


