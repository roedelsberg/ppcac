#!/usr/bin/perl

use strict;

my $usage = "

program transcript_file [60]

the longest transcript above a min protein length

";
die $usage unless( @ARGV > 0 );


my $file = shift @ARGV;

my $protein = 0;
my $min_length = 120;   # nucleotides
my $fold_length = 1;


my $prefix = $file;
$prefix =~ s/_.+//g;
print STDERR "$prefix\n";
my %genome = readFasta( $file );

my $seq_id = 1;
foreach my $g ( keys %genome ){
    my ($new_id, $protein) = denovo_orf( $genome{ $g } , $g  , $min_length, $fold_length);	
   
    if( $new_id ne ""){
	$g =~ s/_Confidence.+//g;
	#print ">".$prefix."$seq_id\n";
	#print ">".$prefix."-$new_id:$g\n";
	#print toLines( $protein );
	$seq_id++;
	
    } 
}



sub denovo_orf{
	my ( $seq , $refid, $min_length , $fold_length) = @_;
	my %L;
        my $id = "";
	foreach my $ori ( "f", "r"){
	    for( my $i = 0; $i < 3 ; $i++ ){

                my $cut = reading_frame( $seq, $i, $ori );  ## gets the sequence in right orientation and with the starting nucleotide
                my @starts = potential_starts( $cut) ;	
	
                my $N = length( $cut );
                        
                        foreach my $s( @starts ){
                                my $n = 0;
                                my $seq = "";
				my $e = $s;
                                for( $e; $e < $N; $e +=3 ){   ## translate until a stop codon is aquired
                                        my $codon = substr( $cut, $e, 3 );
                                        $seq = substr( $cut, $s, $e+3-$s );
                                        
                                        if( $codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA"){                     
                                                last;
                                                }
                                        }
                                $n = length($seq);


                                $id = "$refid:".($s+$i).":".($e+$i+1+3).":+";  ## i is the offset from reading frame
				$id = "$refid:".($s+$i).":".($e+$i+1+3).":-" if( $ori eq "r"); #	print STDERR "$id\n" if( $s+$i==2);
				$L{$id} = $seq ;           
				#print STDERR "$i\t$s\t$id\t".length($seq)."\n";
                        }
                
           }     
        }
	#print STDERR scalar (keys %L)."\n";

        my @sorted_by_length = sort{ (-1) *( length($L{$a}) <=> length($L{$b})) } keys %L;
       
	my $first      = $sorted_by_length[0];
	my $second    =    $sorted_by_length[1];
	my $n1 = length( $L{$first} );
	my $n2 = length( $L{$second} );
	
	if( $n1 / $n2 >= $fold_length &&  $n1 >= $min_length ){ 
	print ">".$refid."\n";
	my $protein =  translate( $L{$first} );
	print  toLines($protein);
	}


}



sub potential_starts
## gets all potential starting positions, includin 0 for a truncated sequence
{
	my ($seq) = @_;
	my @starts = ( 0 ) ;

	for( my $k = 0; $k < length($seq) ; $k+=3 ){   ### check all potential starts
	    push( @starts, $k) if( substr( $seq, $k, 3 ) eq "ATG");	
	}			

	
	return @starts;
}



sub reading_frame
## gets the sequence in right orientation and with the starting nucleotide
{
	
my ( $seq, $i, $st ) = @_;
my $cut;
if( $st eq "f"){
		$cut = substr( $seq , $i );	
		}
else{
	$cut = substr( revert( $seq ), $i );	
	}

return $cut;
}


sub readFasta{
	
	my ( $database ) = @_;
	
	my %genome;
	open(DATA,"<$database") || die "Unable to open $database!\n";
	
	my $sequence = "";
	my $firstLine = <DATA>;    
	my $name = getFastaName( $firstLine ); 
	while (<DATA>){
		if( $_ =~ m/\>/ ){                      #  IF THE LINE STARTS WITH ">", SAVE THE SEQUENCE
		$genome{ $name } = $sequence; 	      
		$name = getFastaName($_); $sequence = "";    #  INITIALIZE NAME AND SEQUENCE
		}
		else{
			my $line = uc $_;  chomp $line;    # CONVERT TO UPPER CASE
			$sequence .= $line;
		}
	}
	$genome{ $name } = $sequence;    # OTHERWISE WE OMIT THE LAST SEQUENCE
	close (DATA); 

	return %genome;
}




sub getFastaName{
 my ($sequence) = @_;
 my @data = split ( /\s/ , $sequence );
 @data = split ( /\>/ , $data[0] );
 print "no FASTA format found\n" if ( !$sequence =~ /\n/ );
 return $data[1];
}




sub revert{
    my ($dna) = @_ ;
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    return reverse($dna);
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

sub translate{
    my ($seq) = @_;
    my $cdr3 = "";
    my %translater = genetic_code();
    my $okay =1;
    for(my $i = 0; $i <= length($seq)-3; $i = $i + 3 ){
        my $codon = substr($seq, $i, 3 );
        my $aa = $translater{ $codon };
        if( $aa ne "" ){
        	$cdr3 .=  $aa;
        }
    }
    return ( $cdr3);

}


sub genetic_code{
    
    my %code = (
        "TTT" => "F",
        "TTC" => "F",
        "TTA" => "L",
        "TTG" => "L",
        "TCC" => "S",
        "TCT" => "S",
        "TCA" => "S",
        "TCG" => "S",
        "TAT" => "Y",
        "TAC" => "Y",
        "TAA" => "*",
        "TAG" => "*",
        "TGA" => "*",
        "TGT" => "C",
        "TGC" => "C",
        "TGG" => "W",
        "CTT" => "L",
        "CTC" => "L",
        "CTA" => "L",
        "CTG" => "L",
        "CCA" => "P",
        "CCC" => "P",
        "CCG" => "P",
        "CCT" => "P",
        "CAT" => "H",
        "CAC" => "H",
        "CAA" => "Q",
        "CAG" => "Q",
        "CGT" => "R",
        "CGG" => "R",
        "CGC" => "R",
        "CGA" => "R",
        "ATT" => "I",
        "ATC" => "I",
        "ATA" => "I",
        "ATG" => "M",
        "ACT" => "T",
        "ACG" => "T",
        "ACC" => "T",
        "ACA" => "T",
        "AAT" => "N",
        "AAC" => "N",
        "AAA" => "K",
        "AAG" => "K",
        "AGT" => "S",
        "AGC" => "S",
        "AGA" => "R",
        "AGG" => "R",
        "GTT" => "V",
        "GTC" => "V",
        "GTA" => "V",
        "GTG" => "V",
        "GCT" => "A",
        "GCC" => "A",
        "GCA" => "A",
        "GCG" => "A",
        "GAT" => "D",
        "GAC" => "D",
        "GAA" => "E",
        "GAG" => "E",
        "GGT" => "G",
        "GGC" => "G",
        "GGA" => "G",
        "GGG" => "G"
        );

    return %code;
}

        
