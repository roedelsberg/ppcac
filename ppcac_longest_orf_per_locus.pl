#!/usr/bin/perl

use strict;

my $table;


die "program features.loc ORF_length.stats exon_number.txt\n" if( @ARGV < 2);

my $orf_length;
my $table;
my $forbidden;
my $good_ids;  ### list of transcripts that has been chosen 
my $exon_number;
my $min_exon_number = 3;
my $min_orf_length = 60;
my $query = "XXXYXXXYYXYXYYX";


#####################################
#
#  Read file with exon numbers
#
####################################

open( IN3, $ARGV[2] );
while(<IN3>){
    my ($id, $n_exons) = split;
    $exon_number->{ $id } = $n_exons;
}
close(IN3);


#####################################
#
#  Read ORF file with peptide lengths
#
####################################

open( IN2, $ARGV[1] );
while(<IN2>){
    my ($iid, $length ) = split;
    #my $iid = $id."_". $isoform;
    $orf_length->{ $iid } = $length;
    $forbidden->{$iid} = 1 if( $length < $min_orf_length );
}
close(IN2);



#####################################
#
#  Read exon coordinates
#  single exon -> all transcripts sharing this exon
#
####################################

open( IN, $ARGV[0] );
while(<IN>){
    my ($id, $locus, $type ) = split;
    ### require only one exon boundary to be identical
    my ($contig, $s, $e, $o ) = split(/:/, $locus );
    my $win = 100;
    $s = int($s/$win)*$win;
     $e = int($e/$win)*$win;
    my $loc1 = "$contig:$s:$s:$o";
    my $loc2 = "$contig:$e:$e:$o";
    $id =~ s/:\d+\b//g;  ## truncate the identifier to match the gene ID
    if( $type eq "exon" && ! exists $forbidden->{$id} ){	
	#$table = update( $table, $locus, $id, 1 );
	$table = update( $table, $loc1, $id, 1 );
	$table = update( $table, $loc2, $id, 1 );
    }
    
}
close(IN);

#####################################
#
# First handle all cases of shared exons
#
#####################################
##### the order of loci  has to be sorted by coverage of the exon
my $exon_cov;
foreach my $locus( keys %$table ){
    my $hash = $table->{$locus};
    $exon_cov->{ $locus } = scalar( keys %$hash );
}
my @loci = sort{ (-1) * ($exon_cov->{$a} <=> $exon_cov->{$b}); }keys %$exon_cov;


foreach my $locus(@loci ) { 
    my $hash = $table->{$locus};
    
    if( scalar( keys %$hash ) > 1 ){   ## look at exons that are shared by multiple transcripts
	my @ids = sort{ (-1) * ($orf_length->{ $a } <=> $orf_length->{ $b }); } keys %$hash;  ## sort by ORF lengths
	my $forbid = 0;
	foreach my $id ( @ids ){
	    if( $good_ids->{$id} ){   ## ignore all 
		$forbid = 1;
	    }
	}
	
	foreach my $id ( @ids ){
	  
	    
	    if( $forbid ){
		$forbidden->{$id} = 1;
	    }
	    else{
		if( $exon_number->{$id} < $min_exon_number ){
		    $forbidden->{$id} = 1;
		} 
		else{    
		    unless(exists $forbidden->{$id} ){ ## once the longest ORF has been chosen disable all other ORFs/transcripts sharing this exon
			print "$id\n";
			$good_ids->{$id}  = 1; 
			$forbid = 1;
			$forbidden->{$id} = 1;
		    }
		}
	    }	    	    
	}		
    }
}


foreach my $locus( keys %$table ){
     my $hash = $table->{$locus};
    if( scalar( keys %$hash ) == 1 ){  ## now deal with exons that are in only a single transcript
	my $forbid = 0;
	
	foreach my $id ( keys %$hash ){
	
	   
	    if( $exon_number->{$id} < $min_exon_number ){
		    $forbidden->{$id} = 1;
	    }
	    else{
		unless(exists $forbidden->{$id} ){
		    unless( exists $forbidden->{$id}  ){
			print STDERR "$id\n"  if( $id=~ /$query/ );
			print "$id\n" ;
			$good_ids->{$id}  = 1; 
			$forbidden->{$id} = 1;
		    }
		}
	    }
	}
    }

}





sub update{
        my ($hash, $q, $t , $v ) = @_;
        if( exists $hash->{$q} ){
                $hash->{$q}->{$t} = $v ;
        } 
        else{
                my %new_hash = ( $t => $v );
                 $hash->{$q} = \%new_hash;
        }
        return $hash;
}

