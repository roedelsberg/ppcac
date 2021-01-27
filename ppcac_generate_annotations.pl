#!/usr/bin/perl

use strict;

my $sequences;
my $mrna2cds;
my $minimal_exon_number = 0;
my $min_cds_length = 0;

die "program gtf genome \n


\n" if( @ARGV != 2);



my $genome_fasta =  $ARGV[1] ;
my %genome = read_genome( $genome_fasta );


my $annotation = read_gtf( $ARGV[0] , \%genome );
foreach my $mrna( keys %$annotation ){
	my $cds =  $annotation->{ $mrna };
	my $N = 0;
	foreach my $exon( keys %$cds ){

		my ($chr, $sta, $end, $ori) = split( /:/, $exon );	
		$N+= $end - $sta;
	}

}



my $correct = 0;
my $false = 0;


open( OUT2, ">RefSeq_mrnas.fa");
open( OUT, ">proteins.fa");
print STDERR (keys %$annotation)."\n";
foreach my $mrna( keys %$annotation){

	        my $seq = uc($sequences->{$mrna});
                
	        my $exons = $annotation->{ $mrna };
		my @Ex = sort{ compareLocs( $a , $b);} keys %$exons; # A := exons
		my ($chr, $start, $end, $ori)  = split(/:/, $Ex[0]);
		my $dna = $genome{$chr};		
	        my $reconstructed = get_mrna( $exons, $dna );
	#print STDERR "$mrna \t".@Ex."\n".toLines($reconstructed)."\n";	
	if( @Ex > $minimal_exon_number ){   
	       	my ($newid, $prot) = new_orf( $reconstructed, $mrna,  $min_cds_length );
		if( $newid ne ""){
		    if( $ori eq "-" || $ori eq "+"){
			my ($buffer_gff, $buffer_seq) = parse_annoations( $mrna, $newid , $dna, $reconstructed, $exons);
			if( $buffer_gff ne "" && ($buffer_gff =~ /exon/ )){
			    $correct++;
			    print $buffer_gff;
			    print OUT2 ">$mrna\n".toLines(  $reconstructed)."\n";
			    print OUT ">$mrna\n".toLines(  $prot )."\n";
			}
			else{
			    print STDERR "$mrna\t$newid\t\nExpected:\n".toLines( $prot )."\n";
			    print STDERR "Observed:\n".toLines( $buffer_seq )."\n";
			    print STDERR "mRNA:\n".toLines( $reconstructed )."\n";
			    $false++;
			    #exit 0;
			}
		    }
		}

	    }
	else{   
	    # single exon gene
            # get CDS 2 * longer than second best
	    # if positive, use the same pipeline
	    # if negative, revert mrna, exonerate, annotate

	    my ($newid, $prot) = denovo_orf( $reconstructed, $mrna );
	    if( $newid ne ""){
		my $ori = get_strand( $newid );
		if( $ori eq "+"){
		    my ($buffer_gff, $buffer_seq) = parse_annoations( $mrna, $newid , $dna, $reconstructed, $exons);
		    if( $buffer_gff ne "" && ($buffer_gff =~ /exon/ )){		      
			    $correct++;
			    print $buffer_gff;
			    print OUT2 ">$mrna\n".toLines(  $reconstructed)."\n";
			    print OUT ">$mrna\n".toLines(  $prot )."\n";
		    }
		    else{
			 $false++;
		    }
		}
		else{        ## 
		    my $tmp = revert( $reconstructed);
	            my $exons = reannotate( $mrna, $tmp,  $genome_fasta  );

		    if( scalar keys %$exons > 0 ){
			 my $reconstructed = get_mrna( $exons, $dna );
			 my ($newid, $prot) = new_orf( $reconstructed, $mrna,  $min_cds_length );
			 if( $newid ne ""){
			     my ($buffer_gff, $buffer_seq) = parse_annoations( $mrna, $newid , $dna, $reconstructed, $exons);
			     if( $buffer_gff ne "" && ($buffer_gff =~ /exon/ ) ){ 
				 print $buffer_gff;
				 print OUT2 ">$mrna\n".toLines(  $reconstructed)."\n";
				 print OUT ">$mrna\n".toLines(  $prot )."\n";
				 $correct++;
			     }
			     else{
				 $false++;
			     }
			 }

		     }		  
		}
	    }
	}

}
close( OUT );
close( OUT2 );

print STDERR "correctly_read $correct $false\n";




###############################################
###############################################
##
##         subroutines 
##
###############################################
###############################################

sub get_strand{
    my ($id ) = @_;
    my ($mrna,$s,$e,$o) = split(/:/, $id );
    return $o;
}

sub get_mrna{
    my ( $exons, $dna ) = @_;
    my @A = sort{ compareLocs( $a , $b);} keys %$exons;
    my $buffer ="";
    foreach my $exon ( @A ){
	    my ($contig, $s, $e, $o) = split(/:/, $exon);
	    if( $o eq "+"){
		$buffer .=  substr( $dna, $s-1, $e-$s+1);
	    }
	    else{
		$buffer .=  revert(substr( $dna, $s-1, $e-$s+1));
	    }
    }
    return $buffer;
}


sub parse_annoations{
	my ( $mrna, $protein_id, $dna, $seq, $exons ) = @_; 

         my @A = sort{ compareLocs( $a , $b);} keys %$exons;
         my ($tmp,$cds_start,$cds_end,$other) = split(/:/, $protein_id);
	#print STDERR (($cds_end-$cds_start+1)/3)."\n";
	 #my $cds_length = $cds_end - $cds_start+1;
         my ($chr, $start, $end, $ori)  = split(/:/, $A[0]);
         my $buffer="";
         my $phase=0;
		if( $ori eq "+" ){
		
			my $tlength = 0;
			my $state = "5utr";
	  		my $buffer_seq = "";
                        my $buffer_gff = "";
 			my ($first, $last) = (0,0 );  ## first and last position of the gene

			foreach my $exon ( @A ){   ## search cds start
                        	my ($contig, $s, $e, $o) = split(/:/, $exon);
				$first = $s if( $first == 0 || $s < $first);  ## update the first position of the gene
                                $last = $e   if( $last == 0 || $e > $last);   ## update the last position of the gene
				last if( $s == $e);      ### this is weird but, but it occured in CUFF_3761_4

				if( $state eq "5utr"){
					if( $tlength + $e-$s+1  >= $cds_start ){
						my $xstart = $s + ($cds_start - $tlength);  ## -1 
						 my $rest = $xstart - $s;
						$buffer_gff .= "$chr\tCoding_transcript\tfive_prime_UTR\t$s\t".($xstart-1)."\t.\t$ori\t.\tParent=$mrna\n"  if( $rest > 0 );
						$state = "CDS";
						if(  $tlength + $e-$s+1  >=  $cds_end ){ ## single exon gene -> go into CDS mode
						   
						    $tlength += $rest;
						    $s = $xstart;
						}
						else{	
						   
						    my $exon_seq = substr( $dna, $xstart-1, $e-$xstart+1);

						    if(  $tlength + $e-$s+1  <  $cds_end ){ 
						
							$buffer_gff .= "$chr\tCoding_transcript\texon\t$xstart\t$e\t.\t$ori\t0\tParent=$mrna\n" unless($xstart > $e);  ## first exon has a phase of zero , unless statement add 2019-05
							$buffer_seq .=  translate( $exon_seq);
						    }
						    $tlength += $e-$s+1;
						    $phase = set_phase( $exon_seq  );
						    
					 	    next; #unless(  $tlength + $e-$s+1  >= $cds_end );
						}
					} 
					else{
					    
						$buffer_gff .= "$chr\tCoding_transcript\tfive_prime_UTR\t$s\t$e\t.\t$ori\t.\tParent=$mrna\n";
					 	$tlength += $e-$s+1;
						}	
                			}
				if( $state eq "CDS" ){
					if(  $tlength + $e-$s+1  >=  $cds_end ){
						my $xend = $s + ($cds_end - $tlength);
						 my $exon_seq = substr( $dna, $s+$phase-1, $xend-$s+1-$phase);
                                                $buffer_gff .=  "$chr\tCoding_transcript\texon\t$s\t$xend\t.\t$ori\t$phase\tParent=$mrna\n";
						my $rest = $e-$xend;
						$buffer_gff .=  "$chr\tCoding_transcript\tthree_prime_UTR\t".($xend+1)."\t$e\t.\t$ori\t.\tParent=$mrna\n" if( $rest > 0 ); 
						$buffer_seq .= translate(  $exon_seq );
						$state="3utr";	
						next;
					}
					else{
						my $exon_seq = substr( $dna, $s +$phase-1, $e-$s+1-$phase);
                                                $tlength += $e-$s+1;
						$buffer_gff .=  "$chr\tCoding_transcript\texon\t$s\t$e\t.\t$ori\t$phase\tParent=$mrna\n";
						$buffer_seq .=   translate($exon_seq);
						$phase =  set_phase( $exon_seq );
						last if( length($seq) < $tlength+3);  ###  CUFF_3761_4
						}
					}
				if($state eq "3utr" ){
					$buffer_gff .=  "$chr\tCoding_transcript\tthree_prime_UTR\t$s\t$e\t.\t$ori\t.\tParent=$mrna\n";
		     
					}
				}
					
	 		my $n_stop = stops( $buffer_seq );
                        if( $n_stop <= 0 ){
			                $buffer_gff = "$chr\tCoding_transcript\tmRNA\t$first\t$last\t.\t$ori\t.\tID=$mrna\n". $buffer_gff;
					return ( $buffer_gff , $buffer_seq);
                                 }
                        else{
			    print STDERR "$mrna\n".length($seq)."\n$buffer_gff\n";
			    return("", substr($buffer_seq, 0, 600) );
                                }
			}
		else{                                    ### minus strand
			my $tlength = 0;
                        my $state = "5utr";
			my $buffer_seq = "";
			my $buffer_gff = "";
			my ($first, $last) = (0,0 );
                        foreach my $exon ( @A ){   ## search cds start
                                my ($contig, $s, $e, $o) = split(/:/, $exon);
			#	$s++;                                              ## special for cufflinks minus strand
			#	$e++;
				$first = $s if( $first == 0 || $s < $first);
				$last = $e if( $last == 0 || $e > $last);
                                if( $state eq "5utr"){
                                        if( $tlength + $e-$s+1  >= $cds_start ){
                                                my $xend = $e - ($cds_start - $tlength);
                                                my $rest = $e - $xend;
                                                $buffer_gff .= "$chr\tCoding_transcript\tfive_prime_UTR\t".($xend+1)."\t$e\t.\t$ori\t.\tParent=$mrna\n"  if( $rest > 0 );
                                               
                                                $state = "CDS";
                                                if(  $tlength + $e-$s+1  >= $cds_end ){   ## 5utr and 3 utr may be in same exon
							$e = $xend;
							$tlength += $rest; 
						        ##$buffer_gff .= "5'UTR TO END $tlength\n";
						} 
						else{	### first coding exon will always has a phase of 0
						        my $exon_seq =  revert(substr( $dna, $s-1, $xend-$s+1));
							$tlength += $e-$s+1;
							 if( $xend != $s-1  ){   ## remove artificial exons of size -1
							     $buffer_seq .=  translate( $exon_seq ); 
							     $buffer_gff .= "$chr\tCoding_transcript\texon\t$s\t".($xend)."\t.\t$ori\t0\tParent=$mrna\n";
							 }
							$phase =  set_phase( $exon_seq );
							next;
							}
                                        }
                                        else{
					    $tlength += $e-$s+1;
                                                $buffer_gff .= "$chr\tCoding_transcript\tfive_prime_UTR\t$s\t$e\t.\t$ori\t.\tParent=$mrna\n";
                                                
                                                }
                                        }
				  if( $state eq "CDS" ){
						 if( $tlength + $e-$s+1  >= $cds_end ){	
							my $xstart = $e - ($cds_end - $tlength);
							my $exon_seq = revert(substr( $dna, $xstart-1, $e-$xstart+1)) ;
							$exon_seq = substr( $exon_seq, $phase, length($exon_seq) - $phase );
							$buffer_gff .= "$chr\tCoding_transcript\texon\t$xstart\t$e\t.\t$ori\t$phase\tParent=$mrna\n";
                                                	my $rest = $xstart-$s;
                                                	$buffer_gff .= "$chr\tCoding_transcript\tthree_prime_UTR\t$s\t".($xstart-1)."\t.\t$ori\t.\tParent=$mrna\n" if( $rest > 0 ); 
                                                	 $buffer_seq .= translate( $exon_seq);
							$state="3utr";
                                                	next;

							
							}
						else{

							my $exon_seq = revert(substr( $dna, $s-1, $e-$s+1)) ;
							$exon_seq = substr( $exon_seq, $phase, length($exon_seq) - $phase );
							$tlength += $e-$s+1;
							$buffer_gff .= "$chr\tCoding_transcript\texon\t$s\t$e\t.\t$ori\t$phase\tParent=$mrna\n";
							$buffer_seq .= translate($exon_seq);
                                                	$phase = set_phase( $exon_seq );
                                                	

						}
					}	
				    if( $state eq "3utr"){
						 $buffer_gff .= "$chr\tCoding_transcript\tthree_prime_UTR\t$s\t$e\t.\t$ori\t.\tParent=$mrna\n";
					}		
				}##  end of foreach exon
						
			my $n_stop = stops( $buffer_seq );
			if( $n_stop <= 0 ){
			        $buffer_gff = "$chr\tCoding_transcript\tmRNA\t$first\t$last\t.\t$ori\t.\tID=$mrna\n". $buffer_gff;
				return ( $buffer_gff , $buffer_seq);
				}
			else{
			         print STDERR "$buffer_gff\n";
			        return ("", substr($buffer_seq, 0, 300)) ;
				}	
			}



}	

sub stops{
	my ($string) = @_;
	my @A = split( // , $string );
	my $n = 0;
	for( my $i = 0; $i < @A-1 ;$i++ ){
		$a = $A[$i];
		$n++ if( $a eq "*");
		}
	return $n;	

}

sub compareLocs{
	my ($a,$b) = @_;
	my ($seqA, $startA, $stopA, $oriA ) = split( /:/, $a );
        my ($seqB, $startB, $stopB, $oriB ) = split( /:/, $b );
#	print STDERR "$a vs $b\n";        

        if( $seqA ne $seqB ){
                return $seqA cmp $seqB;
                }
        else{
		if( $oriA eq "+"){
	                return $startA <=> $startB;             
			}
		else{
			return (-1) * ( $startA <=> $startB );
			}
                }
        


}

sub read_gtf{
	my ($file , $genome) = @_;
	print STDERR "Reading $file\n";
	my $genes;
	my $alignment_count;
	open( IN , $file );
	while( <IN> ){
		my ($chr, $source, $type, $start, $end, $scor, $ori, $phase,$descr  ) = split( /\t/, $_);
		#print STDERR $_;
		next unless( exists $genome->{$chr} );
	        if( $type eq "gene"){
		    my @B = split(/;\s*/, $descr );
		    my $attributes;
		    foreach my $b( @B ){
			my ($key, $value) = split( /=/, $b);
			$value =~ s/\"//g;
			$attributes->{$key} = $value;
		    }
		    my $mrna = $attributes->{ "ID" };
		    $mrna =~ s/\s+//g;
		    $alignment_count->{ $mrna }++;
		    #print STDERR "$type\t$mrna\n";
                }
	        if( $type eq "CDS" ){
				my @B = split(/;\s*/, $descr );
				my $attributes;
				foreach my $b( @B ){
					my ($key, $value) = split( /=/, $b);
					$value =~ s/\"//g;
					$attributes->{$key} = $value;
				}
				my $mrna = $attributes->{ "Parent" };
#				$mrna =~ s/\./_/g;
				 $mrna =~ s/\s+//g;
			
				my $exon = "$chr:$start:$end:$ori";
				#print STDERR "$mrna $exon\n";
				$genes = update( $genes, $mrna, $exon, 1);
			}
		}
	close( IN);
	foreach my $key ( keys %$alignment_count ){
	    if( $alignment_count->{ $key} > 1  ){
		delete   $genes->{ $key};
	    }
	}
	return $genes;
}

sub read_single_gff{
	my ($file ) = @_;
	my $exons;
	open( IN , $file );
	while( <IN> ){
		my ($chr, $source, $type, $start, $end, $scor, $ori, $phase,$descr  ) = split( /\t/, $_);
		if( $type eq  "CDS" ){
				my @B = split(/;/, $descr );
				my $attributes;
				foreach my $b( @B ){
					my ($key, $value) = split( /=/, $b);
					$attributes->{$key} = $value;
				}
				my $mrna = $attributes->{ "Parent" };
				my $exon = "$chr:$start:$end:$ori";
				$exons->{$exon} = 1;
			}
		}
	close( IN);
	return  $exons;
}

sub set_phase{
    
    my ($exon_seq) =@_;
    my	$phase = length( $exon_seq )%3;
        $phase = 3 - $phase;
        $phase = 0 if( $phase == 3 );
    return $phase;
}

sub reannotate{
my ($mrna, $seq, $genome ) = @_;

my $file = "mrna.fa";
open( OUT,">$file");
print OUT ">$mrna\n".toLines($seq);
close( OUT );#				
# my $genome = "exspectatus_allpath1_gap_closed.fa";
my $call = "exonerate  -m est2genome --showtargetgff TRUE --maxintron 20000 --bestn 1 --dnawordlen 20  --showalignment FALSE --showvulgar false -q $file  -t $genome > exonerate.gff";
system("$call");
system("cu_exonerate2gff.pl  exonerate.gff  >  exonerate.gff3");
system($call);
return read_single_gff( "exonerate.gff3" );

}

sub new_orf{
	my ( $seq , $refid, $min_length ) = @_;
	my %L;
        my $id = "";
        for( my $i = 0; $i < 3 ; $i++ ){  ## only on one strand

                my $cut = reading_frame( $seq, $i, "f" ); ## cut the sequence according to the reading frame
                my @starts = potential_starts( $cut) ;
		push(  @starts , 0 );#  if( $i == 0);  ## modified on the 2013-08-20 , a 5' truncated gene has to be considered at any position
	#	print STDERR "@starts\n";
                my $N = length( $cut );
                        
                        foreach my $s( @starts ){
                                my $n = 0;
                                my $seq = "";
				my $e = $s;
                                for( $e; $e < $N; $e +=3 ){   ## get length of transcript
                                        my $codon = substr( $cut, $e, 3 );
                                        $seq = substr( $cut, $s, $e+3-$s );
                                        $n = length($seq);
                                        if( $codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA"){                     
                                                last;
                                                }
                                        }
                                
                                if($n > $min_length ){      ## update hits

                                        $id = "$refid:".($s+$i).":".($e+$i+2).":+";  ## i is the offset from reading frame  ## +2 in order to include the last codon with inclusive coordinates
				#	print STDERR "$id\n" if( $s+$i==2);
                                        if( !exists $L{$id} ){
                                                $L{$id} = $seq ;                                        
                                                }
                                        else{
                                                $L{$id} = $seq  if( exists $L{$id} < $n );
                                                }
                                                
                                        }
                        }
                
                
        }
        my @sorted_by_length = sort{ (-1) *( length($L{$a}) <=> length($L{$b})) } keys %L;
#	print STDERR length($L{$sorted_by_length[0]})."\t".  length($L{$sorted_by_length[1]})."\n";
        if( @sorted_by_length  > 0 &&  length($L{$sorted_by_length[0]})> 1*  length($L{$sorted_by_length[1]}) ){ 
	  my $id = $sorted_by_length[0];	
	  my $protein =  translate( $L{$id} );
	   return ($id, $protein);	
        }
	return ("","");
}

sub denovo_orf{
	my ( $seq , $refid) = @_;
	print STDERR "denovo_orf\n";
	my %L;
        my $id = "";
        my $min_length = 180;
	foreach my $ori ( "f", "r"){
	    for( my $i = 0; $i < 3 ; $i++ ){

                my $cut = reading_frame( $seq, $i, $ori );
                my @starts = potential_starts( $cut) ;
		push(  @starts , 0 ) if( $i == 0);
                my $N = length( $cut );
                        
                        foreach my $s( @starts ){
                                my $n = 0;
                                my $seq = "";
				my $e = $s;
                                for( $e; $e < $N; $e +=3 ){   ## get length of transcript
                                        my $codon = substr( $cut, $e, 3 );
                                        $seq = substr( $cut, $s, $e+3-$s );
                                        $n = length($seq);
                                        if( $codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA"){                     
                                                last;
                                                }
                                        }
                                
                               # if($n > $min_length ){      ## update hits

                                        $id = "$refid:".($s+$i).":".($e+$i+1+3).":+";  ## i is the offset from reading frame
					$id = "$refid:".($s+$i).":".($e+$i+1+3).":-" if( $ori eq "r"); #	print STDERR "$id\n" if( $s+$i==2);
                                        if( !exists $L{$id} ){
                                                $L{$id} = $seq ;                                        
                                                }
                                        else{
                                                $L{$id} = $seq  if( exists $L{$id} < $n );
                                                }
                                                
                                #        }
                        }
                
           }     
        }
        my @sorted_by_length = sort{ (-1) *( length($L{$a}) <=> length($L{$b})) } keys %L;
        if( @sorted_by_length  > 0 &&  length($L{$sorted_by_length[0]})> 2 *  length($L{$sorted_by_length[1]}) &&
			length($L{$sorted_by_length[0]}) >= $min_length 
	 ){ ## best > 2 * second longest
	  my $id = $sorted_by_length[0];	
	  my $protein =  translate( $L{$id} );
	   return ($id, $protein);	
        }
	return ("","");
}



sub store_seq{
	my ($sequences,$id,$current_seq) = @_;
	 $sequences->{$id} = $current_seq;
#         print STDERR ">$id\n".toLines( $current_seq)."\n";
        ($current_seq, $id) = ("","");
	return ($sequences,$id,$current_seq);
}


sub getFastaName{
 my ($sequence) = @_;
 my @data = split ( /\s/ , $sequence );
 @data = split ( /\>/ , $data[0] );
 print "no FASTA format found\n" if ( !$sequence =~ /\n/ );
 return $data[1];
}

sub update{
        my ($hash, $q, $t , $v) = @_;
        if( exists $hash->{$q} ){
                $hash->{$q}->{$t} =$v ;
        } 
        else{
                my %new_hash = ( $t => $v );
                 $hash->{$q} = \%new_hash;
        }
        return $hash;
}


sub toLines{
        
        my ($dna) = @_;
#	print STDERR "tL: $dna\n";
        my $l = length( $dna );
        my $i = 0;
        my $ret = "";
        my $d = 60; 
        while( $i < $l) {
                my $line =  substr($dna, $i, 60). "\n";
#		print STDERR "tL:".$line;
                $ret .= $line;
                $i = $i + $d;
        } 
        return $ret;
        
        
}


sub potential_starts{
#  
#  start with first position
#  start with each methionine
#  start after each stop codon
#
        my ($seq) = @_;
        my @starts ;
#        my $last = 0;
        for( my $k = 0; $k < length($seq) ; $k+=3 ){   ### check all potential starts
	    my $codon =  substr( $seq, $k, 3 );
	    push( @starts , $k ) if( $codon eq "ATG" );

	 #   if( $codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA"){   
	#	my $j = $k + 3;
	#	if( $j < length($seq) - 3){
	#	    my $next_codon = substr( $seq, $j, 3 );
	#	     if( $next_codon ne "TAG" && $next_codon ne "TAA" && $next_codon ne "TGA"){   
	#		 push( @starts , $j );
	#	     }
	#	}
	#    }
	}
        
        return @starts;
}

sub reading_frame{
        
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



sub revert{
    my ($dna) = @_ ;
    #print STDERR "Reverting $dna ".reverse($dna)."\n";
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    return reverse($dna);
}


sub translate{
    my ($seq) = @_;
    my $cdr3 = "";
    my %translater = genetic_code();
    my $okay =1;
    for(my $i = 0; $i <= length($seq)-3; $i = $i + 3 ){
        my $codon = substr($seq, $i, 3 );
        my $aa = "X";
	$aa = $translater{ $codon } if( exists $translater{ $codon } );
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

        

sub read_genome{
    

    my ($database) = @_;
    print STDERR "Reading $database\n";
        my %genome;
        open(DATA,"<$database") || die "Unable to open $database!\n";
        
        my $sequence = "";
        my $firstLine = <DATA>;    
        my $name = getFastaName( $firstLine ); 
        while (<DATA>){
                if( $_ =~ m/\>/ ){                      #  IF THE LINE STARTS WITH ">", SAVE THE SEQUENCE
                $genome{ $name } = $sequence;         
                # print STDERR $name." ".substr( $sequence, 0, 10)."\n";  # TEST
                $name = getFastaName($_); $sequence = "";    #  INITIALIZE NAME AND SEQUENCE
                }
                else{
                        my $line = uc $_;  chomp $line;    # CONVERT TO UPPER CASE
                        $sequence .= $line;
                }
        }
        $genome{ $name } = $sequence;    # OTHERWISE WE OMIT THE LAST SEQUENCE
        #print STDERR $name." ".substr( $sequence, 0, 10)."\n";  # TEST
        close (DATA); 
        return %genome;
}

################################
#                              #
#   get Fasta Name             #
#                              #
################################


sub getFastaName{
 my ($sequence) = @_;
 my @data = split ( /\s/ , $sequence );
 @data = split ( /\>/ , $data[0] );
 print "no FASTA format found\n" if ( !$sequence =~ /\n/ );
 return $data[1];
}


