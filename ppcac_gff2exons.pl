#!/usr/bin/perl

use strict;


=head

	Script for extracting exons and exon boundaries from the Pristionchus gff file
	plotting exon number

=cut

die "program gff3\n" if( @ARGV < 1 );
# $file = "/Users/roedel/Data/Annotations/pristionchus_proteomics_processed_for_ccodon.gff3";

my $file = $ARGV[0];

my $exons_per_gene;
my $strandof ;
my $mrna_start;
my $mrna_end;
my $contig_of;
my $segments = $file;
$segments =~ s/\.\S+//g;
$segments .= ".loc";

open(IN, $file);
while(<IN>){
	my ($chr,$source,$class, $start, $end, $score, $strand, $offset, $comment ) = split( /\t/ );
	
	if( $class eq "exon" ){  ## CDS (C.elegans) or exon (most others)
		my $hash = parseComment( $comment );
		my $mrna_id = $hash->{"Parent"};
		#my $mrna_id = $hash->{"ID"};  # C_briggsae
		($mrna_start, $mrna_start, $contig_of, $strandof , $exons_per_gene) = update_mrna( $mrna_start, $mrna_start, $exons_per_gene, $mrna_id, $start, $end , $contig_of, $strandof, $chr,  $strand);	
	}
	
}
close( IN );

open( OUT, ">$segments");
foreach my $mrna ( keys %$exons_per_gene ){    ## intron
	my @elist = keys %{$exons_per_gene->{ $mrna }};
	my @sorted = sort{ compareLocs($a,$b) } @elist;
	
	my $coding_length = 0;	
	my $chr = $contig_of->{ $mrna };	
	my $s   = $mrna_start->{ $mrna };
	my $e   = $mrna_end->{ $mrna };
	print OUT "$mrna\t$chr:$s:$e:".$strandof->{$mrna}."\tgenomic\n";
	my ($chr,$s,$e) = split(/:/, $sorted[0]); 
	my $last_end = $e+1;
	print OUT $mrna.":$s\t$chr:$s:$e:$strandof->{$mrna}\texon\n";
	$coding_length += $e-$s+1;
	for( my $i = 1; $i < @sorted; $i++ ){
		my $exon = $sorted[$i];
		my ($chr,$s,$e) = split( /:/, $exon );
		my $intron_end = $s-1;
		print OUT "$mrna:$last_end\t$chr:$last_end:$intron_end:".$strandof->{$mrna}."\tintron\n";
		print OUT $mrna.":$s\t$chr:$s:$e:$strandof->{$mrna}\texon\n";
		$last_end = $e+1;
		$coding_length += $e-$s+1;
	}	
	print "$mrna\t".@sorted."\t$coding_length\t".(abs($mrna_end->{$mrna}-$mrna_start->{$mrna})+1)."\n";
	
}
close( OUT );
print STDERR "Locations written in $segments\n";

 sub update_mrna{
		my  ($mrna_start, $mrna_start, $exons_per_gene, $mrna_id, $start, $end , $contig_of, $strandof, $chr,  $strand) = @_;
		$mrna_start->{$mrna_id} = $start if( !(exists $mrna_start->{$mrna_id} ) || $start < $mrna_start->{$mrna_id} );
		$mrna_end->{$mrna_id}  = $end    if( !(exists $mrna_end->{$mrna_id} )   ||  $end > $mrna_end->{$mrna_id} ); 
		$exons_per_gene = update( $exons_per_gene,  $mrna_id , "$chr:$start:$end", 1 );	
		unless( exists $contig_of->{$mrna_id} ){
				$contig_of->{$mrna_id} = $chr;
				$strandof->{ $mrna_id} = $strand;
		}	
		return ($mrna_start, $mrna_start, $contig_of, $strandof, $exons_per_gene);
}



sub update{
        my ($hash, $q, $t , $value) = @_;
        if( exists $hash->{$q} ){
                $hash->{$q}->{$t} = $value ;
        } 
        else{
                my %new_hash = ( $t => $value );
                 $hash->{$q} = \%new_hash;
        }
        return $hash;
}




sub parseComment{
	my ($text) = @_;
	my @A = split( /;/, $text );
	my $hash;
	foreach my $a( @A ){
		my ($key,$value) = split( /=/, $a );
		#print STDERR "$key\t$value\n";
		chomp $value;
		$hash->{$key} = $value;
		
	}
	
	return $hash;
}


sub compareLocs{
	my ($a,$b) = @_;
	my ($seqA, $startA, $stopA, $oriA ) = split( /:/, $a );
	my ($seqB, $startB, $stopB, $oriB ) = split( /:/, $b );
	
	if( $seqA ne $seqB ){
		return $seqA cmp $seqB;
		}
	else{
		return $startA <=> $startB;		
		}
	
}

