#!/usr/bin/perl




use strict;


my $usage = "program fasta ";
my $description = "Counts the sequence length and GC content

seqID N_bp\tGC\tN_called



 ";
if( scalar(@ARGV) == 0 ){ print "$description\n\n$usage\n"; exit 0;}




my $database = $ARGV[0];

my $speciesPrefix = "";


#Fasta:print;

###################################################################################
#                                                                                 #
#  READ GENOME SEQUENCE FILES INTO %GENOME                                       #
#                                                                                 #
###################################################################################


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

my $number_of_seqs = keys( %genome );
my $sequence = "undef";
my $previous = "undef";

my @A = sort{ (-1) * ( length($genome{$a}) <=> length($genome{$b})) } keys %genome;

my $N = 0;
my $M = 0;
foreach my $id ( @A ){
	
	my $gc = GC( $genome{$id});
	my $string = $genome{$id};
	$M  += length($string);
	$string =~ tr/N//d;
	$N += length($string);
	print "$id\t".length( $genome{$id})."\t".$gc."\t".length($string)."\n";
	
}
print STDERR "$ARGV[0]\t$N\t$M\n";



############################################################
#
# Subroutines
# 
############################################################
sub countStr{
	my ($query, $string) = @_;
	my $counts = 0;
	for( my $i = 0; $i < length( $string ) - length( $query ) +1 ; $i++ ){
		my $subs = substr( $string, $i, length( $query ));	
		$counts++ if( $subs eq $query);
	}
	return $counts;
}


sub revert{
    my ($dna) = @_ ;
    $dna =~ tr/ACGTacgt/TGCAtgca/;
    
    if(0){
    #print "reverse: $dna\n";
    my @bases = split( //, $dna );
    my @revers;
    for( my $i = $#bases ; $i >= 0 ; $i-- ){
	push(@revers, @bases[ $i ]);
    }
    $"=""; # changes the format of the next statment
    return "@revers";}
    return reverse($dna);
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

sub GC{
        my ($dna) = @_;
        $dna =~ s/N//g;
        my $n = length( $dna );
        $dna =~ s/A//g; 
        $dna =~ s/T//g; 
        return 0.5 if( $n == 0 );
        return sprintf( "%.3f", length($dna)/$n );
        
}

