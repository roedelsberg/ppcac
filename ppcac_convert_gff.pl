#!/usr/bin/perl

use strict;

die "program exonerate_output\n" if(@ARGV < 1 );

my $name;
my $file = shift( @ARGV );
my $mapping = read_hash(  shift( @ARGV ) );
my $sequence_count;

open(IN, $file);
while(<IN>){
	my ($chr,$source,$class, $start, $end, $score, $strand, $offset, $comment ) = split( /\t+/ );
	my $hash = parseComment( $comment );
	#$chr= "Contig35";
	$source = "est2genome";
	my $general = "$chr\t$source\t$class\t$start\t$end\t$score\t$strand\t$offset";
	if( $class eq "gene" ){	    
	    $name = $hash->{"sequence"};
	    if( exists $mapping->{$name} ){
		$name = $mapping->{ $name };
#		print STDERR $name ." replaced\n";	
	}
	    if( exists $sequence_count->{$name} ){
		my $number = $sequence_count->{$name} + 1;
		$sequence_count->{$name}++;
		$name = $name . "_loc$number";	
	#	##print STDERR "new $name\n";
	#	
	    }
	    else{
		my $number = 1;
		$sequence_count->{$name} = $number;
                $sequence_count->{$name}++;
		$name = $name . "_loc$number";
	    }
	    print "$general\tID=".$name.";Name=".$name."\n";
	    
	}	
	 if( $class eq "exon"){
#	if( $class eq "exon" || $class eq "utr5" || $class eq "utr3" ){
		$general =~ s/$class/CDS/g;
		print "$general\tParent=".$name."\n";
	}



}
close(IN);




sub parseComment{
	my ($text) = @_;
	my @A = split( /\s*;\s*/, $text );
	my $hash;
	foreach my $a( @A ){
		my ($key,$value) = split( / /, $a );
		chomp $value;
		$hash->{$key} = $value;
		
	}
	
	return $hash;
}


sub read_hash{
	my ($file) = @_;
	my $mapping;
	open(IN2, $file);
	while(<IN2>){
		my ($query, $target) = split;
		$mapping->{ $query} = $target unless ( $target eq "" || $query eq "");
	}
	close( IN2 );
	return $mapping;

}
