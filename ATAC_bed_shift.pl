#!/usr/bin/perl -w
use strict;
use warnings;

# shifte + 4 bp and − 5 bp for positive and negative strand, respectively
# remove reads from the mitochondrial genome

if(@ARGV !=2) {
    print STDERR "Usage: ATAC_bed_shift.pl in.bed out.bed\n";
    exit(0);
}
my ($inf, $outf)=@ARGV;
open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
	chomp;
	my @info=split(/\t/,$_);
	if($info[0] =~ /^chrM/){
		next;
	}
	if($info[5] eq "+"){
		$info[1]=$info[1]+4;
		print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\n";
	}elsif($info[5] eq "-"){
		$info[2]=$info[2]-5;
		print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$info[5]\n";			
	}
}
close IN;
close OUT;
