#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV !=3) {
    print STDERR "Usage: blacklist_remove_bed.pl blacklist.bed inf outf\n";
    exit(0);
}

my ($black, $inf, $outf)=@ARGV;

my %hash;
open(IN, $black) or die "Cannot open $black\n";
while(<IN>){
	chomp;
	my @info=split(/\t/, $_);
	for(my $i=$info[1]+1; $i<=$info[2]; $i++){
		$hash{$info[0]}{$i}=1;
	}
}
close IN;

my $count=0;
open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
	chomp;
	if($_=~/^Chr|^chrom/){
		print OUT "$_\n";
	}else{
		my @info=split(/\t/, $_);
		my $ss=$info[1]+1;
		my $ee=$info[2];
		my $mid=int(($ss+$ee)/2);
		if(!exists $hash{$info[0]}{$ss} && !exists $hash{$info[0]}{$ee} && !exists $hash{$info[0]}{$mid}){
			print OUT "$_\n";
		}else{
			$count++;
		}
	}
}
close OUT;
close IN;
print STDERR "$count peaks were removed because they overlaps with blacklist regions.\n";