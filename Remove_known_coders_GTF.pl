#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 3) {
    print STDERR "Usage: Remove_known_coders_GTF.pl gencode.vM25.annotation.gtf merge.gtf output.gtf\n";
    exit(0);
}

my ($gtf, $inf, $outf)=@ARGV;

my %hash;
open(IN, $gtf) or die $!;
while(<IN>){
	next if $_=~/^#/;
	chomp;
	my @info=split(/\t/, $_);
	if($info[2] eq "transcript" && $info[8]=~/transcript_id "([^"]+)";.+transcript_type "([^"]+)";/){
		my ($tid, $type)=($1, $2);
		if($type=~/protein_coding/){
			$hash{$tid}=2;
		}else{
			$hash{$tid}=1;
		}
	}
}
close IN;

open(IN, $inf) or die $!;
open(OUT, ">$outf") or die $!;
while(<IN>){
	next if $_=~/^#/;
	my @info=split(/\t/, $_);
	if($info[8]=~/transcript_id "([^"]+)";/){
		my $tid=$1;
		if(exists $hash{$tid} && $hash{$tid}==2){
			next;
		}else{
			print OUT $_;
		}
	}
}
close IN;
close OUT;
