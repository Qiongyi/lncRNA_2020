#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: lncRNA_randomize_genomic_location.pl chromosome_size(mm10.fasta.fai) GTF inf(lncRNA_list) outf\n";
    exit(0);
}
my ($size, $gtf, $inf, $outf)=@ARGV;


# to record the size of each chromosome
my %size;
open(IN, $size) or die "Cannot open $size\n";
while(<IN>){
	my @chr=split(/\t/, $_);
	$size{$chr[0]}=$chr[1];
}
close IN;

# read the GTF file
my %hash;
open(IN, $gtf) or die "Cannot open $gtf\n";
while(<IN>){
	my @info=split(/\t/,$_);
	if($info[8]=~/transcript_id \"([^"]+)\";/){
		my $t_id=$1;
		$hash{$t_id}.=$_;
	}
}
close IN;

open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
	if($_=~/^MSTR|^ENSMUST/){
		my @info=split(/\t/, $_);
		my @gtf=split(/\n/, $hash{$info[0]});
		my @tt=split(/\t/, $gtf[0]);
		
		my $len=$tt[4]-$tt[3]+1;
		my $max=$size{$tt[0]}-$len;
		my $start_locus=int(rand($max));

		my $shift=$start_locus-$tt[3];
		foreach my $line (@gtf){
			my @tmp=split(/\t/, $line);
			$tmp[3]=$tmp[3]+$shift;
			$tmp[4]=$tmp[4]+$shift;
			my $new_line=join"\t",@tmp;
			print OUT "$new_line\n";
		}
	}
}
close IN;
close OUT;
