#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 4) {
    print STDERR "Usage: Wei_polish1.pl Nucleus.lncRNACapture.gtf GENCODE.gtf Nucleus.lncRNACapture.genomic.xls Nucleus.lncRNACapture.genomic.polish.xls\n";
    exit(0);
}

my ($gtf, $gencode_gtf, $inf, $outf)=@ARGV;

my %coor; # to record the genomic coordinates;
my %len; # to record the length of the transcript;
open(IN, $gtf) or die $!;
while(<IN>){
	next if $_=~/^#/;
	chomp;
	my @info=split(/\t/, $_);
	if($info[2] eq "transcript" && $info[8]=~/transcript_id "([^"]+)";/){
		$coor{$1}="$info[0]:$info[3]-$info[4]:$info[6]";
	}elsif($info[2] eq "exon" && $info[8]=~/transcript_id "([^"]+)";/){
		$len{$1}+=$info[4]-$info[3]+1;
	}
}
close IN;

my %type;
open(IN, $gencode_gtf) or die $!;
while(<IN>){
	next if $_=~/^#/;
	chomp;
	my @info=split(/\t/, $_);
	if($info[2] eq "transcript" && $info[8]=~/transcript_id "([^"]+)";.+gene_type "([^"]+)"/){
		$type{$1}=$2;
	}elsif($info[2] eq "transcript"){
		print STDERR "It seems there are different patterns in the 9th column for gene_type. Please check.\n$_\n"; exit;
	}
}
close IN;


open(IN, $inf) or die $!;
open(OUT, ">$outf") or die $!;
while(<IN>){
	my @info=split(/\t/, $_);
	my $tmp=join"\t", @info[3..$#info];
	if($_=~/^transcriptIDs/){
		$tmp=~s/\.ctab//g;
		print OUT "$info[0]\t$info[1]\t$info[2]\tgeneType\ttranscriptLength\tgenomic_coordinates(chr:start-end:strand)\t$tmp";
	}elsif($_=~/\w/){
		if(!exists $type{$info[0]}){
			$type{$info[0]}=".";
		}
		print OUT "$info[0]\t$info[1]\t$info[2]\t$type{$info[0]}\t$len{$info[0]}\t$coor{$info[0]}\t$tmp";
	}
}
close IN;
close OUT;
