#!/usr/bin/perl -w
use strict;
use warnings;

##
##      Program name: gtf_addGeneTranscriptEnties.pl
##      Author: Qiongyi Zhao 
##      Email: q.zhao@uq.edu.au
##      This script is used to add gene entries (column 3) to the GTF file (e.g. from StringTie merge) for some tools, such as deepTools
## 		
##
###############################################################################################################

if(@ARGV != 2) {
    print STDERR "Usage: gtf_addGeneEnties.pl gtf outf\n";
    exit(0);
}

my ($gtf, $outf)= @ARGV;

my %gg;
open(IN, $gtf) or die "Cannot open $gtf\n"; 
while(<IN>){
	next if $_=~/^#/;
	my @info=split(/\t/,$_);
	if($info[2] eq "exon" && $info[8]=~/gene_id "([^"]+)".+ transcript_id "([^"]+)";/){
		my ($g_id, $t_id)=($1, $2);
#		print STDERR "~~~$info[3]~~~$info[4]~~~\n"; sleep 1;
#		$tt{$t_id}.="$info[3]\t$info[4]\t";
		$gg{$g_id}.="$info[3]\t$info[4]\t";
	}
}
close IN;


my %hash2; # for gene
open(OUT, ">$outf") or die $!;
open(IN, $gtf) or die "Cannot open $gtf\n"; 
while(<IN>){
	next if $_=~/^#/;
	my @info=split(/\t/,$_);
	if($info[2] eq "transcript" && $info[8]=~/gene_id "([^"]+)".+ transcript_id "([^"]+)";/){
		my ($g_id, $t_id)=($1, $2);
		if(!exists $hash2{$g_id}){
			$gg{$g_id}=~s/\t$//;
			my @unsorted=split(/\t/, $gg{$g_id});
			my @sorted = sort { $a <=> $b } @unsorted;
			my $anno;
			if($info[8]=~/(gene_id[^;]+)/){
				$anno.="$1; ";
			}
			if($info[8]=~/(gene_name[^;]+)/){
				$anno.="$1; ";
			}			
			print OUT "$info[0]\t$info[1]\tgene\t$sorted[0]\t$sorted[$#sorted]\t$info[5]\t$info[6]\t$info[7]\t$anno\n";
			$hash2{$g_id}=1;
		}
		print OUT $_;
	}else{
		print OUT $_;
	}
}
close IN;
close OUT;
