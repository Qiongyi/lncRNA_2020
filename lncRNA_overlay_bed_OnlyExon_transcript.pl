#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 5) {
    print STDERR "Usage: lncRNA_overlay_bed_OnlyExon_transcript.pl column_name H3K27ac_mm10.bed gtf Nucleus.lncRNACapture.genomic.polish.xls outf\n";
    exit(0);
}

my ($column_name, $bed, $gtf, $inf, $outf) = @ARGV;

my %hash;

open(IN, $bed) or die "Cannot open $bed\n";
while(<IN>){
    next if $_=~/^#/;
    next if $_=~/^chrom/;
    my @info=split(/\t/, $_);
    for(my $i=$info[1]+1; $i<=$info[2]; $i++){
        $hash{$info[0]}{$i}=1;
    }
}
close IN;

my %gene;
open(IN, $gtf) or die "Cannot open $gtf\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[2] eq "exon" && $info[8]=~/transcript_id "([^"]+)"/){
        my $gene_id=$1;
        for(my $i=$info[3]; $i<=$info[4]; $i+=20){
            if(exists $hash{$info[0]}{$i}){
              $gene{$gene_id}=1;
              last;
            }
        }
        unless(exists $gene{$gene_id}){
          if(exists $hash{$info[0]}{$info[4]}){
            $gene{$gene_id}=1;
          }
        }
    }
}
close IN;

my $count=0;
my $all=0;
open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
    next if $_=~/^#/;
    chomp;
    my @info=split(/\t/, $_);
    if($info[0]=~/^transcriptIDs/){
      print OUT "$_\t$column_name\n";
    }elsif(exists $gene{$info[0]}){
      print OUT "$_\tT\n";
      $count++;
      $all++;
    }else{
      print OUT "$_\tF\n";
      $all++;
    }
}
close OUT;
close IN;

print STDERR "$count/$all lncRNAs exons overlay with $bed\n\n";

