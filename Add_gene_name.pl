#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 3) {
    print STDERR "Usage: Add_gene_name.pl ref.gtf inf outf\n";
    exit(0);
}

my ($ref, $inf, $outf) = @ARGV;

my %hash; # save the gene name

open(IN, $ref) or die "Cannot open $ref\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[8]=~/gene_id "([^"]+)".+gene_name "([^"]+)"/){
        $hash{$1}=$2;
    }
}
close IN;

open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
    next if $_=~/^#/;
    chomp;
    my @info=split(/\t/, $_);
    if($info[0]=~/Row\.names/){
      my $tmp=join"\t",@info[1..$#info];
      print OUT "gene_id\tgene_name\t$tmp\taverageCPM_EXT\taverageCPM_RC\n";
    }else{
      
      my $ave1=sprintf("%.2f", ($info[4]+$info[5]+$info[6])/3);
      my $ave2=sprintf("%.2f", ($info[7]+$info[8]+$info[9])/3);
      unless($ave1>0){
        $ave1=0;
      }
      unless($ave2>0){
        $ave2=0;
      }
      for(my $i=4; $i<=$#info; $i++){
        if($info[$i]>0){
          $info[$i]=sprintf("%.2f", $info[$i]);
        }
      }
      my $tmp=join"\t",@info[1..$#info];
      if(exists $hash{$info[0]}){
        print OUT "$info[0]\t$hash{$info[0]}\t$tmp\t$ave1\t$ave2\n";
      }else{
        print OUT "$info[0]\t$info[0]\t$tmp\t$ave1\t$ave2\n";
      }
    }
}
close IN;
close OUT;
