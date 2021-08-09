#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 5) {
    print STDERR "Usage: lncRNA_overlay_gtf_ballgown.pl distance_cutoff(e.g. 10000) ref_gene.gtf lncRNA.gtf RCvsEXT.edgeR.xls outf\n";
    exit(0);
}

my ($cutoff, $ref, $gtf, $inf, $outf) = @ARGV;

# to consider both strands from lncRNA for proximal protein coding genes

### for overlapping with coding genes
my %hash;
open(IN, $ref) or die "Cannot open $ref\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[8]=~/gene_id "([^"]+)"; transcript_id "NM_.+/){
      my $chr=$info[0].$info[6];
      for(my $i=$info[3]; $i<=$info[4]; $i++){
        $hash{$chr}{$i}=$1;
      }
    }
}
close IN;
print STDERR "$ref 1st pass finished!\n";

### for proximal
my %pro; # to record the Proximal sites of reference coding genes (eg. upstream 10kb)
open(IN, $ref) or die "Cannot open $ref\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[8]=~/gene_id "([^"]+)"; transcript_id "NM_.+/){
      my $chr=$info[0].$info[6];
      if($info[6] eq "+"){
        for(my $i=$info[3]-$cutoff; $i<$info[3]; $i++){
          if(!exists $hash{$chr}{$i}){
            $pro{$chr}{$i}=$1;
          }
        }
      }else{
        for(my $i=$info[4]+1; $i<=$info[4]+$cutoff; $i++){
          if(!exists $hash{$chr}{$i}){
            $pro{$chr}{$i}=$1;
          }
        }
      }
    }
}
close IN;
print STDERR "$ref 2nd pass finished!\n";

my %tag; # 1 for overlapping, 2 for proximal same strand with lncRNA, 3 for proximal the other strand;
my %gene;  # for proximal protein coding gene at the same strand of lncRNA
my %antigene; # for proximal protein coding gene at the other strand of lncRNA
my %anti; # for antisense protein coding gene overlapping with the lncRNA
my %atg; # for antisense protein coding gene overlapping with the lncRNA, and the start codon ATG is embeded in the lncRNA exons.
open(IN, $gtf) or die "Cannot open $gtf\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[2] eq "transcript" && $info[8]=~/transcript_id "([^"]+)"/){
        my $gene_id=$1;
        if($info[6] eq "."){$info[6]="+";}
        my $chr=$info[0].$info[6];
        for(my $i=$info[3]; $i<=$info[4]; $i+=200){
            if(exists $hash{$chr}{$i}){
              $tag{$gene_id}=1;
              $gene{$gene_id}=$hash{$chr}{$i};
              last;
            }
        }
        unless(exists $tag{$gene_id}){
          if(exists $hash{$chr}{$info[4]}){
            $tag{$gene_id}=1;
            $gene{$gene_id}=$hash{$chr}{$info[4]};
          }
        }
        if(exists $tag{$gene_id}){
          next;
        }

        # Proximal
        if($info[6] eq "+"){
            if(exists $pro{$chr}{$info[4]}){
              $tag{$gene_id}.=2;
              $gene{$gene_id}=$pro{$chr}{$info[4]};
            }
        }elsif($info[6] eq "-"){
            if(exists $pro{$chr}{$info[3]}){
              $tag{$gene_id}.=2;
              $gene{$gene_id}=$pro{$chr}{$info[3]};
            }
        }
      

        # antisense
        if($info[6] eq "+"){
          $chr=$info[0]."-";
        }else{
          $chr=$info[0]."+";
        }
        for(my $i=$info[3]; $i<=$info[4]; $i+=200){
            if(exists $hash{$chr}{$i}){
              $anti{$gene_id}=$hash{$chr}{$i};
              last;
            }
        }
        unless(exists $anti{$gene_id}){
          if(exists $hash{$chr}{$info[4]}){
            $anti{$gene_id}=$hash{$chr}{$info[4]};
          }
        }

        # Proximal antisense strand
        if($info[6] eq "+"){
            if(exists $pro{$chr}{$info[3]}){
              $tag{$gene_id}.=3;
              $antigene{$gene_id}=$pro{$chr}{$info[3]};
            }
        }elsif($info[6] eq "-"){
            if(exists $pro{$chr}{$info[4]}){
              $tag{$gene_id}.=3;
              $antigene{$gene_id}=$pro{$chr}{$info[4]};
            }
        }
    }
}
close IN;

my $count_proximal=0;
my $count_intragenic=0;
my $count_extragenic=0;
my $count_antisense=0;
my $all=0;
open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
    next if $_=~/^#/;
    chomp;
    my @info=split(/\t/, $_);
    if($info[0]=~/^transcriptIDs/){
      print OUT "$_\tProximal/Intragenic/Extragenic\tProximal/Intragenic_refgene\tProximal_antisense_refgene\tAntisense_refgene\n";
    }else{
      if(!exists $anti{$info[0]}){
        $anti{$info[0]}="-";
      }else{
        $count_antisense++;
      }

      if(exists $tag{$info[0]} && $tag{$info[0]}=~/1/){
        print OUT "$_\tIntragenic\t$gene{$info[0]}\t-\t$anti{$info[0]}\n";
        $count_intragenic++;
        $all++;
      }elsif(exists $tag{$info[0]} && $tag{$info[0]}=~/2/){
        if(!exists $antigene{$info[0]}){
            $antigene{$info[0]}="-";
        }
        print OUT "$_\tProximal\t$gene{$info[0]}\t$antigene{$info[0]}\t$anti{$info[0]}\n";
        $count_proximal++;
        $all++;
      }elsif(exists $tag{$info[0]} && $tag{$info[0]}=~/3/){
        if(!exists $gene{$info[0]}){
            $gene{$info[0]}="-";
        }
        print OUT "$_\tProximal\t$gene{$info[0]}\t$antigene{$info[0]}\t$anti{$info[0]}\n";
        $count_proximal++;
        $all++;
      }else{
        print OUT "$_\tExtragenic\t-\t-\t$anti{$info[0]}\n";
        $count_extragenic++;
        $all++;
      }
    }
}
close OUT;
close IN;
=cut

my $count_proximal=445;
my $all=24330;
my $count_intragenic=7194;
my $count_extragenic=16691;
my $count_antisense=7250;
=cut

print STDERR "Proximal: $count_proximal/$all\t".sprintf("%.2f%%", $count_proximal/$all*100)."\n";

print STDERR "Intragenic: $count_intragenic/$all\t".sprintf("%.2f%%", $count_intragenic/$all*100)."\n";

print STDERR "Extragenic: $count_extragenic/$all\t".sprintf("%.2f%%", $count_extragenic/$all*100)."\n";

print STDERR "Antisense: $count_antisense/$all\t".sprintf("%.2f%%", $count_antisense/$all*100)."\n";


