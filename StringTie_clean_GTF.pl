#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 3) {
    print STDERR "Usage: StringTie_clean_GTF.pl ref.gtf inf outf\n";
    exit(0);
}

my ($ref, $inf, $outf) = @ARGV;

my %gene;
my %transcript;
open(IN, $ref) or die "Cannot open $ref\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    if($info[8]=~/gene_id "([^"]+)"/){
        my $id=$1;
        $gene{$id}=1;
    }
    if($info[8]=~/transcript_id "([^"]+)"/){
        my $id=$1;
        $transcript{$id}=1;
    }
}
close IN;

open(IN, $inf) or die "Cannot open $inf\n";
open(OUT, ">$outf") or die "Cannot open $outf\n";
while(<IN>){
    next if $_=~/^#/;
    my @info=split(/\t/, $_);
    my $id1;
    my $id2;
    if($info[8]=~/gene_id "([^"]+)"/){
        $id1=$1;
    }
    if($info[8]=~/transcript_id "([^"]+)"/){
        $id2=$1;
    }
    if(exists $gene{$id1} && exists $transcript{$id2}){
        print OUT "$_";
    }
}
close IN;
close OUT;
