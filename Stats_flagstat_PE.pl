#!/usr/bin/perl -w
use strict;
use warnings;
if(@ARGV != 2) {
    print STDERR "Usage: Stats_flagstat_PE.pl flatstat_dir outf_stats\n";
    exit(0);
}

my ($indir, $outf)= @ARGV;

open(OUT, ">$outf") or die $!;
print OUT "Sample name\tTotal reads\tMapped reads\t%\tProperly paired mapped reads\t%\n";
opendir(INDIR, $indir) or die "Cannot open dir $indir\n";
my @files = readdir(INDIR);

@files = sort {lc($a) cmp lc($b)} @files;

foreach my $stat (@files) {
		if($stat =~ /(.+)\.flagstat$/){
			my $sn=$1;
			my $pp;
			my $map;
			my $total;
			print OUT "$sn\t";
			open(IN, "$indir/$stat") or die "Cannot open file $indir/$stat\n";
			while(<IN>){
				if($_=~/^(\d+) .+ paired in sequencing/){
					$total=$1;
					print OUT "$1\t";
				}elsif($_=~/^(\d+) .+ properly paired \(([^\s]+)/){
					$pp = "$1\t$2";
				}elsif($_=~/^(\d+) .+with itself and mate mapped/){
					$map=$1;
				}elsif($_=~/^(\d+) .+singletons /){
					$map+=$1;
					my $ratio=sprintf("%.2f%%", $map/$total*100);
					print OUT "$map\t$ratio\t$pp\n";
				}

			}
			close IN;
		}
}
closedir INDIR;
close OUT;
