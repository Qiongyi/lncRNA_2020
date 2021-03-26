#!/usr/bin/perl -w
use warnings;
use strict;
if(@ARGV < 2) {
    print STDERR "Usage: pipeline_alignment_stats.pl outf indir1 indir2 ...\n";
  	exit(0);
}
my ($outf, @indir)=@ARGV;

my %total;
my %hash;
my ($cutadapt, $rrna, $phix, $sequin, $th, $hisat2)=(0, 0, 0, 0, 0, 0);

for my $indir (@indir){
	if($indir=~/cutadapt/i){
		$cutadapt=1;
		&get_cutadapt_stats($indir, "cutadapt");
	}elsif($indir=~/rRNA_bowtie2/i){
		$rrna=1;
		&get_bowtie2_stats($indir, "rRNA");
	}elsif($indir=~/PhiX_bowtie2/i){
		$phix=1;
		&get_bowtie2_stats($indir, "PhiX");
	}elsif($indir=~/Sequins_bowtie2/i){
		$sequin=1;
		&get_bowtie2_stats($indir, "sequins");
	}elsif($indir=~/tyrosine_hydroxylase_bowtie2/i){
		$th=1;
		&get_bowtie2_stats($indir, "tyrosine_hydroxylase");
	}elsif($indir=~/HISAT2/){
		$hisat2=1;
		&get_hisat2_stats($indir);
	}
}

open(OUT, ">$outf") or die "Cannot open $outf\n";
my $header="sample_name\ttotal_reads";

if($cutadapt){
	$header.="\treads_after_cutadapt\treads_after_cutadapt%";
}

if($rrna){
	$header.="\tmap_to_rRNA\tmap_to_rRNA%";
}

if($phix){
	$header.="\tmap_to_PhiX\tmap_to_PhiX%";
}

if($sequin){
	$header.="\tmap_to_sequins\tmap_to_sequins%";
}

if($th){
	$header.="\tmap_to_tyrosine_hydroxylase\tmap_to_tyrosine_hydroxylase%";
}

if($hisat2){
	$header.="\tclean_reads\tmapped_reads\tmapped_reads%\tconcordantly_uniquely_mapped_reads\tconcordantly_uniquely_mapped_reads%";
}

print OUT "$header\n";



foreach my $sn (sort keys %total){
	my $line="$sn\t$total{$sn}";
	if($cutadapt){
		my $type="cutadapt";
		$line.="\t$hash{$sn}{$type}";
	}

	if($rrna){
		my $type="rRNA";
		$line.="\t$hash{$sn}{$type}";
	}
	if($phix){
		my $type="PhiX";
		$line.="\t$hash{$sn}{$type}";
	}
	if($sequin){
		my $type="sequins";
		$line.="\t$hash{$sn}{$type}";
	}
	if($th){
		my $type="tyrosine_hydroxylase";
		$line.="\t$hash{$sn}{$type}";
	}	
	if($hisat2){
		my $type="hisat2";
		$line.="\t$hash{$sn}{$type}";
	}
	print OUT "$line\n";
}
close OUT;


sub get_cutadapt_stats{
	my ($indir, $type)=@_;
	opendir(INDIR, "$indir/qlog") or die "Cannot open dir $indir/qlog\n";
	while(my $file = readdir(INDIR)){
		if($file=~/(.+)\.cutadapt\.out/){
			my $sn=$1;
			my $num;
			open(IN, "$indir/qlog/$file") or die "Cannot open file $indir/qlog/$file\n";
			while(<IN>){
				if($_=~/^Total read pairs processed:\s+([\d,]+)/){
					$num=$1;
					$num=~s/,//g;
					$num=$num*2;
				}elsif($_=~/^Pairs written \(passing filters\):\s+([\d,]+)/){
					my $clean=$1;
					$clean=~s/,//g;
					$clean=$clean*2;
					my $ratio=sprintf("%.4f%%", $clean/$num*100);
					$hash{$sn}{$type}="$clean\t$ratio";
					$total{$sn}=$num;
					last;
				}
			}
			close IN;
		}
	}
	closedir(INDIR);
}

sub get_bowtie2_stats{
	my ($indir, $type)=@_;
	opendir(INDIR, "$indir/qlog") or die "Cannot open dir $indir/qlog\n";
	while(my $file = readdir(INDIR)){
		if($file=~/(.+)\.bowtie2\.out/){
			my $sn=$1;
			my $num;
			open(IN, "$indir/qlog/$file") or die "Cannot open file $indir/qlog/$file\n";
			while(<IN>){
				if($_=~/^(\d+) reads; of these:/){
					$num=$1*2;
				}elsif($_=~/(\d+) \([0-9\.]+\%\) aligned 0 times/){
					my $map=$num-$1;
					my $ratio=sprintf("%.4f%%", $map/$num*100);
					$hash{$sn}{$type}="$map\t$ratio";
					if(!exists $total{$sn}){
						$total{$sn}=$num;
					}
					last;
				}
			}
			close IN;
		}
	}
	closedir(INDIR);
}

sub get_hisat2_stats{
	my $indir=shift;
	my $type="hisat2";
	opendir(INDIR, "$indir/qlog") or die "Cannot open dir $indir/qlog\n";
	while(my $file = readdir(INDIR)){
		if($file=~/(.+)\.hisat2\.out/){
			my $sn=$1;
			my $num;
			my $con_map;
			open(IN, "$indir/qlog/$file") or die "Cannot open file $indir/qlog/$file\n";
			while(<IN>){
				if($_=~/^(\d+) reads; of these:/){
					$num=$1*2;
				}elsif($_=~/(\d+) \([0-9\.]+\%\) aligned concordantly exactly 1 time/){
					$con_map=$1*2;
					#print STDERR "conmap:$con_map\n";
				}elsif($_=~/(\d+) \([0-9\.]+\%\) aligned 0 times/){
					my $map=$num-$1;
					my $ratio=sprintf("%.4f%%", $map/$num*100);
					my $ratio2=sprintf("%.4f%%", $con_map/$num*100);
					$hash{$sn}{$type}="$num\t$map\t$ratio\t$con_map\t$ratio2";
					#print STDERR "map:$map\n";
					if(!exists $total{$sn}){
						$total{$sn}=$num;
					}
					last;
				}
			}
			close IN;
		}
	}
	closedir(INDIR);
}



