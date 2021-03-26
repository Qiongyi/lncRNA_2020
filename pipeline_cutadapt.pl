#!/usr/bin/perl -w
use warnings;
use strict;
if(@ARGV != 3) {
    print STDERR "Usage: pipeline_cutadapt.pl sample_info(space/tab delimited) fq_indir outdir\n";
  	exit(0);
}
my ($sn, $indir, $outdir)=@ARGV;


### output dirs
my $out_qlog = $outdir."/qlog";
my $out_cmd = $outdir."/commands";

unless(-e $out_qlog){system("mkdir -p $out_qlog");}
unless(-e $out_cmd){system("mkdir -p $out_cmd");}

=sample_info.txt
####################
-q 20 --minimum-length=50
AH1_combined    AH1
AH2_combined    AH2
AP1_combined    AP1
AP2_combined    AP2
EC1_combined    EC1
EC2_combined    EC2
OH1_combined    OH1
OH2_combined    OH2
OP1_combined    OP1
OP2_combined    OP2
#####################
=cut


my $option=" ";
open(IN, $sn) or die "cannot open $sn\n";
while(<IN>){
	if($_=~/^(-.+)/){
			$option=$1;
			next;
	}elsif($_=~/^\w/){
		chomp;
		my @info=split(/\s+/, $_);
		my ($r1, $r2);
		# supported file name patterns
		if(-e "$indir/$info[0]"."_R1.fastq.gz"){
			$r1="$indir/$info[0]"."_R1.fastq.gz";
			$r2="$indir/$info[0]"."_R2.fastq.gz";
		}elsif(-e "$indir/$info[0]"."_read1.fastq.gz"){
			$r1="$indir/$info[0]"."_read1.fastq.gz";
			$r2="$indir/$info[0]"."_read2.fastq.gz";			
		}elsif(-e "$indir/$info[0]"."_1.fq.gz"){
			$r1="$indir/$info[0]"."_1.fq.gz";
			$r2="$indir/$info[0]"."_2.fq.gz";
		}elsif(-e "$indir/$info[0]"."_R1.fastq"){
			$r1="$indir/$info[0]"."_R1.fastq";
			$r2="$indir/$info[0]"."_R2.fastq";
		}else{
			print STDERR "Warning: there are no fastq files for sample $info[0] - skipped\n";
			next;
		}

		my $outf="$out_cmd/$info[1].cutadapt.sh";
		open(OUT, ">$outf") or die "cannot open $outf\n";
		print OUT "#!/bin/bash\n";
		print OUT "#SBATCH --job-name=$info[1].cutadapt\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=10g\n";
		print OUT "#SBATCH --cpus-per-task=1\n";
		print OUT "#SBATCH --error=$out_qlog/$info[1].cutadapt.out\n";
		print OUT "#SBATCH --output=$out_qlog/$info[1].cutadapt.out\n";		

		my $dest_read1="$outdir/$info[1]"."_read1.fastq";
		my $dest_read2="$outdir/$info[1]"."_read2.fastq";
###		Libraries using SMARTer® Stranded Total RNA-Seq Kit v2 need to trim 3bp at the begining of read2 using the $option in the first line of sample info:
###		IMPORTANT: When performing paired-end sequencing, the first three nucleotides of the second sequencing read (Read 2) are derived from the Pico v2 SMART Adapter. 
###		These three nucleotides must be trimmed prior to mapping.
		if($option!~/\-q/){
			$option.=" -q 20,20";
		}
		if($option!~/\-\-minimum\-length/){
			$option.=" --minimum-length=36";
		}
		my $cmd=qq(cutadapt -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA $option -o $dest_read1 -p $dest_read2 $r1 $r2);
		print OUT "$cmd\n";
		close OUT;
		# provide the job command and the job name
		&submit_job($outf, "$info[1].cutadapt");
	}
}
close IN;


### submit jobs to wiener server
sub submit_job{
	my ($cmd, $jobname)=@_;
	my $run_cmd=qq(sbatch $cmd);
	my $jobnum = `$run_cmd`;
	if($jobnum=~/Submitted batch job (\d+)/){ $jobnum=$1;}
	else{print STDERR "Error: Job number is not captured:\n$jobnum\n"; exit;}
	print STDERR "Job number $jobnum ($jobname) submitted\n\n";
	return $jobnum;
}
