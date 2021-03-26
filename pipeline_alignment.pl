#!/usr/bin/perl -w
use warnings;
use strict;
if(@ARGV != 5) {
    print STDERR "Usage: pipeline_alignment.pl alignment_tool(bowtie2/hisat2/bwa) config.txt ref_index fq_indir outdir\n";
  	exit(0);
}
my ($tool, $config, $ref, $indir, $outdir)=@ARGV;


my %ev;


### to check if the splicesites.txt exists or not in the index folder, if not, you need to supply it in the config.txt file
my $splicesites;
if($tool eq "hisat2"){
	my $dir;
	if($ref=~/(.+)\//){
		$dir=$1;
	}else{
		$dir=".";
	}
	if(-e "$dir/splicesites.txt"){
		$ev{"splicesites"}="$dir/splicesites.txt";
	}else{
		print STDERR "Warning: splicesites.txt is not available in $dir. Please generate splicesites.txt using the below command (example for mouse mm10):\npython /illumina/tools/HISAT2/hisat2-2.1.0/hisat2_extract_splice_sites.py mm10.gtf > $dir/splicesites.txt\n\n";
		print STDERR "\n\nYou need to supply the splicesites file in the config file.\n"; exit;
	}
}


### read the config file for alignment options
open(IN, $config) or die "Cannot open $config\n";
while(<IN>){
	if($_=~/^\w/){
		chomp;
		my @info=split(/:/, $_);
		$ev{$info[0]}=$info[1];
	}
}
close IN;

### output dirs
my $out_qlog = $outdir."/qlog";
my $out_cmd = $outdir."/commands";
my $out_bam = $outdir."/BAM";
my $out_tmp = $outdir."/tmp";

unless(-e $out_qlog){system("mkdir -p $out_qlog");}
unless(-e $out_cmd){system("mkdir -p $out_cmd");}
unless(-e $out_bam){system("mkdir -p $out_bam");}
unless(-e $out_tmp){system("mkdir -p $out_tmp");}

opendir(INDIR, $indir) or die "Cannot open dir $indir\n";
while ( my $fastq = readdir(INDIR) ) {
	if($fastq=~/(.+)_read1\.fastq$/){
		my $sn=$1;
		if($tool eq "bowtie2"){
			&run_bowtie2($sn);
		}elsif($tool eq "hisat2"){
			&run_hisat2($sn);
		}elsif($tool eq "bwa"){
			&run_bwa($sn);
		}else{
			print STDERR "Error: wrong alignment tool: $tool. It should be bowtie2/hisat2/bwa. Please specify the alignment tool and run again\n"; exit;
		}
	}elsif($fastq=~/(.+)_read1\.fastq\.gz$/){
		my $sn=$1;
		if($tool eq "bowtie2"){
			&run_bowtie2_gz($sn);
		}elsif($tool eq "hisat2"){
			&run_hisat2_gz($sn);
		}else{
			print STDERR "Error: wrong alignment tool: $tool. It should be bowtie2/hisat2. Please specify the alignment tool and run again\n"; exit;
		}
	}
}
closedir INDIR;


sub run_bowtie2{
	my $sn=shift;
	my $tool="bowtie2";
	my $outf="$out_cmd/$tool.$sn.sh";
	open(OUT, ">$outf") or die "cannot open $outf\n";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH --job-name=$tool.$sn\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=20g\n";
		print OUT "#SBATCH --cpus-per-task=8\n";
		print OUT "#SBATCH --error=$out_qlog/$sn.$tool.out\n";
		print OUT "#SBATCH --output=$out_qlog/$sn.$tool.out\n";
		my $r1="$indir/$sn"."_read1.fastq";
		my $r2="$indir/$sn"."_read2.fastq";
		my $cmd=qq($ev{"bowtie2"} $ev{"bowtie2_option"} -x $ref -1 $r1 -2 $r2 -S $out_bam/$sn.sam);
		print OUT "$cmd\n";
		close OUT;
		my $jobname="$tool.$sn";
		&submit_job($outf, $jobname);
}

sub run_bowtie2_gz{
	my $sn=shift;
	my $tool="bowtie2";
	my $outf="$out_cmd/$tool.$sn.sh";
	open(OUT, ">$outf") or die "cannot open $outf\n";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH --job-name=$tool.$sn\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=20g\n";
		print OUT "#SBATCH --cpus-per-task=8\n";
		print OUT "#SBATCH --error=$out_qlog/$sn.$tool.out\n";
		print OUT "#SBATCH --output=$out_qlog/$sn.$tool.out\n";
		my $r1="$indir/$sn"."_read1.fastq.gz";
		my $r2="$indir/$sn"."_read2.fastq.gz";
		my $cmd=qq($ev{"bowtie2"} $ev{"bowtie2_option"} -x $ref -1 $r1 -2 $r2 -S $out_bam/$sn.sam);
		print OUT "$cmd\n";
		close OUT;
		my $jobname="$tool.$sn";
		&submit_job($outf, $jobname);
}

sub run_hisat2{
	my $sn=shift;
	my $tool="hisat2";
	my $outf="$out_cmd/$tool.$sn.sh";
	open(OUT, ">$outf") or die "cannot open $outf\n";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH --job-name=$tool.$sn\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=20g\n";
		print OUT "#SBATCH --cpus-per-task=8\n";
		print OUT "#SBATCH --error=$out_qlog/$sn.$tool.out\n";
		print OUT "#SBATCH --output=$out_qlog/$sn.$tool.out\n";
		my $r1="$indir/$sn"."_read1.fastq";
		my $r2="$indir/$sn"."_read2.fastq";
		my $cmd;
		$cmd=qq($ev{"hisat2"} $ev{"hisat2_option"} --known-splicesite-infile $ev{"splicesites"} --novel-splicesite-outfile $outdir/$sn.novelsplicesites.txt -x $ref -1 $r1 -2 $r2 -S $out_tmp/$sn.sam);
		print OUT "$cmd\n";
		
		if(exists $ev{"samtools"}){
			#$cmd=qq($ev{"samtools"} sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
			$cmd=qq(
$ev{"samtools"} fixmate -m $out_tmp/$sn.sam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} sort -o $out_bam/$sn.positionsort.bam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} index $out_bam/$sn.positionsort.bam &
$ev{"samtools"} markdup -r $out_bam/$sn.positionsort.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} view -f2 -q 20 -h -o $out_tmp/$sn.rmdup.Q20.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} sort -o $out_bam/$sn.rmdup.Q20.sort.bam $out_tmp/$sn.rmdup.Q20.bam
$ev{"samtools"} index $out_bam/$sn.rmdup.Q20.sort.bam
$ev{"samtools"} flagstat $out_bam/$sn.positionsort.bam > $out_bam/$sn.positionsort.flagstat
$ev{"samtools"} flagstat $out_bam/$sn.rmdup.Q20.sort.bam > $out_bam/$sn.rmdup.Q20.sort.flagstat);
		}else{
			print STDERR "samtools is not defined in the config file\n"; exit;
			#$cmd=qq(samtools sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
		}
		print OUT "$cmd\n";

		# $cmd=qq($ev{"samtools"}/samtools index $outdir/$sn.namesort.bam);
		# print OUT "$cmd\n";		
		close OUT;
		my $jobname="$tool.$sn";
		&submit_job($outf, $jobname);
}

sub run_hisat2_gz{
	my $sn=shift;
	my $tool="hisat2";
	my $outf="$out_cmd/$tool.$sn.sh";
	open(OUT, ">$outf") or die "cannot open $outf\n";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH --job-name=$tool.$sn\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=20g\n";
		print OUT "#SBATCH --cpus-per-task=8\n";
		print OUT "#SBATCH --error=$out_qlog/$sn.$tool.out\n";
		print OUT "#SBATCH --output=$out_qlog/$sn.$tool.out\n";
		my $r1="$indir/$sn"."_read1.fastq.gz";
		my $r2="$indir/$sn"."_read2.fastq.gz";
		my $cmd;
		$cmd=qq($ev{"hisat2"} $ev{"hisat2_option"} --known-splicesite-infile $ev{"splicesites"} -x $ref -1 $r1 -2 $r2 -S $out_tmp/$sn.sam);
		print OUT "$cmd\n";
		if(exists $ev{"samtools"}){
			#$cmd=qq($ev{"samtools"} sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
			$cmd=qq(
$ev{"samtools"} fixmate -m $out_tmp/$sn.sam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} sort -o $out_tmp/$sn.positionsort.bam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} markdup -r $out_tmp/$sn.positionsort.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} view -f2 -q 20 -h -o $out_tmp/$sn.rmdup.Q20.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} sort -o $out_bam/$sn.rmdup.Q20.sort.bam $out_tmp/$sn.rmdup.Q20.bam
$ev{"samtools"} index $out_bam/$sn.rmdup.Q20.sort.bam
$ev{"samtools"} flagstat $out_bam/$sn.rmdup.Q20.sort.bam > $out_bam/$sn.rmdup.Q20.sort.flagstat);
		}else{
			print STDERR "samtools is not defined in the config file\n"; exit;
			#$cmd=qq(samtools sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
		}
		print OUT "$cmd\n";
		# $cmd=qq($ev{"samtools"}/samtools index $outdir/$sn.namesort.bam);
		# print OUT "$cmd\n";		
		close OUT;
		my $jobname="$tool.$sn";
		&submit_job($outf, $jobname);
}

sub run_bwa{
	my $sn=shift;
	my $tool="bwa";
	my $outf="$out_cmd/$tool.$sn.sh";
	open(OUT, ">$outf") or die "cannot open $outf\n";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH --job-name=$tool.$sn\n";
		print OUT "#SBATCH --nodes=1\n";
		print OUT "#SBATCH --ntasks=1\n";
		print OUT "#SBATCH --mem=20g\n";
		print OUT "#SBATCH --cpus-per-task=8\n";
		print OUT "#SBATCH --error=$out_qlog/$sn.$tool.out\n";
		print OUT "#SBATCH --output=$out_qlog/$sn.$tool.out\n";
		my $r1="$indir/$sn"."_read1.fastq";
		my $r2="$indir/$sn"."_read2.fastq";
		my $cmd;
		$cmd=qq($ev{"bwa"} mem $ev{"bwa_option"} -R "\@RG\\tID:BWAmem\\tSM:$sn\\tPL:Illumina\\tLB:QBI" $ref $r1 $r2 > $out_tmp/$sn.sam);
		print OUT "$cmd\n";
		if(exists $ev{"samtools"}){
			#$cmd=qq($ev{"samtools"} sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
			$cmd=qq(
$ev{"samtools"} flagstat $out_tmp/$sn.sam > $out_tmp/$sn.flagstat &
$ev{"samtools"} fixmate -m $out_tmp/$sn.sam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} sort -m 5G -o $out_tmp/$sn.positionsort.bam $out_tmp/$sn.fixmate.bam
$ev{"samtools"} markdup -r $out_tmp/$sn.positionsort.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} view -f2 -q 20 -h -o $out_tmp/$sn.rmdup.Q20.bam $out_tmp/$sn.rmdup.bam
$ev{"samtools"} sort -o $out_bam/$sn.rmdup.Q20.sort.bam $out_tmp/$sn.rmdup.Q20.bam
$ev{"samtools"} index $out_bam/$sn.rmdup.Q20.sort.bam
$ev{"samtools"} flagstat $out_bam/$sn.rmdup.Q20.sort.bam > $out_bam/$sn.rmdup.Q20.sort.flagstat);

		}else{
			print STDERR "samtools is not defined in the config file\n"; exit;
			#$cmd=qq(samtools sort -n -o $outdir/$sn.namesort.bam $outdir/$sn.sam);
		}
		print OUT "$cmd\n";
		# $cmd=qq($ev{"samtools"}/samtools index $outdir/$sn.namesort.bam);
		# print OUT "$cmd\n";		
		close OUT;
		my $jobname="$tool.$sn";
		&submit_job($outf, $jobname);
}

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
