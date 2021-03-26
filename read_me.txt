
####################################################################################################################
############                        1. lncRNA Capture-Seq data analysis                                 ############
####################################################################################################################

### 1) use cutadapt (v1.17) to clip low-quality nucleotides and adaptor sequences
# Bases lower than a defined Phred quality threshold (default: 20) at the 3′ end were trimmed off from each read using cutadapt (http://code.google.com/p/cutadapt/). Next, known Illumina primers and adaptor sequences were clipped off from each read by cutadapt, which computes sensitive semi-global alignments of all the reads against all the primer/adaptor sequences, allowing gapped and mismatched alignments.

pipeline_cutadapt.pl sample_info.txt /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/ori /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/cutadapt


### 2) mapping against the mouse genome (mm10) using HISAT2 (v2.1.0), SAM to BAM, sort BAM and index the BAM files
indir=/illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab
pipeline_alignment.pl hisat2 config.txt /illumina/reference/mm10/HISAT2_index/mm10 $indir/cutadapt $indir/HISAT2


### 3) mapping stats
pipeline_alignment_stats.pl mapping_stats.xls ./cutadapt ./HISAT2


### 4) run StringTie
for i in "EXT1" "EXT2" "EXT3" "RC4" "RC5" "RC6"
do
nohup /illumina/tools/StringTie/stringtie-2.0.3.Linux_x86_64/stringtie /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/HISAT2/$i.rmdup.Q20.sort.bam -p 2 -o $i.eG.gtf -e -b /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/StringTie/$i.ctab -G /illumina/reference/GENCODE/gencode.vM24.long_noncoding_RNAs.gtf &
done


### 5) get count matrix using the prepDE.py script
/illumina/tools/StringTie/stringtie-1.3.4d/prepDE.py -l 140 -g gene_count_matrix.csv -t transcript_count_matrix.csv --input=sample_list.txt


### 6) run edgeR
# see "edgeR_command.R"

### 7) add RefSeq gene name and calculate the average CPM for RC and EXT
Add_gene_name.pl EXT1.eG.gtf RCvsEXT.edgeR RCvsEXT.edgeR.xls




####################################################################################################################
############                           2.  ATAC-Seq data analysis                                       ############
####################################################################################################################

indir=/illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab

### 1) cut adaptor
pipeline_cutadapt.pl sample_info.txt /illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab/ori/ $indir/cutadapt


### 2) mapping
# mapping using BWA, mm10
# build index
# /illumina/reference/mm10
# /illumina/bwa/bwa-0.7.17/bwa index mm10.fasta
pipeline_alignment.pl bwa config.txt /illumina/reference/mm10/mm10.fasta $indir/cutadapt $indir/BWA


### 3) get mapping stats
indir=/illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab/BWA
Stats_flagstat_PE.pl $indir/tmp bwa_mapping_stats.xls
Stats_flagstat_PE.pl $indir/BAM bwa_mapping_stats_Q20.xls


### 4) prepare bed files for macs2 peak calling
# (using one sample as an example)

# bam2bed
/illumina/tools/bedtools/bedtools2.27.1/bin/bedtools bamtobed -i /illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab/BWA/BAM/Ep1.rmdup.Q20.sort.bam > Ep1.rmdup.bed

# shift + 4 bp and − 5 bp for positive and negative strand, respectively and remove reads from the mitochondrial genome
ATAC_bed_shift.pl Ep1.rmdup.bed Ep1.final.bed


### 5) macs2 call peaks
# (using one sample as an example)
/illumina/tools/python3.6.8_virtual_environment/bin/macs2 callpeak -t Ep1.final.bed -n Ep1 --shift -75 --extsize 150 --nomodel -B --SPMR -g mm --keep-dup all


### 6) merge multiple replicates
/illumina/tools/bedtools/bedtools2.27.1/bin/multiIntersectBed -header -i Ep1_peaks.narrowPeak Ep3_peaks.narrowPeak > Epos.multiIntersectBed.bed
awk '$4==2' Epos.multiIntersectBed.bed > Epos.bed


####################################################################################################################
############         3. customised analyses that combine the results of multiple datasets              ############
####################################################################################################################

### 1) liftover H3K27ac peaks (potential enhancers) from GEO GSE60192
# Download the enchancer region from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60192
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60192/suppl/GSE60192_01_04_B1B2merged_un_H3K27ac_ab4729_p1e-5_peaks_genomewide.bed.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60192/suppl/GSE60192_03_05_B1B2merged_KCl_H3K27ac_ab4729_p1e-5_peaks_genomewide.bed.gz
# liftover from mm9 to mm10
liftOver GSE60192_01_04_B1B2merged_un_H3K27ac_ab4729_p1e-5_peaks_genomewide.bed /illumina/reference/liftOver_chain/mm9ToMm10.over.chain un_H3K27ac_mm10.bed unMapped_un_H3K27ac.txt
liftOver GSE60192_03_05_B1B2merged_KCl_H3K27ac_ab4729_p1e-5_peaks_genomewide.bed /illumina/reference/liftOver_chain/mm9ToMm10.over.chain KCl_H3K27ac_mm10.bed unMapped_KCl_H3K27ac.txt


### 2) extract ATAC-Seq peaks that overlap with H3K27ac peaks
/illumina/tools/bedtools/bedtools2.27.1/bin/bedtools intersect -wa -a Epos.bed -b ../KCl_H3K27ac_mm10.bed -wa > EposATAC.KCl.OverlapPeaks.xls
/illumina/tools/bedtools/bedtools2.27.1/bin/bedtools intersect -wa -a Eneg.bed -b ../un_H3K27ac_mm10.bed -wa > EnegATAC.un.OverlapPeaks.xls

### 3) Overlay lnRNA with both ATAC and H3K27ac peaks 
lncRNA_overlay_bed_OnlyExon.pl un_H3K27ac.ATAC EnegATAC.un.OverlapPeaks.xls /illumina/reference/GENCODE/gencode.vM24.long_noncoding_RNAs.gtf RCvsEXT.edgeR.xls tmp.xls
lncRNA_overlay_bed_OnlyExon.pl KCl_H3K27ac.ATAC EposATAC.KCl.OverlapPeaks.xls /illumina/reference/GENCODE/gencode.vM24.long_noncoding_RNAs.gtf tmp.xls RCvsEXT.edgeR.K27ac.ATAC.xls


### 4) category lncRNAs to  intergenic, intronic, antisense
lncRNA_overlay_gtf_ex.pl 10000 mm10_20190124_gene_region.gtf /illumina/reference/GENCODE/gencode.vM24.long_noncoding_RNAs.gtf RCvsEXT.edgeR.K27ac.ATAC.xls RCvsEXT.edgeR.enhancer.genomic.xls







