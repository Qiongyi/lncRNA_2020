
####################################################################################################################
############                        1. lncRNA Capture-Seq data analysis                                 ############
####################################################################################################################

### 1) use cutadapt (v1.17) to clip low-quality nucleotides and adaptor sequences
# Bases lower than a defined Phred quality threshold (default: 20) at the 3′ end were trimmed off from each read using cutadapt (http://code.google.com/p/cutadapt/). Next, known Illumina primers and adaptor sequences were clipped off from each read by cutadapt, which computes sensitive semi-global alignments of all the reads against all the primer/adaptor sequences, allowing gapped and mismatched alignments.

pipeline_cutadapt.pl sample_info.txt /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/ori /illumina/Data/others/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/cutadapt


### 2) mapping against the mouse genome (mm10) using HISAT2 (v2.1.0); convert “SAM” files to “BAM” files, remove duplicate reads, sort and index the “BAM” files. To avoid the artefact signals potentially introduced by misalignments, we only kept properly PE aligned reads with mapping quality at least 20 for downstream analyses.

pipeline_alignment.pl hisat2 config.txt /illumina/reference/mm10/HISAT2_index/mm10 ./cutadapt ./HISAT2


### 3) mapping stats
pipeline_alignment_stats.pl mapping_stats.xls ./cutadapt ./HISAT2


### 4) run StringTie
# 4.1) perform reference-guided transcriptome assembly for each sample
for i in "EXT1" "EXT2" "EXT3" "RC1" "RC2" "RC3"
do
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie ./HISAT2/$i.rmdup.Q20.sort.bam -p 2 -o $i.gtf -e -b ./StringTie/$i.ctab -G /illumina/reference/GENCODE/gencode.vM25.annotation.gtf
done

# 4.2) generate a non-redundant set of transcripts using the StringTie merge mode
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie --merge -G /clusterdata/uqqzhao/reference/GENCODE/gencode.vM25.annotation.gtf -o merge_with_genecode.gtf EXT1.gtf EXT2.gtf EXT3.gtf RC1.gtf RC2.gtf RC3.gtf &

# remove known protein-coding transcripts based on GENCODE annotation
Remove_known_coders_GTF.pl /clusterdata/uqqzhao/reference/GENCODE/gencode.vM25.annotation.gtf merge_with_genecode.gtf merge.RemoveKnownCoders.gtf

# add gene entries
gtf_addGeneEnties.pl merge.RemoveKnownCoders.gtf merge.gtf

# 4.3) quantitate the transcript-level expression for each sample, with the option of “-e -G merged.gtf”
for i in "EXT1" "EXT2" "EXT3" "RC1" "RC2" "RC3"
do
/clusterdata/uqqzhao/tools/StringTie/stringtie-2.1.4.Linux_x86_64/stringtie ./HISAT2/$i.rmdup.Q20.sort.bam -p 2 -o $i.eG.gtf -e -b $i.ctab -G merge.gtf &
done

# clean StringTie GTF (remove transcripts that are not listed in "merge.gtf" - note this is a potential bug that sometimes exists in the StringTie output file.)
for i in "EXT1" "EXT2" "EXT3" "RC1" "RC2" "RC3"
do
StringTie_clean_GTF.pl merge.gtf $i.eG.gtf $i.clean.gtf
done


### 5) get count matrix using the prepDE.py script
prepDE.py -l 140 -g gene_count_matrix.csv -t transcript_count_matrix.csv --input=sample_list.txt


### 6) run ballgown
# see the "ballgown_command.R" script
# we then focused on transcripts with FPKM>1 in at least one group
awk '$1=="transcriptIDs" || ($NF+$(NF-1)+$(NF-2))>3 || (($NF-4)+$(NF-5)+$(NF-6))>3' Nucleus.lncRNACapture.ballgown.xls > Nucleus.lncRNACapture.FPKM.xls


### 7) assign lncRNAs to Proximal/Intragenic/Extragenic categories
lncRNA_overlay_gtf_ballgown.pl 10000 /clusterdata/uqqzhao/reference/mm10/mm10_20190124_gene_region.gtf merged.gtf Nucleus.lncRNACapture.FPKM.xls Nucleus.lncRNACapture.genomic.xls &


### 8) add transcript length, the gene type and transcipt type if known, and edit the column IDs (remove ".ctab" etc.)
Wei_polish1.pl merged.gtf /clusterdata/uqqzhao/reference/GENCODE/gencode.vM25.annotation.gtf Nucleus.lncRNACapture.genomic.xls Nucleus.lncRNACapture.genomic.polish.xls





####################################################################################################################
############                           2.  ATAC-Seq data analysis                                       ############
####################################################################################################################

indir=/illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab

### 1) cut adaptor
pipeline_cutadapt.pl sample_info.txt $indir/ori/ $indir/cutadapt


### 2) mapping
# mapping using BWA, mm10
pipeline_alignment.pl bwa config.txt /illumina/reference/mm10/mm10.fasta $indir/cutadapt $indir/BWA


### 3) get mapping stats
indir=/illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab/BWA
Stats_flagstat_PE.pl $indir/tmp bwa_mapping_stats.xls
Stats_flagstat_PE.pl $indir/BAM bwa_mapping_stats_Q20.xls

# drop sample Ep4 due to extremely low sequencing yields.

### 4) prepare bed files for macs2 peak calling
# (using one sample as an example)

# bam2bed
/clusterdata/uqqzhao/tools/bedtools/bedtools2.27.1/bin/bedtools bamtobed -i /illumina/Data/others/202009_Xiang_ATAC_FACS_BredyLab/BWA/BAM/Ep1.rmdup.Q20.sort.bam > Ep1.rmdup.bed

# shift + 4 bp and − 5 bp for positive and negative strand, respectively and remove reads from the mitochondrial genome
ATAC_bed_shift.pl Ep1.rmdup.bed Ep1.final.bed


### 5) macs2 call peaks
# (using one sample as an example)
macs2 callpeak -t Ep1.final.bed -n Ep1 --shift -75 --extsize 150 --nomodel -B --SPMR -g mm --keep-dup all


### 6) merge multiple replicates and select peaks supported by all replicates for each group
/clusterdata/uqqzhao/tools/bedtools/bedtools2.27.1/bin/multiIntersectBed -header -i Ep1_peaks.narrowPeak Ep3_peaks.narrowPeak > Epos.multiIntersectBed.bed
/clusterdata/uqqzhao/tools/bedtools/bedtools2.27.1/bin/multiIntersectBed -header -i Rp1_peaks.narrowPeak Rp2_peaks.narrowPeak Rp3_peaks.narrowPeak > Rpos.multiIntersectBed.bed

awk '$4==2{print $1"\t"$2"\t"$3}' Epos.multiIntersectBed.bed > ./Epos.con.bed
awk '$4==3{print $1"\t"$2"\t"$3}' Rpos.multiIntersectBed.bed > ./Rpos.con.bed

### 7) remove peaks in mm10 blacklist
blacklist_remove_bed.pl /clusterdata/uqqzhao/reference/blacklist/Blacklist-2.0/lists/mm10-blacklist.v2.bed Epos.con.bed Epos.removeBlacklist.bed
blacklist_remove_bed.pl /clusterdata/uqqzhao/reference/blacklist/Blacklist-2.0/lists/mm10-blacklist.v2.bed Rpos.con.bed Rpos.removeBlacklist.bed



####################################################################################################################
############         3. customised analyses that combine the results of multiple datasets              ############
####################################################################################################################


### 1) overlay analysis between lncRNAs and Arc+ RC & Arc+ EXT ATAC peaks
lncRNA_overlay_bed_OnlyExon_transcript.pl ATAC_RCpos Rpos.removeBlacklist.bed merged.gtf Nucleus.lncRNACapture.genomic.polish.xls tmp1
lncRNA_overlay_bed_OnlyExon_transcript.pl ATAC_EXTpos Epos.removeBlacklist.bed merged.gtf tmp1 Nucleus.lncRNACapture.genomic.ATAC.xls


### 2) overlay analysis between lncRNAs and CBP dataset
# download "GSM530174_CBP_Millipore_KCl_B1_E120.bedgraph.gz" from GEO GSM530174 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM530174)
# mm9 to mm10 liftover
gunzip GSM530174_CBP_Millipore_KCl_B1_E120.bedgraph.gz
liftOver GSM530174_CBP_Millipore_KCl_B1_E120.bedgraph /clusterdata/uqqzhao/reference/liftOver_chain/mm9ToMm10.over.chain GSM530174_CBP_Millipore_KCl_B1_E120_mm10.bed unMapped.tmp

# overlay analysis
lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 CBP GSM530174_CBP_Millipore_KCl_B1_E120_mm10.bed merged.gtf Nucleus.lncRNACapture.genomic.ATAC.xls Nucleus.lncRNACapture.genomic.ATAC.CBP.xls


### 3) overlay analysis between lncRNAs and H3K4ME1/H3K27AC datasets
# download the ".bw" files from GEO GSE74964 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74964)
# BigWig to bedGraph
for i in "GSM1939127_SHC-01H-CA1-NEU-H3K4ME1-1" "GSM1939128_SHC-01H-CA1-NEU-H3K4ME1-2" "GSM1939159_SHC-01H-CA1-NEU-H3K27AC-1" "GSM1939160_SHC-01H-CA1-NEU-H3K27AC-2"
do
/clusterdata/uqqzhao/tools/UCSC_binary_utilities/bigWigToBedGraph $i.bw $i.mm10.bedGraph &
done

# call peaks using macs2
# /clusterdata/uqqzhao/illumina/201912_Wei_Nuclear_LncRNA_CapSeq_BredyLab/enhancer/bonn_data/
for i in "GSM1939127_SHC-01H-CA1-NEU-H3K4ME1-1" "GSM1939128_SHC-01H-CA1-NEU-H3K4ME1-2" "GSM1939159_SHC-01H-CA1-NEU-H3K27AC-1" "GSM1939160_SHC-01H-CA1-NEU-H3K27AC-2"
do
nohup macs2 bdgpeakcall -i ${i}.mm10.bedGraph --o-prefix ${i} --no-trackline &
done


# merge two replicates and keep peaks that are supported by both replicates
i=GSM1939127_SHC-01H-CA1-NEU-H3K4ME1-1
j=GSM1939128_SHC-01H-CA1-NEU-H3K4ME1-2
N=H3K4ME1_SHC-01H-CA1
sort -k1,1 -k2,2n ${i}_c5.0_l200_g30_peaks.narrowPeak > ${i}_mm10.sort.bed
sort -k1,1 -k2,2n ${j}_c5.0_l200_g30_peaks.narrowPeak > ${j}_mm10.sort.bed
/clusterdata/uqqzhao/tools/bedtools/bedtools2.27.1/bin/multiIntersectBed -header -i ${i}_mm10.sort.bed ${j}_mm10.sort.bed > ${N}_B1B2merged_mm10.tmp
awk '$4==2{print $1"\t"$2"\t"$3}' ${N}_B1B2merged_mm10.tmp > ${N}_B1B2merged_mm10.bed

i=GSM1939159_SHC-01H-CA1-NEU-H3K27AC-1
j=GSM1939160_SHC-01H-CA1-NEU-H3K27AC-2
N=H3K27AC_SHC-01H-CA1
sort -k1,1 -k2,2n ${i}_c5.0_l200_g30_peaks.narrowPeak > ${i}_mm10.sort.bed
sort -k1,1 -k2,2n ${j}_c5.0_l200_g30_peaks.narrowPeak > ${j}_mm10.sort.bed
/clusterdata/uqqzhao/tools/bedtools/bedtools2.27.1/bin/multiIntersectBed -header -i ${i}_mm10.sort.bed ${j}_mm10.sort.bed > ${N}_B1B2merged_mm10.tmp
awk '$4==2{print $1"\t"$2"\t"$3}' ${N}_B1B2merged_mm10.tmp > ${N}_B1B2merged_mm10.bed


# overlay analysis
lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 H3K27ac H3K27AC_SHC-01H-CA1_B1B2merged_mm10.bed merged.gtf Nucleus.lncRNACapture.genomic.ATAC.CBP.xls Bonn.tmp1

lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 H3K4me1 H3K4ME1_SHC-01H-CA1_B1B2merged_mm10.bed merged.gtf Bonn.tmp1 Nucleus.lncRNACapture.genomic.ATAC.CBP.H3K27ac.H3K4me1.xls

# note "Nucleus.lncRNACapture.genomic.ATAC.CBP.H3K27ac.H3K4me1.xls" is the Supplemental Table 1.




####################################################################################################################
############         4. The overlay analysis: real data vs. random data                                 ############
####################################################################################################################

# Please note this analysis was performed to address one of the reviewer's question regarding Figure 1B (The venn diagram): “Are these regions more likely to overlap than similar sized genomic regions sampled by chance?”

# To address this question, we randomized the genomic location of same number of “lncRNAs” (with same size and same number of exons), and then performed the exact same overlay analysis with ATAC-Seq, CBP ChIP-Seq, H3K4me1 ChIP-Seq, and H3K27ac ChIP-Seq peaks. We performed the random relocation of lncRNAs three times independently, and identified 32, 39 and 33 “lncRNAs” that overlapped with all four enhancer markers, respectively. These numbers of overlapping "lncRNAs" sampled by chance are dramatically less than 434 overlapping lncRNAs identified in our study. Three random datasets were stored in "random1_lncRNA", "random2_lncRNA" and "random3_lncRNA" folders.


### 1) Randomized the genomic location of “lncRNAs”. The script will take real genomic location of lncRNAs from the input file "lncRNA_2109.xls", randomly relocate it to the same chromosome, and generate a new GTF file "lncRNA_random.gtf".

lncRNA_randomize_genomic_location.pl mm10.fasta.fai merged.gtf lncRNA_2109.xls lncRNA_random.gtf

### 2) overlay with ATAC-Seq peaks

lncRNA_overlay_bed_OnlyExon_transcript.pl ATAC_EXTpos /clusterdata/uqqzhao/illumina/202009_Xiang_ATAC_FACS_BredyLab/results_using_all_replicates/supported_by_all_replicates/Epos.removeBlacklist.bed lncRNA_random.gtf lncRNA_2109.xls lncRNA_random.ATAC.xls

### 3) overlay with CBP peaks
lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 CBP /clusterdata/uqqzhao/reference/Mouse_enhancer_markers/Greenberg_data/GSM530174_CBP_Millipore_KCl_B1_E120_mm10.bed lncRNA_random.gtf lncRNA_random.ATAC.xls lncRNA_random.ATAC.CBP.xls

### 4) overlay with H3K27ac peaks
lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 H3K27AC /clusterdata/uqqzhao/reference/Mouse_enhancer_markers/Bonn_data/H3K27AC_SHC-01H-CA1_B1B2merged_mm10.bed lncRNA_random.gtf lncRNA_random.ATAC.CBP.xls lncRNA_random.ATAC.CBP.H3K27AC.xls

### 5) overlay with H3K4Me1 peaks
lncRNA_overlay_bed_OnlyExon_transcript_TSS.pl 500 H3K4ME1 /clusterdata/uqqzhao/reference/Mouse_enhancer_markers/Bonn_data/H3K4ME1_SHC-01H-CA1_B1B2merged_mm10.bed lncRNA_random.gtf lncRNA_random.ATAC.CBP.H3K27AC.xls lncRNA_random.ATAC.CBP.H3K27AC.H3K4ME1.xls
