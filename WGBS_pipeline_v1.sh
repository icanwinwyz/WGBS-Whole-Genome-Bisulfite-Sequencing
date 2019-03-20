#########################################################################
# File Name: WGBS_pipeline_v1.sh
# Author: Di Wu
# email: di.wu@cshs.org
# Created Time: Dec 2018
#########################################################################

#!/bin/sh
#
#$ -v PATH=/common/genomics-core/anaconda2/bin:/opt/sge/bin:/opt/sge/bin/lx-amd64:/opt/sge/bin:/opt/sge/bin/lx-amd64:/opt/sge/bin:/opt/sge/bin/lx-amd64:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/hpc/scripts:/hpc/apps/python/2.7.15/bin
#$ -l mem_free=90G


display_usage() {
	echo -e "NAME:\n  WGBS pipeline v1."
	echo -e "\nDESCRIPTION:\n   This pipeline integrates pre-processing, mapping, post-processing and methylation calling by Biscuit, PicardTools, Samtools and TrimGalore."
	echo -e "\nUsage:\n   bash WGBS_pipeline_v1.sh <Read1.fastq.gz> <Read2.fastq.gz> <Reference> <tmp_folder>"
    echo "Input options:"
    echo "   -h|--help    show this help"
    echo "Input files:"
    echo "   The first input is Read1.fq.gz."
    echo "   The second input is Read2.fq.gz."
    echo "   The third input is ref.fa."
    echo "   The fourth input is a temporary folder. Without setting an additional temporary folder, the system's TMP directory will be filling up and the program will crash. "
    echo "Output files:"
    echo "   bam folder: includes all bam/bai files after mapping, indexing, sorting, deduplication"
    echo "   fastq folder: includes original fastq files"
    echo "   fastqc folder: includes fastqc reports"
    echo "   Trimgalore folder: includes trimmed fastq files, trimming report and fastqc reports"
    echo "   picrad_metrics folder: includes library metrics estimated by PicardTools"
    echo "   results folder: includes multiqc report and bedgrapgh files for CG sites only and for all Cs"
    echo "   others folder: intermediate files"
    echo "   log folder: includes log files during data processing"

	}
# check whether user had supplied -h or --help . If yes display usage
if [[ ( $1 == "--help" ) ||  ( $1 == "-h" ) ]]
then
	display_usage
	exit 0
fi

module load R/3.4.1
module load samtools

R1=$1
R2=$2
ref=$3
tmp_folder=$4
base=$(basename $R1 "_1.fq.gz")
biscuit_path=/common/genomics-core/apps/biscuit/biscuit
#biscuit_path=/common/yeagern/biscuit
samtools_path=/hpc/apps/samtools/1.6/bin/samtools
picardtools_path=/common/genomics-core/bin/picardtools


############## raw reads QC ###################
echo "raw reads QC"
mkdir ${base}"_fastqc"
/common/genomics-core/anaconda2/bin/fastqc -o ./${base}"_fastqc" -f fastq $R1 $R2
/common/genomics-core/bin/wgbs_invert_dups_check_v0.01.pl ${base}"_fastqc"/${base}.InvertedReadPairDups.metric.txt $R1 $R2 > ${base}"_fastqc"/${base}.InvertedReadPairDups.metric.out


############## Trimming adapter, low quality bases and QC ###################
echo "adapter/quality trimming and QC"
mkdir ${base}"_Trimgalore"
mkdir ${base}"_Trimgalore/fastqc"
/common/genomics-core/apps/TrimGalore-0.4.5/trim_galore --retain_unpaired --path_to_cutadapt /common/genomics-core/anaconda2/bin/cutadapt --illumina --paired $R1 $R2 -o ${base}"_Trimgalore" --fastqc_args "-o ./${base}"_Trimgalore/fastqc""



############## Alignment ##################
echo "Alignment to reference genome"
$biscuit_path align -t 10 $ref ${base}"_Trimgalore"/${base}_1_val_1.fq.gz  ${base}"_Trimgalore"/${base}_2_val_2.fq.gz > ${base}.sam
$samtools_path view -bS ${base}.sam -@ 10 > ${base}".bam"
$picardtools_path SortSam I=${base}.bam O=${base}"_sorted.bam" TMP_DIR=$tmp_folder SORT_ORDER=coordinate
$samtools_path index ${base}"_sorted.bam" > ${base}"_sorted.bam.bai"
$samtools_path flagstat ${base}"_sorted.bam" > ${base}"_sorted.bam_flagstat.txt"

##############BAM file#####################
echo "Mark duplicates"
$picardtools_path AddOrReplaceReadGroups I=${base}"_sorted.bam" O=${base}"_sorted_RG.bam" RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit$count RGSM=20 TMP_DIR=$tmp_folder
$picardtools_path MarkDuplicates I=${base}"_sorted_RG.bam" O=${base}"_sorted_RG_markdups.bam" M=${base}"_marked_dups".txt TMP_DIR=$tmp_folder

#ReOrder BAM file to line up ordering used in reference file (usually karyotypic order instead of numeric)
$picardtools_path  ReorderSam I=${base}"_sorted_RG_markdups.bam" O=${base}"_sorted_RG_markdups.reorder.bam" R=$ref TMP_DIR=$tmp_folder

#metrics will point at reordered file instead of original BAM file now
bam_file=${base}"_sorted_RG_markdups.reorder.bam"

#creates index file that allows fast lookup of data within a BAM file ... must be sorted in coordinate order
$picardtools_path  BuildBamIndex I=$bam_file O=$bam_file".bai" TMP_DIR=$tmp_folder


###############Generating some quality control metrics using picardtools##################
#collects info about library construction, including insert size distribution and read orientation of paired-end libraries
$picardtools_path CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT  MINIMUM_PCT=0.05 I=$bam_file O=$bam_file".mdups.bam.CollectInsertSizeMetrics.metric.txt" H=$bam_file".mdups.bam.CollectInsertSizeMetrics.metric.txt.histogram.pdf" TMP_DIR=$tmp_folder

#quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters (specific to Illumina data
$picardtools_path CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT  I=$bam_file O=$bam_file".mdups.bam.CollectAlignmentSummaryMetrics.metric.txt" R=$ref BS=true TMP_DIR=$tmp_folder

#info about relative proportions of G/G nucleotides
$picardtools_path CollectGcBiasMetrics VALIDATION_STRINGENCY=SILENT  I=$bam_file O=$bam_file".CollectGcBiasMetrics.metric.txt" CHART_OUTPUT=$bam_file".CollectGcBiasMetrics.metric.chart.pdf" S=$bam_file".CollectGcBiasMetrics.metric.summary.txt" BS=true TMP_DIR=$tmp_folder REFERENCE_SEQUENCE=$ref

#determine overall quality of given run
$picardtools_path QualityScoreDistribution VALIDATION_STRINGENCY=SILENT  I=$bam_file O=$bam_file".mdups.bam.QualityScoreDistribution.metric.txt" CHART=$bam_file".mdups.bam.QualityScoreDistribution.metric.chart.pdf" TMP_DIR=$tmp_folder

#typically used on single lane run, but can be used on merged BAM
$picardtools_path MeanQualityByCycle VALIDATION_STRINGENCY=SILENT  I=$bam_file O=$bam_file".mdups.bam.MeanQualityByCycle.metric.txt" CHART=$bam_file".mdups.bam.MeanQualityByCycle.metric.chart.pdf" TMP_DIR=$tmp_folder

#estimates number of unique molecules in sequencing library
$picardtools_path EstimateLibraryComplexity VALIDATION_STRINGENCY=SILENT  I=$bam_file O=$bam_file".mdups.bam.EstimateLibraryComplexity.metric.txt" TMP_DIR=$tmp_folder

#Collect nucleotide distribution by cycle
$picardtools_path CollectBaseDistributionByCycle VALIDATION_STRINGENCY=SILENT I=$bam_file O=$bam_file".mdups.bam.CollectBaseDistributionByCycle.txt" CHART=$bam_file".mdups.bam.CollectBaseDistributionByCycle.chart.pdf" R=$ref


##############pileup and get bedgraph###################
# Methylation extraction using Biscuit pileup
$biscuit_path pileup -r  $ref -i $bam_file -q 10 -o $bam_file".out" -w $bam_file"_pileup_stats.out"

#convert vcf to bedgraph using biscuit and calculating methylation ratio and coverage for CpG sites and/or for all Cs
$biscuit_path vcf2bed -t cg -c $bam_file".out" > ${base}"_CpG.bedgraph"
$biscuit_path vcf2bed -t c -c $bam_file".out" > ${base}"_all_Cs.bedgraph"

############### MultiQC report ######################
## conduct MultiQC in current folder
/common/genomics-core/anaconda2/bin/multiqc .

############### orgnization ######################
mkdir bam
mv *.bam *.bai bam/
mkdir fastq
mv *fq.gz fastq/
mkdir metrics
mv *bam.* *marked_dups* metrics/
mkdir results
mv *bedgraph *html *flagstat* results/
mv  others/
mkdir others
mv *sam *.out others/
mkdir log
mv *.e* *.o* log/


echo "Subject: WGBS pipeline is done for ${base}" | sendmail -v di.wu@cshs.org
