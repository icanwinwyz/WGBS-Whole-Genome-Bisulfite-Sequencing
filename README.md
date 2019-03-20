# WGBS(Whole-Genome-Bisulfite-Sequencing)

## NAME:
  WGBS pipeline v1.

## DESCRIPTION:
   This pipeline integrates pre-processing, mapping, post-processing and methylation calling by [Biscuit](https://github.com/zwdzwd/biscuit), [PicardTools](https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary), [Samtools](http://samtools.sourceforge.net) and [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

## Usage:
### Index reference for alignment
```bash
$ biscuit index Reference.fa
```
### Creates a sequence dictionary for a reference sequence  
It is reauired by many processing and analysis steps of Picrdtools
```bash
$ java -jar picard.jar CreateSequenceDictionary REFERENCE=Reference.fa OUTPUT=Reference.dict
```
### Run WGBS pipeline
```bash
$ WGBS_pipeline_v1.sh <Read1.fastq.gz> <Read2.fastq.gz> <Reference.fa> <tmp_folder>
```
### Input options:  
   -h|--help    show this help
   
### Input files:
   The first input is Read1.fq.gz  
   The second input is Read2.fq.gz  
   The third input is ref.fa  
   The fourth input is a temporary folder. Without setting an additional temporary folder, the defacult system's TMP directory will be filling up and the program will crash.  
   
### Output files:  
**bam folder**: includes all bam/bai files after mapping, indexing, sorting, deduplication   
**fastq folder**: includes original fastq files  
**fastqc folder**: includes fastqc reports   
**Trimgalore folder**: includes trimmed fastq files, trimming report and fastqc reports  
**picrad_metrics folder**: includes library metrics estimated by PicardTools  
**results folder**: includes multiqc report and bedgrapgh files for CG sites only and for all Cs  
**others folder**: includes intermediate files  
**log folder**: includes log files during data processing  
