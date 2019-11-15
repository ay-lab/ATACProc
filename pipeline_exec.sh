#!/bin/bash 

#=================
# main executable script of the ATAC seq pipeline
#=================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================

#=================
# script 1 - when fastq files of paired end read are provided as the input
#=================

genome='/home/sourya/genomes/bowtie2_index/hg19/hg19'
dirdata='/home/sourya/test1/'
inpfile1=$dirdata'fastafiles/001_R1.fastq.gz'
inpfile2=$dirdata'fastafiles/001_R2.fastq.gz'
outdir=$dirdata'Sample_TEST_ATAC'
prefix='001'

`pwd`/bin/pipeline.sh -f $inpfile1 -r $inpfile2 -C `pwd`'/configfile_hg19' -n $prefix -g $genome -d $outdir -t 8 -m "16G" -q 20 -D 1 -O 0

#=================
# script 2 - when fastq files of single end end read are provided as the input
#=================

genome='/home/sourya/genomes/bowtie2_index/hg19/hg19'
dirdata='/home/sourya/test2/'
inpfile=$dirdata'merged_inp.fastq.gz'
outdir=$dirdata'Sample_TEST_ATAC'
prefix='002'

`pwd`/bin/pipeline.sh -f $inpfile -C `pwd`'/configfile_hg19' -n $prefix -g $genome -d $outdir -t 8 -m "16G" -q 20 -D 0 -O 0

#=================
# script 3 - when a BAM file is provided as the input
# here reference genome is not used
# however, -w parameter is used to specify the genome for 
# creating UCSC compatible tracks
#=================

dirdata='/home/sourya/test3/'
inpfile=$dirdata'inp.bam'
outdir=$dirdata'Sample_TEST_ATAC'
prefix='003'

`pwd`/bin/pipeline.sh -f $inpfile -C `pwd`'/configfile_hg19' -n $prefix -d $outdir -t 8 -m "16G" -q 20 -D 1 -O 0 -w "hg19"


