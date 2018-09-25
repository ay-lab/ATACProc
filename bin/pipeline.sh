#!/bin/bash

#=================================
# this program denotes a sample pipeline for ATAC-seq data
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: 

1) When raw fastq files (paired end) are provided as input:
./pipeline.sh [-h] [-C configfile] [-f FASTQ1] [-r FASTQ2] [-n PREFIX] [-g BOWTIE2_GENOME] [-a ALIGNVALIDMAX] [-l MAXFRAGLEN] [-d OUTDIR] [-c CONTROLBAM] [-m MAX_MEM] [-O overwrite]
2) When single end fastq files is provided as input
./pipeline.sh [-h] [-C configfile] [-f FASTQ1] [-n PREFIX] [-g BOWTIE2_GENOME] [-a ALIGNVALIDMAX] [-l MAXFRAGLEN] [-d OUTDIR] [-c CONTROLBAM] [-m MAX_MEM] [-O overwrite]
3) When an existing alignment is provided as input
./pipeline.sh [-h]  [-f FASTQ1] [-d OUTDIR] [-c CONTROLBAM] [-w BigWigGenome] [-m MAX_MEM] [-n PREFIX] [-O overwrite]

Example:
1) ./pipeline.sh -f R1.fq.gz -r R2.fq.gz -n 'demo' -g '/home/sample_ATAC/bowtie2_index/hg19/hg19' -a 4 -m "4G" -l 1000 -d '/home/sample_ATAC' -O 1
2) ./pipeline.sh -f inp_align.bam -n 'demo' -a 4 -m "4G" -l 1000 -d '/home/sample_ATAC' -w 'hg19' -O 1

Options:    

  -- required:
  	-C  ConfigFile		 Name of the configuration file for executing ATAC-seq pipeline
	-f  FASTQ1           R1 of pair-end sequencing data  [.fq|.gz|.bz2]. 
						 Or, even an aligned genome (.bam) file can be provided.
	-r  FASTQ2           R2 of pair-end sequencing data [.fq|.gz|.bz2]. If not provided, 
						 the input is assumed to be single ended.
	-g  BOWTIE2_GENOME   Bowtie2 indexed reference genome.
	-n  PREFIX           Prefix of output files.
	-d  OutDir 			 Set the output directory which will contain all the results
	-c  CONTROLBAM		 Control file used to call MACS2. 
						 The control file can be a single file, or can be a 
						 collection of files. 
						 Even it may not be specified at all, in which case 
						 MACS2 operates without any control.
					 	 Control file can be either in BAM or in tag align (.gz) format.
					 	 If a set of control files are provided, user needs 
					 	 to ensure that all of the files follow the same format.
 	-w 	BigWigGenome	 The reference genome which is used to convert BAM file 
 						 to a BigWig file.
 						 If -g option is enabled (i.e. the Bowtie2 index genome is 
						 provided) then this option is not required
 						 as the reference genome will be derived from the genome name 
 						 provided as the Bowtie2 index.
 						 Otherwise, this option needs to be filled with the 
 						 reference genome string (such as 'hg19')
 	-D  DEBUG_TXT		 this flag signifies whether the read count and other 
 						 statistical paramters are computed or not. 
 						 Can be 1 or 0. If 1, the statistics is generated and 
 						 stored in respective files.
	-O 	Overwrite		 this boolean option signifies whether existing output 
						 files would be overwritten (1) or not (0).
						 Default = 0
  -- optional:
	-t  INT              Set number of sorting, Bowtie2 mapping threads [8].
	-m  MAX_MEM          Set max memory of duplication removal [8G].
	-a  ALIGNVALIDMAX	 Set the number of (max) valid alignments which will be searched
	-l  MAXFRAGLEN 		 Set the maximum fragment length to be used for Bowtie2 alignment
	-p  PEAKCALLGENOMESIZE genome size parameter for MACS2 peak calling 
						   ("hs", "mm", "ce", "dm": default is "hs")
	-q  MAPQ_THR		 Quality value threshold, below which the mapped reads 
						 are to be removed (Default 30)

EOF
}

#============================
# this flag signifies whether the read count and other statistical paramters are computed or not
# if 1, then the statistics is generated
DEBUG_TXT=1

# a few bowtie2 related parameters
# the multimapping flag - at most, this no of valid alignments will be searched
MULTIMAP=4

# maximum fragment length considered for valid paired end alignments
MAX_FRAG_LEN=2000

# the number of threads used for execution
THREADS=1

# set the default output directory
OutDir=`pwd`

# maximum memory allotted
MAX_MEM="8G"

# default prefix 
PREFIX=""

# threshold of mapq quality
MAPQ_THR=30

# genome size parameter employed for MACS specific peak calling
PEAKCALLGENOMESIZE="hs"

# # default control bam file
# CONTROLBAM=""

# default fasta formatted input files
FASTQ1=""
FASTQ2=""

# P value used to compute the MACS2 peak
MACS2_P_Val=1e-3

#MACS2_Q_Val=0.05

# Q value thresholds
# -log10(0.05)
Q_Thr1=1.3
# -log10(0.01)
Q_Thr2=2

# Bowtie2 index file of the reference genome
GENOME=""

# Reference genome file to construct the bigwig file from the input BAM file
# if Bowtie2 index file of the reference genome (i.e. the value GENOME) is not provided, 
# this option needs to be provided
BigWigGenome=""

# this boolean option signifies whether existing output 
# files would be overwritten (1) or not (0).
# Default = 0
Overwrite=0

#============================
# reference packages / executables
#============================
# R package installed - executable
RPackageExec='which Rscript'

# python executable
PythonExec=`which python`

# picard_exec='/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

# conversion from one or more BAM files to tag align (.gz) files
TagAlignExec='TagAlign.sh'

# executable to convert the sorted bam file to the bigwig format
BigWigCreateExec='bam_to_bigwig.sh'

# PeakOverlapCode='/home/sourya/proj/Analysis_Scripts_Util/Peak_Intersect.r'

#====================
# utilities to convert the MACS2 detected peak in to big bed format
# useful for displaying in UCSC genome browser
#====================

# # file (SQL) required to convert the narrowPeak file to the bigBed format
# NarrowPeakASFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/narrowPeak.as'
# BigNarrowPeakASFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/bigNarrowPeak.as'
# BroadPeakASFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/broadPeak.as'

# # chromosome size information of the reference hg19 genome
# Refhg19ChrSize='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/chrom_hg19.sizes'

# # chromosome size information of the reference hg38 genome
# Refhg38ChrSize='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/hg38.chrom.sizes'

# # chromosome size information of the reference mm9 genome
# Refmm9ChrSize='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/chrom_mm9.sizes'

# # chromosome size information of the reference mm10 genome
# Refmm10ChrSize='/mnt/BioAdHoc/Groups/vd-vijay/sourya/genomes/chrsize/mm10.chrom.sizes'

#============================

while getopts "C:f:r:n:g:t:m:d:a:l:p:c:q:w:D:O:" opt;
do
	case "$opt" in
		C) ConfigFile=$OPTARG;;
		f) FASTQ1=$OPTARG;;
		r) FASTQ2=$OPTARG;;
		c) CONTROLBAM+=($OPTARG);;	# one or more control files can be provided
		n) PREFIX=$OPTARG;;
		g) GENOME=$OPTARG;;
		d) OutDir=$OPTARG;;
		t) THREADS=$OPTARG;;
		m) MAX_MEM=$OPTARG;;
		a) MULTIMAP=$OPTARG;;
		l) MAX_FRAG_LEN=$OPTARG;;
		p) PEAKCALLGENOMESIZE=$OPTARG;;
		q) MAPQ_THR=$OPTARG;;
		w) BigWigGenome=$OPTARG;;
		D) DEBUG_TXT=$OPTARG;;
		O) Overwrite=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $ConfigFile ]]; then
	echo 'Configuration file is not provided - check the option -C - quit !! '
	exit 1
fi

# error check on no input
if [[ -z $FASTQ1 ]] && [[ -z $FASTQ2 ]]; then
	echo 'User did not provide any fastq file nor any BAM file - quit !!'
	exit 1
fi

# check the BigWig file reference genome information
if [[ -z $BigWigGenome ]]; then
	if [[ -z $GENOME ]]; then
		echo 'User neither provided Bowtie2 reference genome, nor BigWig reference genome - quit !! '
		exit 1
	else
		# #infer the BigWig reference genome from the Bowtie2 reference index genome
		# basically cut the string after the last '/' character
		BigWigGenome="${GENOME##*/}"
		echo 'Inferring BigWig reference genome from the Bowtie2 reference genome index'
		echo 'Genome: '$GENOME
		echo 'BigWigGenome: '$BigWigGenome
	fi
fi

# extract the base name of the bowtie2 reference genome
if [[ ! -z $GENOME ]]; then
	BOWTIE2_GENOME=$(basename ${GENOME})
else
	BOWTIE2_GENOME=$BigWigGenome
fi

# creating the directory which will contain the ATAC seq output
# remove the trailing '/' character, if any
if [[ ${OutDir: -1} == "/" ]]; then
	OutDir=${OutDir%?}
fi
mkdir -p $OutDir
echo '**** OutDir of ATAC Seq pipeline: '$OutDir

echo "BOWTIE2_GENOME: "$BOWTIE2_GENOME
echo "BigWigGenome: "$BigWigGenome


#=================================
# parse the configuration file
#=================================
echo -e "\n ================ Parsing input configuration file ================="

# separator used in the config file
IFS="="
while read -r name value
do
	param=$name
	paramval=${value//\"/}
	if [[ -n $param ]]; then
		if [[ $param != \#* ]]; then
			echo -e "Content of $param is $paramval"
			# if [ $param == "sppexec" ]; then
			# 	sppexec=$paramval
			# fi
			if [ $param == "picardexec" ]; then
				picard_exec=$paramval
			fi
			if [ $param == "HOMERPath" ]; then
				HOMERPath=$paramval
			fi
			if [ $param == "DeepToolsDir" ]; then
				DeepToolsDir=$paramval
			fi		
			if [ $param == "RPackageExec" ]; then
				RPackageExec=$paramval
			fi			
			if [ $param == "NarrowPeakASFile" ]; then
				NarrowPeakASFile=$paramval
			fi			
			if [ $param == "BigNarrowPeakASFile" ]; then
				BigNarrowPeakASFile=$paramval
			fi			
			if [ $param == "BroadPeakASFile" ]; then
				BroadPeakASFile=$paramval
			fi			
			if [ $param == "RefChrSizeFile" ]; then
				RefChrSizeFile=$paramval
			fi
			if [ $param == "RefChrFastaFile" ]; then
				RefChrFastaFile=$paramval
			fi
			if [ $param == "RefChrAnnotFile" ]; then
				RefChrAnnotFile=$paramval
			fi
		fi
	fi
done < $ConfigFile

# if [[ -z sppexec ]]; then
# 	echo 'SPP executable path (from the package phantompeakqualtools by Anshul Kundaje et al.) is not provided - check the configuration file - quit !! '
# 	exit 1
# fi

if [[ -z $picard_exec ]]; then
	echo 'Picard tool executable path is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $HOMERPath ]]; then
	echo 'HOMER executable path is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $DeepToolsDir ]]; then
	echo 'Deeptools executable path is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $RPackageExec ]]; then
	echo 'R executable is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $NarrowPeakASFile ]]; then
	echo 'File to convert narrowPeak to BigBed (NarrowPeakASFile) is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $BigNarrowPeakASFile ]]; then
	echo 'File to convert BignarrowPeak to BigBed (BigNarrowPeakASFile) is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $BroadPeakASFile ]]; then
	echo 'File to convert BroadPeak to BigBed (BroadPeakASFile) is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $RefChrSizeFile ]]; then
	echo 'Reference chromosome size file is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $RefChrFastaFile ]]; then
	echo 'Reference chromosome fasta file is not provided - check the configuration file - quit !! '
	exit 1
fi

if [[ -z $RefChrAnnotFile ]]; then
	echo 'Reference chromosome UCSC annotation file (.gtf format) is not provided - check the configuration file - quit !! '
	exit 1
fi


# code in HOMER package
# which annotates peaks according to different genomic segments
HOMERPeakAnnotExec=$HOMERPath'/annotatePeaks.pl'

# #----------------------------------
# # select the chromosome size file according to the 
# # reference chromosome
# if [[ $BOWTIE2_GENOME = *"hg19"* || $BigWigGenome = *"hg19"* ]]; then
# 	RefChrSizeFile=$Refhg19ChrSize
# elif [[ $BOWTIE2_GENOME = *"hg38"* || $BigWigGenome = *"hg38"* ]]; then
# 	RefChrSizeFile=$Refhg38ChrSize
# elif [[ $BOWTIE2_GENOME = *"mm9"* || $BigWigGenome = *"mm9"* ]]; then
# 	RefChrSizeFile=$Refmm9ChrSize
# elif [[ $BOWTIE2_GENOME = *"mm10"* || $BigWigGenome = *"mm10"* ]]; then
# 	RefChrSizeFile=$Refmm10ChrSize
# fi

#----------------------------------
# select the chromosome size file, reference fasta file, reference UCSC annotation file,
# and the effective genome size (for normalized signal track generation) according to the 
# reference chromosome
# the effective genome size parameters are obtained from
# http://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
if [[ $BOWTIE2_GENOME = *"hg19"* || $BigWigGenome = *"hg19"* ]]; then
	# RefChrSizeFile=$Refhg19ChrSize
	# RefChrFastaFile=$hg19_fastafile
	# RefChrAnnotFile=$hg19_ucsc_annotationfile
	EGS=2864785220
elif [[ $BOWTIE2_GENOME = *"hg38"* || $BigWigGenome = *"hg38"* ]]; then
	# RefChrSizeFile=$Refhg38ChrSize
	# RefChrFastaFile=$hg38_fastafile
	# RefChrAnnotFile=$hg38_ucsc_annotationfile
	EGS=2913022398
elif [[ $BOWTIE2_GENOME = *"mm9"* || $BigWigGenome = *"mm9"* ]]; then
	# RefChrSizeFile=$Refmm9ChrSize
	# RefChrFastaFile=$mm9_fastafile
	# RefChrAnnotFile=$mm9_ucsc_annotationfile
	EGS=2620345972
elif [[ $BOWTIE2_GENOME = *"mm10"* || $BigWigGenome = *"mm10"* ]]; then
	# RefChrSizeFile=$Refmm10ChrSize
	# RefChrFastaFile=$mm10_fastafile
	# RefChrAnnotFile=$mm10_ucsc_annotationfile
	EGS=2652783500
else
	EGS=0
fi

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

# re-assignment of the executable
# conversion from one or more BAM files to tag align (.gz) files
TagAlignExec=`pwd`'/'$TagAlignExec

# re-assignment of the executable
# executable to convert the sorted bam file to the bigwig format
BigWigCreateExec=`pwd`'/'$BigWigCreateExec

#----------------------------------
# output directory for bowtie / alignment
bowtie2_outdir=$OutDir'/Alignment_MAPQ'$MAPQ_THR'/'
mkdir -p $bowtie2_outdir
bowtie2_BAM_prefix=$bowtie2_outdir$PREFIX'.align.sort.MAPQ'$MAPQ_THR
bowtie2_logfile=$bowtie2_outdir$PREFIX'.align.log'

# temporary file names
bowtie2_init_align_samfile=$bowtie2_outdir'Bowtie2_Init_Align.sam'
del_mitch_read_bamfile=$bowtie2_outdir'Bowtie2_del_Mitch.bam'
uniq_mapped_read_bamfile=$bowtie2_outdir'UniqMappedRead.bam'
del_random_read_bamfile=$bowtie2_outdir'Bowtie2_del_Random.bam'


#=========================
# check whether fastq / bam files are provided 
# and whether it is a single end or paired end
# if it is in fastq format then we have to start from the alignment
# else if it is BAM format, we can skip the alignment
#=========================

filebase1=$(basename "$FASTQ1")
if [[ $filebase1 =~ \.fastq.gz$ || $filebase1 =~ \.fq.gz$ || $filebase1 =~ \.fastq$ || $filebase1 =~ \.fq$ ]]; then	
	# boolean flag signifying the use of fastq files
	fastq_input=1
	if [ -z "$FASTQ2" ]; then
	    echo "Single end read fastq file is provided as input"
	    paired_read=0
	else
	    echo "Paired end read fastq file is provided as input"
	    paired_read=1
	fi
elif [[ $filebase1 =~ \.bam$ ]]; then
	# boolean flag signifying the use of fastq files
	fastq_input=0
	# paired read information is obtained by the following command
	# 1 means the input is paired end
	paired_read=`samtools view -c -f 1 $FASTQ1`
else
	echo "Input file is not fasta or bam - error."; 
	exit 1;
fi


#=========================
# process the input alignment(s)
# according to their format
# and also whether they are single end or paired end
#=========================
if [ $fastq_input == 1 ]; then

	# if two fastq files are provided (paired end input)
	# then trim the adapters of the fastq file
	if [ $paired_read == 1 ]; then

		trim_adapter_dir=$OutDir'/Trim_Adapter/'
		mkdir -p $trim_adapter_dir

		# extract only the filenames (excluding all the directory information)
		filebase1=$(basename "$FASTQ1")
		filebase2=$(basename "$FASTQ2")

		# check the extensions of input files
		# accordingly, set the output file names
		if [[ $filebase1 =~ \.fastq.gz$ ]]; then
			file1=${filebase1%.fastq.gz}
			file2=${filebase2%.fastq.gz}
		fi

		if [[ $filebase1 =~ \.fq.gz$ ]]; then
			file1=${filebase1%.fq.gz}
			file2=${filebase2%.fq.gz}
		fi

		if [[ $filebase1 =~ \.fastq$ ]]; then
			file1=${filebase1%.fastq}
			file2=${filebase2%.fastq}
		fi

		if [[ $filebase1 =~ \.fq$ ]]; then
			file1=${filebase1%.fq}
			file2=${filebase2%.fq}
		fi

		if [ ! -f $trim_adapter_dir$file1'.trim.fastq.gz' ] || [ ! -f $trim_adapter_dir$file2'.trim.fastq.gz' ]; then
			# call the adapter trimming function
			# specify input fastq files
			# also mention the output directory to store these files
			$PythonExec ../src/trim_adapters.py -a $FASTQ1 -b $FASTQ2 -d $trim_adapter_dir

			# comment - sourya
			# trim_file1=$file1'.trim.fastq'
			# trim_file2=$file2'.trim.fastq'
			# gzip $trim_file1
			# gzip $trim_file2
			# trim_file1=$trim_file1'.gz'
			# trim_file2=$trim_file2'.gz'

			# mv $trim_file1 $trim_adapter_dir
			# mv $trim_file2 $trim_adapter_dir
			# trim_file1=$trim_adapter_dir$trim_file1
			# trim_file2=$trim_adapter_dir$trim_file2

			# add - sourya
			trim_file1=$trim_adapter_dir$file1'.trim.fastq'
			trim_file2=$trim_adapter_dir$file2'.trim.fastq'
			gzip $trim_file1
			gzip $trim_file2
			trim_file1=$trim_file1'.gz'
			trim_file2=$trim_file2'.gz'

		else
			trim_file1=$trim_adapter_dir$file1'.trim.fastq.gz'
			trim_file2=$trim_adapter_dir$file2'.trim.fastq.gz'
		fi

		echo 'trim_file1: '$trim_file1
		echo 'trim_file2: '$trim_file2
	
	fi 	# end paired end read condition

	# the --mm option is used for memory mapped I/O: fast parallel execution
	# output of Bowtie is a sam file
	# get the uniquely mapped reads by the flag 1804 - indicates discarding any improper mapping
	# source: https://github.com/kundajelab/training_camp/wiki/2.3.-Processing-the-aligned-reads
	# also remove the mitochondrial chromosome (indicated by chrM)

	if [ $DEBUG_TXT == 1 ]; then
		# number of reads of fastq files
		if [ $paired_read == 1 ]; then
			fastq1_read1=`zcat $FASTQ1 | wc -l`
			fastq1_read=`bc <<< "$fastq1_read1 / 4"`
		else 
			fastq1_read1=`zcat $trim_file1 | wc -l`
			fastq1_read=`bc <<< "$fastq1_read1 / 4"`
		fi
	fi

	# bowtie2 alignment
	if [[ ! -f $bowtie2_init_align_samfile || $Overwrite == 1 ]]; then
		if [ $paired_read == 0 ]; then
			bowtie2 -k $MULTIMAP --mm --threads $THREADS -X $MAX_FRAG_LEN -x $GENOME -U $FASTQ1 2>$bowtie2_logfile > $bowtie2_init_align_samfile
		else
			bowtie2 -k $MULTIMAP --mm --threads $THREADS -X $MAX_FRAG_LEN -x $GENOME -1 $trim_file1 -2 $trim_file2 2>$bowtie2_logfile > $bowtie2_init_align_samfile
		fi
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# total no of reads - mapped / unmapped
		nread_tot=`samtools view -S $bowtie2_init_align_samfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# number of mappable reads (excluding the unmapped reads)
	# computed after removing the mitochondrial reads
	if [ $DEBUG_TXT == 1 ]; then
		nread_mappable=`samtools view -S -F 0x04 $bowtie2_init_align_samfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# delete the random stuffs (say chr1_... etc.. and chrY also)
	if [[ ! -f $del_random_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -Sh $bowtie2_init_align_samfile | awk '(substr($1, 1, 1)=="@") || (( $3 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M)$/ ) && ( ( $7 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M)$/ ) || ($7=="=") || ($7=="*") ))' - | samtools view -bhS - > $del_random_read_bamfile
	fi

	# count the number of reads remaining
	if [ $DEBUG_TXT == 1 ]; then
		nread_del_random_stuff=`samtools view $del_random_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# now delete the mitochondrial reads 
	if [[ ! -f $del_mitch_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -h $del_random_read_bamfile | sed '/chrM/d;/random/d;/chrUn/d;/chrY/d' - | samtools view -Shb - > $del_mitch_read_bamfile
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# number of reads after mitochondrial read delete
		nread_del_mit=`samtools view $del_mitch_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# create BAM file consisting of the uniquely mapped reads
	# the flag 1804 = read unmapped, mate unmapped, not primary alignment, read quality low, PCR / optical duplicate
	if [[ ! -f $uniq_mapped_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -hb -F 1804 $del_mitch_read_bamfile > $uniq_mapped_read_bamfile
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# number of uniquely mapped reads
		uniq_mapped_read=`samtools view $uniq_mapped_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# perform quality based thresholding and sorting operation
	if [[ ! -f $bowtie2_BAM_prefix'.bam' || $Overwrite == 1 ]]; then
		# old code - samtools version 1.3
		# samtools view -hb -q $MAPQ_THR $uniq_mapped_read_bamfile | samtools sort - $bowtie2_BAM_prefix
		# new code - samtools version 1.6
		samtools view -hb -q $MAPQ_THR $uniq_mapped_read_bamfile | samtools sort -o $bowtie2_BAM_prefix'.bam' - 
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# count the number of reads after quality thresholding
		nread_qual=`samtools view $bowtie2_BAM_prefix'.bam' | cut -f 1 | sort | uniq | wc -l`
	fi

else

	# convert the input alignment file to match the mapping quality constraints
	# also remove the unnecessary chromosome information
	if [ $DEBUG_TXT == 1 ]; then
		# total no of reads - mapped / unmapped
		nread_tot=`samtools view $FASTQ1 | cut -f 1 | sort | uniq | wc -l`
	fi

	# number of mappable reads (excluding the unmapped reads)
	# computed after removing the mitochondrial reads
	if [ $DEBUG_TXT == 1 ]; then
		nread_mappable=`samtools view -F 0x04 $FASTQ1 | cut -f 1 | sort | uniq | wc -l`
	fi

	# delete the random stuffs (say chr1_... etc.. and chrY also)
	if [[ ! -f $del_random_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -h $FASTQ1 | awk '(substr($1, 1, 1)=="@") || (( $3 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M)$/ ) && ( ( $7 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M)$/ ) || ($7=="=") || ($7=="*") ))' - | samtools view -bhS - > $del_random_read_bamfile
	fi

	# count the number of reads remaining
	if [ $DEBUG_TXT == 1 ]; then
		nread_del_random_stuff=`samtools view $del_random_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# now delete the mitochondrial reads 
	if [[ ! -f $del_mitch_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -h $del_random_read_bamfile | sed '/chrM/d;/random/d;/chrUn/d;/chrY/d' - | samtools view -Shb - > $del_mitch_read_bamfile
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# number of reads after mitochondrial read delete
		nread_del_mit=`samtools view $del_mitch_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# create BAM file consisting of the uniquely mapped reads
	# the flag 1804 = read unmapped, mate unmapped, not primary alignment, read quality low, PCR / optical duplicate
	if [[ ! -f $uniq_mapped_read_bamfile || $Overwrite == 1 ]]; then
		samtools view -hb -F 1804 $del_mitch_read_bamfile > $uniq_mapped_read_bamfile
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# number of uniquely mapped reads
		uniq_mapped_read=`samtools view $uniq_mapped_read_bamfile | cut -f 1 | sort | uniq | wc -l`
	fi

	# perform quality based thresholding and sorting operation
	if [[ ! -f $bowtie2_BAM_prefix'.bam' || $Overwrite == 1 ]]; then
		# old code - samtools version 1.3
		# samtools view -hb -q $MAPQ_THR $uniq_mapped_read_bamfile | samtools sort - $bowtie2_BAM_prefix
		# new code - samtools version 1.6
		samtools view -hb -q $MAPQ_THR $uniq_mapped_read_bamfile | samtools sort -o $bowtie2_BAM_prefix'.bam' - 
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# count the number of reads after quality thresholding
		nread_qual=`samtools view $bowtie2_BAM_prefix'.bam' | cut -f 1 | sort | uniq | wc -l`
	fi

fi

#=========================
# for the generated alignment (with quality thresholding)
# now remove the duplicates, and generate quality statistics
#=========================

# index the sorted file
# provided the index file either does not exist
# or has a modification time earlier than the bam file itself
if [[ ! -f $bowtie2_BAM_prefix'.bam.bai' ]]; then
	samtools index $bowtie2_BAM_prefix'.bam'
elif [[ $bowtie2_BAM_prefix'.bam.bai' -ot $bowtie2_BAM_prefix'.bam' ]]; then
	# here -ot corresponds to "older than"
	samtools index $bowtie2_BAM_prefix'.bam'
fi

# now remove any PCR duplicates using Picard tool
if [[ ! -f $bowtie2_BAM_prefix'.rmdup.bam' || $Overwrite == 1 ]]; then
	java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates INPUT=$bowtie2_BAM_prefix'.bam' OUTPUT=$bowtie2_BAM_prefix'.rmdup.bam' ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$bowtie2_BAM_prefix'.picard_metrics.txt'
fi

if [ $DEBUG_TXT == 1 ]; then
	# count the number of reads after removing mitochondrial reads
	nread_rmdup=`samtools view $bowtie2_BAM_prefix'.rmdup.bam' | cut -f 1 | sort | uniq | wc -l`
fi

if [ $DEBUG_TXT == 1 ]; then
	# note down the read count
	out_readcount_file=$OutDir'/Read_Count_Stat.txt'
	# if [ ! -f $out_readcount_file ]; then
		if [ $fastq_input == 1 ]; then
			echo -e 'TotalRawReads \t TotRead \t NumMappableRead \t nread_del_random \t nread_del_mit \t UniqMappedRead \t ReadQualThr \t rmDupRead' > $out_readcount_file
			echo -e '\n'$fastq1_read'\t'$nread_tot'\t'$nread_mappable'\t'$nread_del_random_stuff'\t'$nread_del_mit'\t'$uniq_mapped_read'\t'$nread_qual'\t'$nread_rmdup >> $out_readcount_file
		else
			echo -e 'TotRead \t NumMappableRead \t nread_del_random \t nread_del_mit \t UniqMappedRead \t ReadQualThr \t rmDupRead' > $out_readcount_file
			echo -e '\n'$nread_tot'\t'$nread_mappable'\t'$nread_del_random_stuff'\t'$nread_del_mit'\t'$uniq_mapped_read'\t'$nread_qual'\t'$nread_rmdup >> $out_readcount_file
		fi
	# fi
fi

#======================
# we use the picard tool to measure the fragment length distribution
# Note: applicable for paired end input data
#======================
if [ $paired_read == 1 ]; then

	picard_insert_metricfile=$bowtie2_outdir$PREFIX'.Picard_insert_size_metrics.MAPQ'$MAPQ_THR'.txt'
	picard_insert_histfile=$bowtie2_outdir$PREFIX'.Picard_insert_size_histogram.MAPQ'$MAPQ_THR'.pdf'

	# picard tool based histogram
	if [ ! -f $picard_insert_metricfile ] || [ ! -f $picard_insert_histfile ]; then 
		java -Xmx$MAX_MEM -jar $picard_exec CollectInsertSizeMetrics I=$bowtie2_BAM_prefix'.rmdup.bam' O=$picard_insert_metricfile H=$picard_insert_histfile M=0.5
		# plotting the ATAC seq representative plot
		# normalized read count vs bp distance
		$PythonExec ../src/PlotSample.py -I $picard_insert_metricfile
	fi
fi

#============================
# create a shifted tagalign file (TN5)
# as mentioned in the following links:
# https://www.biostars.org/p/187204/
# http://seqanswers.com/forums/archive/index.php/t-59219.html
# https://github.com/kundajelab/atac_dnase_pipelines
#============================
# amount of shift done in forward and reverse strand
# for producing the TN5 shift
# this is with respect to the original paper 
# This tagalign file is later used for MACS2 peak calling
fwdshft=4
revshft=5

Shifted_TagAlign_File=$bowtie2_BAM_prefix'.TN5.tagAlign.gz'

if [[ ! -f $Shifted_TagAlign_File || $Overwrite == 1 ]]; then
	# call the utility function to shift the bam file
	# with respect to forward and reverse strands
	# Note: the input files (which can be multiple in the target function)
	# are placed at the last part of a command
	$TagAlignExec -N 0 -f $fwdshft -r $revshft -O $Shifted_TagAlign_File -q $MAPQ_THR -I $bowtie2_BAM_prefix'.rmdup.bam'
fi

# =============================
# get the NRF / library complexity
# check the following link: 
# https://www.encodeproject.org/data-standards/terms/#library
# =============================
# first create a tag align file (different from the earlier created shifted tag align file)
# which will contain the information of genomic locations mapped 
curr_tagalign_file=$OutDir'/temp_tagAlign.gz'
if [[ ! -f $curr_tagalign_file || $Overwrite == 1 ]]; then
	# this file is created from the uniquely mapped reads (before duplicate removal)
	$TagAlignExec -N 0 -f 0 -r 0 -O $curr_tagalign_file -q $MAPQ_THR -I $bowtie2_BAM_prefix'.bam'
fi

# generate a temporary file which will contain the genomic positions
# and the count of reads which map uniquely to these positions
temp_NRF_PBC_file=$OutDir'/temp_NRF_PBC.bed'
if [[ ! -f $temp_NRF_PBC_file || $Overwrite == 1 ]]; then
	zcat $curr_tagalign_file | cut -f1-3 | sort -k1,1 -k2,2n -k3,3n | awk -v OFS='\t' '{a[$1" "$2" "$3]+=1}END{for (i in a){print i,a[i]}}' - > $temp_NRF_PBC_file
fi

if [ $DEBUG_TXT == 1 ]; then

	# file to contain the NRF and other mappability statistics
	# for quality analysis
	out_NRFFile=$OutDir'/out_NRF_MAPQ'$MAPQ_THR'.txt'

	# number of distinct genome position where some read maps uniquely
	uniqgenomepos=`cat $temp_NRF_PBC_file | wc -l`

	# NRF (Non redundant fraction) value
	# Number of distinct uniquely mapping reads (after removing duplicates) / total number of reads
	# other definition is number of distinct genome position for uniquely mapped reads / number of uniquely mapped reads
	# we follow this second definition (ENCODE PAPER)
	NRFval=`bc <<< "scale=3; ($uniqgenomepos * 1.0) / $uniq_mapped_read"`

	# number of genomic locations where exactly one read maps uniquely
	M1=`awk '$4==1' $temp_NRF_PBC_file | wc -l`

	# PCR Bottlenecking Coefficient 1 (PBC1) is computed by considering the 
	# number of genomic locations where exactly one read maps uniquely
	# dividing by the number of distinct genomic locations to which some read maps uniquely
	PBC1=`bc <<< "scale=3; ($M1 * 1.0) / $uniqgenomepos"`

	# number of genomic locations where exactly two reads map uniquely
	M2=`awk '$4==2' $temp_NRF_PBC_file | wc -l`

	# PCR Bottlenecking Coefficient 2 (PBC2)
	# ratio of the following quantities
	if [[ $M2 -gt 0 ]]; then
		PBC2=`bc <<< "scale=3; ($M1 * 1.0) / $M2"`
	else
		PBC2=0
	fi

	# write the NRF statistics
	echo -e 'Unique_Mapped_Read \t Unique_Genome_Pos \t NRF \t M1 \t M2 \t PBC1 \t PBC2' > $out_NRFFile
	echo -e $uniq_mapped_read'\t'$uniqgenomepos'\t'$NRFval'\t'$M1'\t'$M2'\t'$PBC1'\t'$PBC2 >> $out_NRFFile

fi

#======================
# convert the alignment file to the bigwig data format
# for track visualization
#======================
BigWig_outdir=$OutDir'/Out_BigWig'
mkdir -p $BigWig_outdir

# we use sorted (before duplicate removal) bam file
if [[ ! -f $BigWig_outdir'/'$PREFIX'.bw' || $Overwrite == 1 ]]; then
	$BigWigCreateExec -I $bowtie2_BAM_prefix'.bam' -g $BigWigGenome -d $BigWig_outdir -n $PREFIX
fi

#======================
# convert the alignment file to the bigwig data format
# here tracks are normalized with respect to the coverage
# deeptools routines are used
#======================
BigWig_outdir1=$OutDir'/Out_BigWig_NormCov'
mkdir -p $BigWig_outdir1

# we use sorted (before duplicate removal) bam file
if [[ ! -f $BigWig_outdir1'/'$PREFIX'_NormCov.bw' || $Overwrite == 1 ]]; then
	# call the deeptools routine
	if [[ $EGS -gt 0 ]]; then
		$DeepToolsDir'/bamCoverage' -b $bowtie2_BAM_prefix'.bam' -o $BigWig_outdir1'/'$PREFIX'_NormCov.bw' -of bigwig -bs 10 --effectiveGenomeSize $EGS --normalizeUsing RPKM -e 200
	fi
fi

#================
# call the new ATAC seq quality control script with the bam file
# having duplicates
#================
inpfile=$bowtie2_BAM_prefix'.bam'
outdir=$(dirname "${inpfile}")'/ATAC_Seq_QC_New'
mkdir -p $outdir

if [[ ! -f $inpfile'.bai' ]]; then
	samtools index $inpfile
elif [[ $inpfile'.bai' -ot $inpfile ]]; then
	# here -ot corresponds to "older than"
	samtools index $inpfile
fi

$RPackageExec $script_dir'/ATACSeqQC.r' $inpfile $outdir 0


#===============================
# we also create a tag align formatted bed file of the control sample(s)
# we have to check the count and extension of the control samples
#===============================
nctrl=${#CONTROLBAM[@]}
echo 'number of control samples provided: '$nctrl

# after correctly determining the no of control samples
# combine / convert the control samples to tag align gz format
if [[ $nctrl -gt 0 ]]; then
	ctrlextn=$(basename "${CONTROLBAM[0]}")
	if [[ $ctrlextn =~ \.gz$ && $nctrl == 1 ]]; then
		# here one control file is provided and it is already in the tagalign format
		# so just use this file
		Control_TagAlign_File_IDR="${CONTROLBAM[0]}"
	else
		# here multiple control files are provided
		# they may be either in BAM or in tagalign format
		# they need to be combined 
		ControlDir=$(dirname "${CONTROLBAM[0]}")
		Control_TagAlign_File_IDR=$ControlDir'/Control_TagAlign.bed.gz'
		if [ ! -f $Control_TagAlign_File_IDR ]; then
			# we form the command to convert the files into a single tagalign file
			# the command is formed so as to keep the input files (which can be more than one)
			# at the last part
			arguments=' -N 0 -f '$fwdshft' -r '$revshft' -O '$Control_TagAlign_File_IDR
			for (( i=0; i<${nctrl}; i++ ));
			do 
				arguments=$arguments' -I '${CONTROLBAM[i]}
			done
			$TagAlignExec $arguments
		fi
	fi
fi

# =============================
# Now call the MACS2 package for peak calling
# the command has been finalized by the following :
# https://github.com/ParkerLab/bioinf525#sifting
# https://github.com/taoliu/MACS/issues/145
# https://www.biostars.org/p/207318/
# https://www.biostars.org/p/209592/
# =============================

#=================================
# calling the narrow peak version
#=================================

#================================
# default MACS2 parameter (P value based criterion)
#================================
MACS2_outdir_default=$OutDir'/MACS2_Default_Tag'
if [[ $nctrl -gt 0 ]]; then
	MACS2_outdir_default=$MACS2_outdir_default'_with_Control/'
else
	MACS2_outdir_default=$MACS2_outdir_default'_No_Control/'
fi
mkdir -p $MACS2_outdir_default

#**************
# derive the narrowPeaks
#**************
MACS2PeakOutFile=$MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak'
if [[ ! -f $MACS2PeakOutFile || $Overwrite == 1 ]]; then
	MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2 -p '$MACS2_P_Val' --outdir '$MACS2_outdir_default
	
	# this is an alernate command 
	# only invoked when the above comamnd fails
	MACS2_alternate_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2 -p '$MACS2_P_Val' --nomodel --extsize 147 --outdir '$MACS2_outdir_default

	if [[ $nctrl -gt 0 ]]; then
		# include the control bed file also
		MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
		MACS2_alternate_cmd=$MACS2_alternate_cmd' -c '$Control_TagAlign_File_IDR
	fi
	# execute the command
	# Note: the try-catch module
	{
		eval $MACS2_cmd	
	} || {
		eval $MACS2_alternate_cmd
		echo -e 'Default MACS2 command failed. So, executing MACS2 with the default model parameters (--nomodel --extsize 147)' > $MACS2_outdir_default'/MACS2_exception.log'
	}
fi 	# end MACS2 file existence condition


# filter the peaks using Q or P values
# 9th field stores the Q value (in log10 scale)
# 8th field stores the P value (in log10 scale)

# now filter the peaks using Q value of 0.05 and 0.01 (using log10 scales)
QFilt1File=$MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt'
if [[ ! -f $QFilt1File || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr1 '($9 > t)' $MACS2PeakOutFile > $QFilt1File
fi

QFilt2File=$MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt'
if [[ ! -f $QFilt2File || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr2 '($9 > t)' $MACS2PeakOutFile > $QFilt2File
fi

#================================
# annotate the peaks 
#================================

#####################
# peaks with q-value threshold of 0.01
#####################

TempPeakFile=$MACS2_outdir_default'/TempPeak_Q0.01.bed'

# output directory storing the peak annotations
# according to Homer specifications
PeakAnnotateDir=$MACS2_outdir_default'/Peak_Annotate_Q0.01'
mkdir -p $PeakAnnotateDir

# file containing the summarized annotations of the peak files
# according to the HOMER specifications
OutTextFile=$PeakAnnotateDir'/Out_Summary.log'

cat ${QFilt2File} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}' - > ${TempPeakFile}

# apply HOMER specific peak annotation
${HOMERPeakAnnotExec} ${TempPeakFile} ${RefChrFastaFile} -gtf ${RefChrAnnotFile} 2>${OutTextFile} > ${PeakAnnotateDir}'/Annotated_Peak_Q0.01filt.txt'

# also obtain the distance of peaks 
# from the nearest TSS sites
# using the same HOMER annotation
# check : https://www.biostars.org/p/205576/
${HOMERPeakAnnotExec} ${QFilt2File} ${RefChrFastaFile} -size 6000 -hist 50 -bedGraph $BigWig_outdir'/Inp.Sorted.bedGraph' 2>${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist_Summary.log' > ${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist.txt'

# remove the temporary peak file
rm ${TempPeakFile}

# now call an R script
# which takes 
# 1) the summary file 
# and finds the percentage of different annotations for this peak file
# and also 2) the TSS distance from the peaks
# and plots the distribution
$RPackageExec ../Analysis/PeakAnnotateHomerSummary.r $OutTextFile ${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist.txt'

#####################
# peaks with q-value threshold of 0.05
#####################

TempPeakFile=$MACS2_outdir_default'/TempPeak_Q0.05.bed'

# output directory storing the peak annotations
# according to Homer specifications
PeakAnnotateDir=$MACS2_outdir_default'/Peak_Annotate_Q0.05'
mkdir -p $PeakAnnotateDir

# file containing the summarized annotations of the peak files
# according to the HOMER specifications
OutTextFile=$PeakAnnotateDir'/Out_Summary.log'

cat ${QFilt1File} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}' - > ${TempPeakFile}

# apply HOMER specific peak annotation
${HOMERPeakAnnotExec} ${TempPeakFile} ${RefChrFastaFile} -gtf ${RefChrAnnotFile} 2>${OutTextFile} > ${PeakAnnotateDir}'/Annotated_Peak_Q0.05filt.txt'

# also obtain the distance of peaks 
# from the nearest TSS sites
# using the same HOMER annotation
# check : https://www.biostars.org/p/205576/
${HOMERPeakAnnotExec} ${QFilt1File} ${RefChrFastaFile} -size 6000 -hist 50 -bedGraph $BigWig_outdir'/Inp.Sorted.bedGraph' 2>${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist_Summary.log' > ${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist.txt'

# remove the temporary peak file
rm ${TempPeakFile}

# now call an R script
# which takes 
# 1) the summary file 
# and finds the percentage of different annotations for this peak file
# and also 2) the TSS distance from the peaks
# and plots the distribution
$RPackageExec ../Analysis/PeakAnnotateHomerSummary.r $OutTextFile ${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist.txt'

#**************
# derive the Broad Peaks
#**************
MACS2BroadPeakOutFile=$MACS2_outdir_default$PREFIX'.macs2_Broad_peaks.broadPeak'
if [[ ! -f $MACS2BroadPeakOutFile || $Overwrite == 1 ]]; then
	MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -p '$MACS2_P_Val' --broad --outdir '$MACS2_outdir_default
	
	# this is an alernate command 
	# only invoked when the above comamnd fails
	MACS2_alternate_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -p '$MACS2_P_Val' --broad --nomodel --extsize 147 --outdir '$MACS2_outdir_default

	if [[ $nctrl -gt 0 ]]; then
		# include the control bed file also
		MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
		MACS2_alternate_cmd=$MACS2_alternate_cmd' -c '$Control_TagAlign_File_IDR
	fi
	# execute the command
	# Note: the try-catch module
	{
		eval $MACS2_cmd	
	} || {
		eval $MACS2_alternate_cmd
		echo -e 'Default MACS2 command failed. So, executing MACS2 with the default model parameters (--nomodel --extsize 147)' > $MACS2_outdir_default'/MACS2_Broad_exception.log'
	}
fi 	# end MACS2 file existence condition

# filter the peaks using Q or P values
# 9th field stores the Q value (in log10 scale)
# 8th field stores the P value (in log10 scale)

# now filter the peaks using Q value of 0.05 and 0.01 (using log10 scales)
QFilt1FileBroad=$MACS2_outdir_default$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt'
if [[ ! -f $QFilt1FileBroad ]]; then
	awk -v t=$Q_Thr1 '($9 > t)' $MACS2BroadPeakOutFile > $QFilt1FileBroad
fi

QFilt2FileBroad=$MACS2_outdir_default$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.01filt'
if [[ ! -f $QFilt2FileBroad || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr2 '($9 > t)' $MACS2BroadPeakOutFile > $QFilt2FileBroad
fi

#**************
# summary statistics
#**************

if [ $DEBUG_TXT == 1 ]; then
	# get the FRiP measure from this MACS2 output
	FRiP_outfile=$MACS2_outdir_default'out_FRiP.txt'
	if [[ ! -f $FRiP_outfile || $Overwrite == 1 ]]; then
		# number of reads within MACS2 narrow peaks (q threshold = 0.05)
		macs2_nreads_narrowpeak=`samtools view -cL $QFilt1File $bowtie2_BAM_prefix'.rmdup.bam'`
		FRiP_narrowpeak=`bc <<< "scale=3; ($macs2_nreads_narrowpeak * 1.0) / $uniq_mapped_read"`

		# number of reads within MACS2 broad peaks (q threshold = 0.05)
		macs2_nreads_broadpeak=`samtools view -cL $QFilt1FileBroad $bowtie2_BAM_prefix'.rmdup.bam'`
		FRiP_broadpeak=`bc <<< "scale=3; ($macs2_nreads_broadpeak * 1.0) / $uniq_mapped_read"`

		echo -e 'UniqMappedRead\tMappedReadNarrowpeak\tFRiPNarrowPeak\tMappedReadBroadpeak\tFRiPBroadPeak' > $FRiP_outfile
		echo -e '\n'$uniq_mapped_read'\t'$macs2_nreads_narrowpeak'\t'$FRiP_narrowpeak'\t'$macs2_nreads_broadpeak'\t'$FRiP_broadpeak >> $FRiP_outfile
	fi
fi


# print the summary statistics of the aligned map file
OutPeakStatFile=$MACS2_outdir_default'/Peak_Statistics.txt'
if [[ ! -f $OutPeakStatFile || $Overwrite == 1 ]]; then	
	npeakNarrowPeak=`cat $MACS2PeakOutFile | wc -l`
	npeakQ1FiltNarrowPeak=`cat $QFilt1File | wc -l`
	npeakQ2FiltNarrowPeak=`cat $QFilt2File | wc -l`

	npeakBroadPeak=`cat $MACS2PeakOutFile | wc -l`
	npeakQ1FiltBroadPeak=`cat $QFilt1FileBroad | wc -l`
	npeakQ2FiltBroadPeak=`cat $QFilt2FileBroad | wc -l`

	echo -e 'TotNarrowPeak\tNarrowPeak_Q_0.05\tNarrowPeak_Q_0.01\tTotBroadPeak\tBroadPeak_Q_0.05\tBroadPeak_Q_0.01' > $OutPeakStatFile
	echo -e $npeakNarrowPeak'\t'$npeakQ1FiltNarrowPeak'\t'$npeakQ2FiltNarrowPeak'\t'$npeakBroadPeak'\t'$npeakQ1FiltBroadPeak'\t'$npeakQ2FiltBroadPeak >> $OutPeakStatFile
fi


#================================
# MACS2 with parameters recommended from the existing studies
#================================
MACS2_outdir_ext=$OutDir'/MACS2_Ext_Tag'
if [[ $nctrl -gt 0 ]]; then
	MACS2_outdir_ext=$MACS2_outdir_ext'_with_Control/'
else
	MACS2_outdir_ext=$MACS2_outdir_ext'_No_Control/'
fi
mkdir -p $MACS2_outdir_ext

#**************
# derive the narrowPeaks
#**************
MACS2PeakOutFile=$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak'
if [[ ! -f $MACS2PeakOutFile || $Overwrite == 1 ]]; then
	MACS2_cmd="macs2 callpeak -t "$Shifted_TagAlign_File" -f BED -g "$PEAKCALLGENOMESIZE" -n "$PREFIX".macs2 -p "$MACS2_P_Val" --nomodel --nolambda --keep-dup all --shift -100 --extsize 200 --outdir "$MACS2_outdir_ext
	if [[ $nctrl -gt 0 ]]; then
		# include the control bed file also
		MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
	fi
	# execute the command
	eval $MACS2_cmd
fi

# filter the peaks using Q or P values
# 9th field stores the Q value (in log10 scale)
# 8th field stores the P value (in log10 scale)

# now filter the peaks using Q value of 0.05 and 0.01 (using log10 scales)
QFilt1File=$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt'
if [[ ! -f $QFilt1File || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr1 '($9 > t)' $MACS2PeakOutFile > $QFilt1File
fi

QFilt2File=$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt'
if [[ ! -f $QFilt2File || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr2 '($9 > t)' $MACS2PeakOutFile > $QFilt2File
fi

#================================
# annotate the peaks 
#================================

#####################
# peaks with q-value threshold of 0.01
#####################

TempPeakFile=$MACS2_outdir_ext'/TempPeak_Q0.01.bed'

# output directory storing the peak annotations
# according to Homer specifications
PeakAnnotateDir=$MACS2_outdir_ext'/Peak_Annotate_Q0.01'
mkdir -p $PeakAnnotateDir

# file containing the summarized annotations of the peak files
# according to the HOMER specifications
OutTextFile=$PeakAnnotateDir'/Out_Summary.log'

cat ${QFilt2File} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}' - > ${TempPeakFile}

# apply HOMER specific peak annotation
${HOMERPeakAnnotExec} ${TempPeakFile} ${RefChrFastaFile} -gtf ${RefChrAnnotFile} 2>${OutTextFile} > ${PeakAnnotateDir}'/Annotated_Peak_Q0.01filt.txt'

# also obtain the distance of peaks 
# from the nearest TSS sites
# using the same HOMER annotation
# check : https://www.biostars.org/p/205576/
${HOMERPeakAnnotExec} ${QFilt2File} ${RefChrFastaFile} -size 6000 -hist 50 -bedGraph $BigWig_outdir'/Inp.Sorted.bedGraph' 2>${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist_Summary.log' > ${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist.txt'

# remove the temporary peak file
rm ${TempPeakFile}

# now call an R script
# which takes 
# 1) the summary file 
# and finds the percentage of different annotations for this peak file
# and also 2) the TSS distance from the peaks
# and plots the distribution
$RPackageExec ../Analysis/PeakAnnotateHomerSummary.r $OutTextFile ${PeakAnnotateDir}'/Peak_Q0.01filt_TSSDist.txt'

#####################
# peaks with q-value threshold of 0.05
#####################

TempPeakFile=$MACS2_outdir_ext'/TempPeak_Q0.05.bed'

# output directory storing the peak annotations
# according to Homer specifications
PeakAnnotateDir=$MACS2_outdir_ext'/Peak_Annotate_Q0.05'
mkdir -p $PeakAnnotateDir

# file containing the summarized annotations of the peak files
# according to the HOMER specifications
OutTextFile=$PeakAnnotateDir'/Out_Summary.log'

cat ${QFilt1File} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}' - > ${TempPeakFile}

# apply HOMER specific peak annotation
${HOMERPeakAnnotExec} ${TempPeakFile} ${RefChrFastaFile} -gtf ${RefChrAnnotFile} 2>${OutTextFile} > ${PeakAnnotateDir}'/Annotated_Peak_Q0.05filt.txt'

# also obtain the distance of peaks 
# from the nearest TSS sites
# using the same HOMER annotation
# check : https://www.biostars.org/p/205576/
${HOMERPeakAnnotExec} ${QFilt1File} ${RefChrFastaFile} -size 6000 -hist 50 -bedGraph $BigWig_outdir'/Inp.Sorted.bedGraph' 2>${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist_Summary.log' > ${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist.txt'

# remove the temporary peak file
rm ${TempPeakFile}

# now call an R script
# which takes 
# 1) the summary file 
# and finds the percentage of different annotations for this peak file
# and also 2) the TSS distance from the peaks
# and plots the distribution
$RPackageExec ../Analysis/PeakAnnotateHomerSummary.r $OutTextFile ${PeakAnnotateDir}'/Peak_Q0.05filt_TSSDist.txt'

#**************
# derive the Broad Peaks
#**************
MACS2BroadPeakOutFile=$MACS2_outdir_ext$PREFIX'.macs2_Broad_peaks.broadPeak'
if [[ ! -f $MACS2BroadPeakOutFile || $Overwrite == 1 ]]; then
	MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -p '$MACS2_P_Val' --nomodel --nolambda --keep-dup all --broad --shift -100 --extsize 200 --outdir '$MACS2_outdir_ext
	
	if [[ $nctrl -gt 0 ]]; then
		# include the control bed file also
		MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
	fi
	# execute the command
	eval $MACS2_cmd	
fi 	# end MACS2 file existence condition

# filter the peaks using Q or P values
# 9th field stores the Q value (in log10 scale)
# 8th field stores the P value (in log10 scale)

# now filter the peaks using Q value of 0.05 and 0.01 (using log10 scales)
QFilt1FileBroad=$MACS2_outdir_ext$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt'
if [[ ! -f $QFilt1FileBroad || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr1 '($9 > t)' $MACS2BroadPeakOutFile > $QFilt1FileBroad
fi

QFilt2FileBroad=$MACS2_outdir_ext$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.01filt'
if [[ ! -f $QFilt2FileBroad || $Overwrite == 1 ]]; then
	awk -v t=$Q_Thr2 '($9 > t)' $MACS2BroadPeakOutFile > $QFilt2FileBroad
fi

#**************
# summary statistics
#**************

if [ $DEBUG_TXT == 1 ]; then
	# get the FRiP measure from this MACS2 output
	FRiP_outfile=$MACS2_outdir_ext'out_FRiP.txt'
	if [[ ! -f $FRiP_outfile || $Overwrite == 1 ]]; then
		# number of reads within MACS2 narrow peaks (q threshold = 0.05)
		macs2_nreads_narrowpeak=`samtools view -cL $QFilt1File $bowtie2_BAM_prefix'.rmdup.bam'`
		FRiP_narrowpeak=`bc <<< "scale=3; ($macs2_nreads_narrowpeak * 1.0) / $uniq_mapped_read"`

		# number of reads within MACS2 broad peaks (q threshold = 0.05)
		macs2_nreads_broadpeak=`samtools view -cL $QFilt1FileBroad $bowtie2_BAM_prefix'.rmdup.bam'`
		FRiP_broadpeak=`bc <<< "scale=3; ($macs2_nreads_broadpeak * 1.0) / $uniq_mapped_read"`

		echo -e 'UniqMappedRead\tMappedReadNarrowpeak\tFRiPNarrowPeak\tMappedReadBroadpeak\tFRiPBroadPeak' > $FRiP_outfile
		echo -e '\n'$uniq_mapped_read'\t'$macs2_nreads_narrowpeak'\t'$FRiP_narrowpeak'\t'$macs2_nreads_broadpeak'\t'$FRiP_broadpeak >> $FRiP_outfile
	fi
fi

# print the summary statistics of the aligned map file
OutPeakStatFile=$MACS2_outdir_ext'/Peak_Statistics.txt'
if [[ ! -f $OutPeakStatFile || $Overwrite == 1 ]]; then	
	npeakNarrowPeak=`cat $MACS2PeakOutFile | wc -l`
	npeakQ1FiltNarrowPeak=`cat $QFilt1File | wc -l`
	npeakQ2FiltNarrowPeak=`cat $QFilt2File | wc -l`

	npeakBroadPeak=`cat $MACS2PeakOutFile | wc -l`
	npeakQ1FiltBroadPeak=`cat $QFilt1FileBroad | wc -l`
	npeakQ2FiltBroadPeak=`cat $QFilt2FileBroad | wc -l`

	echo -e 'TotNarrowPeak\tNarrowPeak_Q_0.05\tNarrowPeak_Q_0.01\tTotBroadPeak\tBroadPeak_Q_0.05\tBroadPeak_Q_0.01' > $OutPeakStatFile
	echo -e $npeakNarrowPeak'\t'$npeakQ1FiltNarrowPeak'\t'$npeakQ2FiltNarrowPeak'\t'$npeakBroadPeak'\t'$npeakQ1FiltBroadPeak'\t'$npeakQ2FiltBroadPeak >> $OutPeakStatFile
fi

# #================================
# # here analyze the MACS2 output 
# # from the default command
# # and the extsize command
# # compute peak overlap statistics
# #================================

# PeakOverlapOutDir=$OutDir'/Peak_Overlap_Statistics'
# mkdir -p $PeakOverlapOutDir

# if [ ! -f $PeakOverlapOutDir'/Peak_Overlap.txt' ]; then
# 	/home/sourya/R-3.4.3/bin/Rscript $PeakOverlapCode --FileList $MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt':$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt' --Labels 'MACS2_Default':'MACS2_Ext' --OutPrefix $PeakOverlapOutDir'/Peak_Ov' --Dump
# fi

#================================
# MACS2 with parameters recommended from the existing studies
# **** TEST ---- here, sorted alignment file before duplicate removal is used for MACS2
#================================

# comment - sourya - for the moment
if [[ 0 == 1 ]]; then

	MACS2_outdir=$OutDir'/MACS2_NoDupRem_Align_Ext_Tag'
	if [[ $nctrl -gt 0 ]]; then
		MACS2_outdir=$MACS2_outdir'_with_Control/'
	else
		MACS2_outdir=$MACS2_outdir'_No_Control/'
	fi
	mkdir -p $MACS2_outdir

	# first form the command
	MACS2PeakOutFile=$MACS2_outdir$PREFIX'.macs2_peaks.narrowPeak'
	if [[ ! -f $MACS2PeakOutFile || $Overwrite == 1 ]]; then
		MACS2_cmd="macs2 callpeak -t "$bowtie2_BAM_prefix".bam -f BAM -g "$PEAKCALLGENOMESIZE" -n "$PREFIX".macs2 -p "$MACS2_P_Val" --nomodel --nolambda --keep-dup all --shift -100 --extsize 200 --outdir "$MACS2_outdir
		if [[ $nctrl -gt 0 ]]; then
			# include the control bed file also
			MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
		fi
		# execute the command
		$MACS2_cmd
	fi

	# filter the peaks using Q or P values
	# 9th field stores the Q value (in log10 scale)
	# 8th field stores the P value (in log10 scale)

	# now filter the peaks using Q value of 0.05 and 0.01 (using log10 scales)
	QFilt1File=$MACS2_outdir$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt'
	if [[ ! -f $QFilt1File || $Overwrite == 1 ]]; then
		awk -v t=$Q_Thr1 '($9 > t)' $MACS2PeakOutFile > $QFilt1File
	fi

	QFilt2File=$MACS2_outdir$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt'
	if [[ ! -f $QFilt2File || $Overwrite == 1 ]]; then
		awk -v t=$Q_Thr2 '($9 > t)' $MACS2PeakOutFile > $QFilt2File
	fi

	# print the summary statistics of the aligned map file
	OutPeakStatFile=$MACS2_outdir'/Peak_Statistics.txt'
	if [[ ! -f $OutPeakStatFile || $Overwrite == 1 ]]; then	
		npeak=`cat $MACS2PeakOutFile | wc -l`
		npeakQ1Filt=`cat $QFilt1File | wc -l`
		npeakQ2Filt=`cat $QFilt2File | wc -l`
		echo -e 'TotPeak\tPeak_Q_0.05\tPeak_Q_0.01' > $OutPeakStatFile
		echo -e $npeak'\t'$npeakQ1Filt'\t'$npeakQ2Filt >> $OutPeakStatFile
	fi

	if [ $DEBUG_TXT == 1 ]; then
		# get the FRiP measure from this MACS2 output
		FRiP_outfile=$MACS2_outdir'out_FRiP.txt'
		if [[ ! -f $FRiP_outfile || $Overwrite == 1 ]]; then
			macs2_nreads_peak=`samtools view -cL $QFilt1File $bowtie2_BAM_prefix'.rmdup.bam'`
			FRiP_measure=`bc <<< "scale=3; ($macs2_nreads_peak * 1.0) / $uniq_mapped_read"`
			echo -e 'Unique mapped read\tMapped read within peaks\tFRiP' > $FRiP_outfile
			echo -e '\n'$uniq_mapped_read'\t'$macs2_nreads_peak'\t'$FRiP_measure >> $FRiP_outfile
		fi
	fi

fi # end dummy if

#=========================
# processing MACS2 default command output 
# converting narrow peak and broad peak files
# to bigNarrowPeak format
# for display in the UCSC genome browser
#=========================

#****************
# convert the narrowPeak with FDR threshold = 0.05
# in the big bed format
#****************
QFilt1File=$MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt'
MACS2PeakBigBedFile=$MACS2_outdir_default$PREFIX'_bigNarrowPeak_Q0.05filt_MACS2_Default.bb'
tempPeakFile=$MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt_MACS2_Default_reduced'

# first modify the detected peaks
# to be compatible with the processing of Bigbed utility
# clipping the 5th field (score) of the MACS2 detected peaks
# within 1000
# Note: narrowPeak file has 10 fields
cat $QFilt1File | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $tempPeakFile

# now convert the converted peak file to the bigbed format
# check the reference chromosome information
# from either the input Bowtie2 alignment
# or the BigWigGenome

# Note: we have adjusted the value of RefChrSizeFile
# according to the reference genome provided
bedToBigBed -as=${BigNarrowPeakASFile} -type=bed6+4 $tempPeakFile $RefChrSizeFile $MACS2PeakBigBedFile

# delete the temporary peak file
rm $tempPeakFile

#****************
# convert the broadPeak with FDR threshold = 0.05
# in the big bed format
#****************
QFilt1FileBroad=$MACS2_outdir_default$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt'
MACS2BroadPeakBigBedFile=$MACS2_outdir_default$PREFIX'_bigBroadPeak_Q0.05filt_MACS2_Default.bb'
tempBroadPeakFile=$MACS2_outdir_default$PREFIX'.macs2_Broad_peaks.narrowPeak_MACS2_Default_Q0.05filt_reduced'

# first modify the detected peaks
# to be compatible with the processing of Bigbed utility
# clipping the 5th field (score) of the MACS2 detected peaks
# within 1000
# Note: broadPeak file has 9 fields
cat $QFilt1FileBroad | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9}' - > $tempBroadPeakFile

# now convert the converted peak file to the bigbed format
# check the reference chromosome information
# from either the input Bowtie2 alignment
# or the BigWigGenome

# Note: we have adjusted the value of RefChrSizeFile
# according to the reference genome provided
bedToBigBed -as=${BroadPeakASFile} -type=bed6+4 $tempBroadPeakFile $RefChrSizeFile $MACS2BroadPeakBigBedFile

# delete the temporary peak file
rm $tempBroadPeakFile



#=========================
# processing MACS2 extsize command output 
# converting narrow peak and broad peak files
# to bigNarrowPeak format
# for display in the UCSC genome browser
#=========================

#****************
# convert the narrowPeak with FDR threshold = 0.05
# in the big bed format
#****************
QFilt1File=$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt'
MACS2PeakBigBedFile=$MACS2_outdir_ext$PREFIX'_bigNarrowPeak_Q0.05filt_MACS2_Ext.bb'
tempPeakFile=$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt_MACS2_Ext_reduced'

# first modify the detected peaks
# to be compatible with the processing of Bigbed utility
# clipping the 5th field (score) of the MACS2 detected peaks
# within 1000
# Note: narrowPeak file has 10 fields
cat $QFilt1File | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $tempPeakFile

# now convert the converted peak file to the bigbed format
# check the reference chromosome information
# from either the input Bowtie2 alignment
# or the BigWigGenome

# Note: we have adjusted the value of RefChrSizeFile
# according to the reference genome provided
bedToBigBed -as=${BigNarrowPeakASFile} -type=bed6+4 $tempPeakFile $RefChrSizeFile $MACS2PeakBigBedFile

# delete the temporary peak file
rm $tempPeakFile

#****************
# convert the broadPeak with FDR threshold = 0.05
# in the big bed format
#****************
QFilt1FileBroad=$MACS2_outdir_ext$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt'
MACS2BroadPeakBigBedFile=$MACS2_outdir_ext$PREFIX'_bigBroadPeak_Q0.05filt_MACS2_Ext.bb'
tempBroadPeakFile=$MACS2_outdir_ext$PREFIX'.macs2_Broad_peaks.narrowPeak_MACS2_Ext_Q0.05filt_reduced'

# first modify the detected peaks
# to be compatible with the processing of Bigbed utility
# clipping the 5th field (score) of the MACS2 detected peaks
# within 1000
# Note: broadPeak file has 10 fields
cat $QFilt1FileBroad | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9}' - > $tempBroadPeakFile

# now convert the converted peak file to the bigbed format
# check the reference chromosome information
# from either the input Bowtie2 alignment
# or the BigWigGenome

# Note: we have adjusted the value of RefChrSizeFile
# according to the reference genome provided
bedToBigBed -as=${BroadPeakASFile} -type=bed6+4 $tempBroadPeakFile $RefChrSizeFile $MACS2BroadPeakBigBedFile

# delete the temporary peak file
rm $tempBroadPeakFile


# #=============================
# # now we check for individual peaks (detected from MACS2)
# # the number of reads mapped into this peak
# # also we measure the peak length
# #=============================
# python ../src/peak_distribution.py -I $MACS2_outdir2 $PREFIX'.macs2_peaks.broadPeak' -R $bowtie2_BAM_prefix'.rmdup.MAPQ'$MAPQ_THR'.bam'
# python ../src/peak_distribution.py -I $MACS2_outdir1 $PREFIX'.macs2_peaks.narrowPeak' -R $bowtie2_BAM_prefix'.rmdup.MAPQ'$MAPQ_THR'.bam'

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
