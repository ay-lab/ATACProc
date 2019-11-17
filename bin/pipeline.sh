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

# q-value used to compute the MACS2 peak
# we use very liberal threshold at first
# Note: previously we were using p-value as threshold
# but found that q-values are not reported then
MACS2_Q_Val=0.1

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

# default file name corresponding to the blacklisted regions
BlackListFile=""

#============================
# reference packages / executables
#============================
# R package installed - executable
RPackageExec=`which Rscript`

# python executable
PythonExec=`which python`

# picard tool executable
# picard_exec=`type -p picard.jar` #'/share/apps/picard-tools/picard-tools-2.7.1/picard.jar'

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

while getopts "C:f:r:n:g:t:m:d:a:l:c:q:w:D:O:" opt;
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
			# if [ $param == "RPackageExec" ]; then
			# 	RPackageExec=$paramval
			# fi
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
			if [ $param == "BlackListFile" ]; then
				BlackListFile=$paramval
			fi	
			if [ $param == "ATAQVPath" ]; then
				ATAQVExec=$paramval
			fi
			if [ $param == "TSSFile" ]; then
				TSSFile=$paramval
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

# if [[ -z $RPackageExec ]]; then
# 	echo 'R executable is not provided - check the configuration file - quit !! '
# 	exit 1
# fi

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

if [[ -z $ATAQVExec ]]; then
	echo 'Path of ataqv executable (ParkerLab - QC of ATAC-seq data) is not provided - check the configuration file - quit !! '
	exit 1
fi

# code in HOMER package
# which annotates peaks according to different genomic segments
HOMERPeakAnnotExec=$HOMERPath'/annotatePeaks.pl'
HOMERMotifExec=$HOMERPath'/findMotifsGenome.pl'

#----------------------------------
# select the reference genome, peak calling genome
# and the effective genome size 
# http://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
if [[ $BOWTIE2_GENOME = *"hg19"* || $BigWigGenome = *"hg19"* ]]; then	
	RefGenome='hg19'
	EGS=2864785220
	PEAKCALLGENOMESIZE='hs'
elif [[ $BOWTIE2_GENOME = *"hg38"* || $BigWigGenome = *"hg38"* ]]; then
	RefGenome='hg38'
	EGS=2913022398
	PEAKCALLGENOMESIZE='hs'
elif [[ $BOWTIE2_GENOME = *"mm9"* || $BigWigGenome = *"mm9"* ]]; then
	RefGenome='mm9'
	EGS=2620345972
	PEAKCALLGENOMESIZE='mm'
elif [[ $BOWTIE2_GENOME = *"mm10"* || $BigWigGenome = *"mm10"* ]]; then
	RefGenome='mm10'
	EGS=2652783500
	PEAKCALLGENOMESIZE='mm'
else
	RefGenome='hg19'	# default
	EGS=0
	PEAKCALLGENOMESIZE='hs'
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
	java -Xmx$MAX_MEM -jar $picard_exec MarkDuplicates -INPUT $bowtie2_BAM_prefix'.bam' -OUTPUT $bowtie2_BAM_prefix'.rmdup.bam' -ASSUME_SORTED true -REMOVE_DUPLICATES true -VALIDATION_STRINGENCY LENIENT -METRICS_FILE $bowtie2_BAM_prefix'.picard_metrics.txt'
fi

if [[ ! -f $bowtie2_BAM_prefix'.rmdup.bam.bai' ]]; then
	samtools index $bowtie2_BAM_prefix'.rmdup.bam'
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
			echo -e 'TotRead \t NumMappableRead \t nread_del_random \t nread_del_mit \t UniqMappedRead \t ReadQualThr  \t rmDupRead' > $out_readcount_file
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

# #=========================
# # add - sourya
# # check if there is any blacklist file corresponding to the current reference genome is provided
# # in such a case, discard the reads which overlap with the blacklisted region
# #=========================
# # if blacklisted genome region file is provided
# if [[ ! -z $BlackListFile && -f $BlackListFile ]]; then
# 	if [ ! -f $bowtie2_BAM_prefix'_possible_overlap_blacklist.bam' ]; then
# 		# -v option reports only those entries which do not have overlap with the bed file
# 		bedtools intersect -v -abam $bowtie2_BAM_prefix'.bam' -b $BlackListFile > $bowtie2_BAM_prefix'_remove_blacklist.bam'
# 		# store the old copy, and also create samtools index
# 		mv $bowtie2_BAM_prefix'.bam' $bowtie2_BAM_prefix'_possible_overlap_blacklist.bam'
# 		samtools index $bowtie2_BAM_prefix'_possible_overlap_blacklist.bam'
# 		# rename the new copy (without blacklisted regions) as the copy to be processed further
# 		mv $bowtie2_BAM_prefix'_remove_blacklist.bam' $bowtie2_BAM_prefix'.bam'
# 	fi
# 	# get the number of reads after removing blacklisted region
# 	numReadBlackList=`samtools view $bowtie2_BAM_prefix'.bam' | cut -f 1 | sort | uniq | wc -l`
# else
# 	numReadBlackList=$nread_qual
# fi

#============================
# use deeptools command to create a de-duplicated BAM file 
# and mapping quality > specified mapping quality threshold
# and without any blacklisted genome read
# and also forward and reverse strands are shifted according to the ATAC seq shift
#============================
ShiftedBAMFile=$bowtie2_BAM_prefix'_TN5_Shift.bam'
tempfile=$bowtie2_BAM_prefix'_TN5_Shift_temp.bam'
if [[ ! -f $ShiftedBAMFile || $Overwrite == 1 ]]; then
	if [[ ! -z $BlackListFile && -f $BlackListFile ]]; then
		$DeepToolsDir'/alignmentSieve' --bam $bowtie2_BAM_prefix'.rmdup.bam' --outFile $tempfile --numberOfProcessors $THREADS --ATACshift --ignoreDuplicates --minMappingQuality $MAPQ_THR --blackListFileName $BlackListFile
	else
		$DeepToolsDir'/alignmentSieve' --bam $bowtie2_BAM_prefix'.rmdup.bam' --outFile $tempfile --numberOfProcessors $THREADS --ATACshift --ignoreDuplicates --minMappingQuality $MAPQ_THR 
	fi
	# sort the file
	samtools sort -o $ShiftedBAMFile $tempfile
	if [[ ! -f $ShiftedBAMFile'.bai' ]]; then
		samtools index $ShiftedBAMFile
	fi
	rm $tempfile
fi

#============================
# use deeptools command on the above generated shifted BAM file
# to get the nucleosome free reads BAM file (fragment length between 0 to 100)
# mononucleosome reads BAM file (fragment length between 180 to 247)
# dinucleosome reads BAM file (fragment length between 315 to 473)
# trinucleosome reads BAM file (fragment length between 558 to 615)
# these read lengths are obtained from 
# https://rdrr.io/github/jianhong/ATACseqQC/man/splitGAlignmentsByCut.html
#============================
File_NFR=$bowtie2_outdir'/NucleosomeFree.bam'
File_1N=$bowtie2_outdir'/mononucleosome.bam'
File_2N=$bowtie2_outdir'/dinucleosome.bam'
File_3N=$bowtie2_outdir'/trinucleosome.bam'
tempfile=$bowtie2_outdir'/temp_nuc.bam'
if [[ ! -f $File_NFR ]]; then
	$DeepToolsDir'/alignmentSieve' --bam $ShiftedBAMFile --outFile $tempfile --numberOfProcessors $THREADS --minFragmentLength 0 --maxFragmentLength 100
	samtools sort -o $File_NFR $tempfile
fi
if [[ ! -f $File_1N ]]; then
	$DeepToolsDir'/alignmentSieve' --bam $ShiftedBAMFile --outFile $tempfile --numberOfProcessors $THREADS --minFragmentLength 180 --maxFragmentLength 247
	samtools sort -o $File_1N $tempfile
fi
if [[ ! -f $File_2N ]]; then
	$DeepToolsDir'/alignmentSieve' --bam $ShiftedBAMFile --outFile $tempfile --numberOfProcessors $THREADS --minFragmentLength 315 --maxFragmentLength 473
	samtools sort -o $File_2N $tempfile
fi
if [[ ! -f $File_3N ]]; then
	$DeepToolsDir'/alignmentSieve' --bam $ShiftedBAMFile --outFile $tempfile --numberOfProcessors $THREADS --minFragmentLength 558 --maxFragmentLength 615
	samtools sort -o $File_3N $tempfile
fi
if [[ -f $tempfile ]]; then
	rm $tempfile
fi

# combine File_1N, File_2N and File_3N
# to get alignments of one or more nucleosomes (+1N)
File_Nucleosome_Merge=$bowtie2_outdir'/Merged_nucleosome.bam'
if [[ ! -f $File_Nucleosome_Merge ]]; then
	samtools merge $File_Nucleosome_Merge $File_NFR $File_1N $File_2N $File_3N 
	samtools index $File_Nucleosome_Merge
fi

#============================
# create a shifted tagalign file (TN5) from the shifted BAM file
# here no further shift operation is required
# as mentioned in the following links:
# https://www.biostars.org/p/187204/
# http://seqanswers.com/forums/archive/index.php/t-59219.html
# https://github.com/kundajelab/atac_dnase_pipelines
#============================
# # amount of shift done in forward and reverse strand
# # for producing the TN5 shift
# # this is with respect to the original paper 
# # This tagalign file is later used for MACS2 peak calling
# fwdshft=4
# revshft=5

# Shifted_TagAlign_File=$bowtie2_BAM_prefix'.TN5.tagAlign.gz'
Shifted_TagAlign_File=$bowtie2_BAM_prefix'_TN5_Shift.bed'

if [[ ! -f $Shifted_TagAlign_File || $Overwrite == 1 ]]; then

	# old command - sourya
	# call the utility function to shift the bam file
	# with respect to forward and reverse strands
	# Note: the input files (which can be multiple in the target function)
	# are placed at the last part of a command	
	# $TagAlignExec -N 0 -f $fwdshft -r $revshft -O $Shifted_TagAlign_File -q $MAPQ_THR -I $bowtie2_BAM_prefix'.rmdup.bam'

	# new command - sourya
	# just convert the shifted BAM file in BEDPE format - use deeptools command
	# compatible for MACS2
	$DeepToolsDir'/alignmentSieve' --bam $ShiftedBAMFile --outFile $Shifted_TagAlign_File --BED --numberOfProcessors $THREADS
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
	# other definition is number of distinct genome position for 
	# uniquely mapped reads / number of uniquely mapped reads
	# we follow this second definition (ENCODE PAPER)
	# discriminate between single end and paired end reads
	if [ $paired_read == 1 ]; then
		# one uniq_mapped_read corresponds to two distinct genome positions
		NRFval=`bc <<< "scale=3; ($uniqgenomepos * 0.5) / $uniq_mapped_read"`
	else
		NRFval=`bc <<< "scale=3; ($uniqgenomepos * 1.0) / $uniq_mapped_read"`
	fi
	
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
for MACS2COMMANDTYPE in 'Default' 'ExtSize'; do

	# MACS2COMMANDTYPE = 'Default'
	# calling MACS2 peaks with de-duplicated reads and with default commands
	# MACS2COMMANDTYPE = 'ExtSize'
	# calling MACS2 peaks with de-duplicated reads and with ExtSize commands

	echo 'Deriving MACS2 peak type : '$MACS2COMMANDTYPE
	
	# current MACS2 peak output directory name
	# depends on using the control BAM files
	if [[ $MACS2COMMANDTYPE == 'Default' ]]; then
		CURR_MACS2_OUTDIR=$OutDir'/MACS2_Default_Tag'
	elif [[ $MACS2COMMANDTYPE == 'ExtSize' ]]; then
		CURR_MACS2_OUTDIR=$OutDir'/MACS2_Ext_Tag'
	# else
	# 	CURR_MACS2_OUTDIR=$OutDir'/MACS2_NoDupRem_Align_Ext_Tag'
	fi
	if [[ $nctrl -gt 0 ]]; then
		CURR_MACS2_OUTDIR=$CURR_MACS2_OUTDIR'_with_Control/'
	else
		CURR_MACS2_OUTDIR=$CURR_MACS2_OUTDIR'_No_Control/'
	fi
	mkdir -p $CURR_MACS2_OUTDIR

	# now classify according to the peak type : narrow peak or broad peak
	for PEAKTYPE in 'narrow' 'broad'; do

		# if [[ $MACS2COMMANDTYPE == 'noDupRem' && $PEAKTYPE == 'broad' ]]; then
		# 	continue
		# fi

		# current MACS2 output peak file
		if [[ $PEAKTYPE == 'narrow' ]]; then
			MACS2PeakOutFile=$CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak'
		else
			MACS2PeakOutFile=$CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak'
		fi

		# check if the peak file already exists
		if [[ ! -f $MACS2PeakOutFile || $Overwrite == 1 ]]; then
			if [[ $MACS2COMMANDTYPE == 'Default' && $PEAKTYPE == 'narrow' ]]; then
				# main MACS2 command
				MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2 -q '$MACS2_Q_Val' --keep-dup all --call-summits --outdir '$CURR_MACS2_OUTDIR
				# this is an alernate command 
				# only invoked when the above comamnd fails
				MACS2_alternate_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2 -q '$MACS2_Q_Val' --nomodel --extsize 147 --keep-dup all --call-summits --outdir '$CURR_MACS2_OUTDIR
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
					echo -e 'Default MACS2 command failed. So, executing MACS2 with the default model parameters (--nomodel --extsize 147)' > $CURR_MACS2_OUTDIR'/MACS2_exception.log'
				}

			elif [[ $MACS2COMMANDTYPE == 'Default' && $PEAKTYPE == 'broad' ]]; then
				# main MACS2 command
				MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -q '$MACS2_Q_Val' --keep-dup all --broad --outdir '$CURR_MACS2_OUTDIR	
				# this is an alernate command 
				# only invoked when the above comamnd fails
				MACS2_alternate_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -q '$MACS2_Q_Val' --keep-dup all --broad --nomodel --extsize 147 --outdir '$CURR_MACS2_OUTDIR
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
					echo -e 'Default MACS2 command failed. So, executing MACS2 with the default model parameters (--nomodel --extsize 147)' > $CURR_MACS2_OUTDIR'/MACS2_Broad_exception.log'
				}

			elif [[ $MACS2COMMANDTYPE == 'ExtSize' && $PEAKTYPE == 'narrow' ]]; then
				# main MACS2 command
				MACS2_cmd="macs2 callpeak -t "$Shifted_TagAlign_File" -f BED -g "$PEAKCALLGENOMESIZE" -n "$PREFIX".macs2 -q "$MACS2_Q_Val" --nomodel --nolambda --shift -100 --extsize 200 --keep-dup all --call-summits --outdir "$CURR_MACS2_OUTDIR
				if [[ $nctrl -gt 0 ]]; then
					# include the control bed file also
					MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
				fi
				# execute the command
				eval $MACS2_cmd

			elif [[ $MACS2COMMANDTYPE == 'ExtSize' && $PEAKTYPE == 'broad' ]]; then
				# main MACS2 command
				MACS2_cmd='macs2 callpeak -t '$Shifted_TagAlign_File' -f BED -g '$PEAKCALLGENOMESIZE' -n '$PREFIX'.macs2_Broad -q '$MACS2_Q_Val' --nomodel --nolambda --broad --shift -100 --extsize 200 --keep-dup all --outdir '$CURR_MACS2_OUTDIR
				if [[ $nctrl -gt 0 ]]; then
					# include the control bed file also
					MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
				fi
				# execute the command
				eval $MACS2_cmd	

			# elif [[ $MACS2COMMANDTYPE == 'noDupRem' && $PEAKTYPE == 'narrow' ]]; then
			# 	# main MACS2 command
			# 	MACS2_cmd="macs2 callpeak -t "$bowtie2_BAM_prefix".bam -f BAM -g "$PEAKCALLGENOMESIZE" -n "$PREFIX".macs2 -q "$MACS2_Q_Val" --nomodel --nolambda --shift -100 --extsize 200 --keep-dup all --call-summits --outdir "$CURR_MACS2_OUTDIR
			# 	# if [[ $nctrl -gt 0 ]]; then
			# 	# 	# include the control bed file also
			# 	# 	MACS2_cmd=$MACS2_cmd' -c '$Control_TagAlign_File_IDR
			# 	# fi
			# 	# execute the command
			# 	$MACS2_cmd

			fi

		fi 	# end check MACS2 peak file existence

		# filter the peaks using Q value of 0.05 and 0.01 (using log10 scales) - 9th field
		for FDRpct in 5 1; do
			#============
			# q-value based filtering
			#============
			QFiltFile=$MACS2PeakOutFile'_Q0.0'$FDRpct'filt'
			if [ $FDRpct == 5 ]; then
				QValThr=$Q_Thr1
			else
				QValThr=$Q_Thr2
			fi
			# sourya - currently commented
			if [[ ! -f $QFiltFile || $Overwrite == 1 ]]; then
				awk -v t=$QValThr '($9 > t)' $MACS2PeakOutFile > $QFiltFile
			fi

			#============
			# also annotate the peaks using HOMER
			#============
			TempPeakFile=$CURR_MACS2_OUTDIR'/TempPeak_Q0.0'$FDRpct'.bed'
			# output directory storing the peak annotations
			PeakAnnotateDir=$CURR_MACS2_OUTDIR'/Peak_Annotate_Q0.0'$FDRpct
			mkdir -p $PeakAnnotateDir

			cat ${QFiltFile} | awk '{print $1"\t"$2"\t"$3"\t"$4"\t0"}' - > ${TempPeakFile}

			# file containing the summarized annotations of the peak files
			OutTextFile=$PeakAnnotateDir'/Out_Summary.log'
			
			if [[ ! -f $OutTextFile || ! -f ${PeakAnnotateDir}'/Annotated_Peak_Q0.0'$FDRpct'filt.txt' || ! -f ${PeakAnnotateDir}'/Peak_Q0.0'$FDRpct'filt_TSSDist.txt' ]]; then

				# apply HOMER specific peak annotation
				${HOMERPeakAnnotExec} ${TempPeakFile} ${RefChrFastaFile} -gtf ${RefChrAnnotFile} 2>${OutTextFile} > ${PeakAnnotateDir}'/Annotated_Peak_Q0.0'$FDRpct'filt.txt'
			
				# also obtain the distance of peaks from the nearest TSS sites 
				# using the HOMER annotation
				# check : https://www.biostars.org/p/205576/
				${HOMERPeakAnnotExec} ${QFiltFile} ${RefChrFastaFile} -size 6000 -hist 50 -bedGraph $BigWig_outdir'/Inp.Sorted.bedGraph' 2>${PeakAnnotateDir}'/Peak_Q0.0'$FDRpct'filt_TSSDist_Summary.log' > ${PeakAnnotateDir}'/Peak_Q0.0'$FDRpct'filt_TSSDist.txt'

				# remove the temporary peak file
				rm ${TempPeakFile}
			
				# now call an R script which takes the summary file 
				# and finds the percentage of different annotations for this peak file
				# and also 2) the TSS distance from the peaks, and plots the distribution
				$RPackageExec ../Analysis/PeakAnnotateHomerSummary.r $OutTextFile ${PeakAnnotateDir}'/Peak_Q0.0'$FDRpct'filt_TSSDist.txt'

			fi

			#===============
			# create big bed file of the peaks for visualization in UCSC genome browser
			#===============
			if [[ $MACS2COMMANDTYPE == 'Default' || $MACS2COMMANDTYPE == 'ExtSize' ]]; then
				if [[ $PEAKTYPE == 'narrow' ]]; then
					
					MACS2PeakBigBedFile=$CURR_MACS2_OUTDIR$PREFIX'_bigNarrowPeak_Q0.0'$FDRpct'filt_MACS2_'$MACS2COMMANDTYPE'.bb'

					if [[ ! -f $MACS2PeakBigBedFile ]]; then
					
						tempPeakFile=$CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak_Q0.0'$FDRpct'filt_MACS2_'$MACS2COMMANDTYPE'_reduced'
					
						# first modify the detected peaks to be compatible with the processing of Bigbed utility
						# clipping the 5th field (score) of the MACS2 detected peaks within 1000
						# Note: narrowPeak file has 10 fields
						cat $QFiltFile | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $tempPeakFile
					
						# now convert the converted peak file to the bigbed format
						# check the reference chromosome information
						# from either the input Bowtie2 alignment or the BigWigGenome
						# Note: we have adjusted the value of RefChrSizeFile
						# according to the reference genome provided
						bedToBigBed -as=${BigNarrowPeakASFile} -type=bed6+4 $tempPeakFile $RefChrSizeFile $MACS2PeakBigBedFile
						
						# delete the temporary peak file
						rm $tempPeakFile

					fi

				else					
					
					MACS2BroadPeakBigBedFile=$CURR_MACS2_OUTDIR$PREFIX'_bigBroadPeak_Q0.0'$FDRpct'filt_MACS2_'$MACS2COMMANDTYPE'.bb'
					
					if [[ ! -f $MACS2BroadPeakBigBedFile ]]; then

						tempBroadPeakFile=$CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_MACS2_'$MACS2COMMANDTYPE'_Q0.0'$FDRpct'filt_reduced'
						
						# first modify the detected peaks to be compatible with the processing of Bigbed utility
						# clipping the 5th field (score) of the MACS2 detected peaks within 1000
						cat $QFiltFile | awk '{if ($1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|Y)$/) {print $0}}' - | awk 'function min(p,q) {return p < q ? p : q} {print $1"\t"$2"\t"$3"\t"$4"\t"min($5,1000)"\t"$6"\t"$7"\t"$8"\t"$9}' - > $tempBroadPeakFile
						# now convert the converted peak file to the bigbed format
						# check the reference chromosome information
						# from either the input Bowtie2 alignment or the BigWigGenome
						# Note: we have adjusted the value of RefChrSizeFile
						# according to the reference genome provided
						bedToBigBed -as=${BroadPeakASFile} -type=bed6+4 $tempBroadPeakFile $RefChrSizeFile $MACS2BroadPeakBigBedFile
						# delete the temporary peak file
						rm $tempBroadPeakFile

					fi					
				fi
			fi

		done 	# end FDRpct loop

	done 	# end PEAKTYPE loop

	#================
	# summary statistics
	#================
	if [ $DEBUG_TXT == 1 ]; then
		# get the FRiP measure from this MACS2 output
		FRiP_outfile=$CURR_MACS2_OUTDIR'out_FRiP.txt'

		if [[ ! -f $FRiP_outfile || $Overwrite == 1 ]]; then
			# number of reads within MACS2 narrow peaks (q threshold = 0.05)
			
			# old code - comment - sourya
			# -c (count option) returns incorrect numbers for paired end reads
			# macs2_nreads_narrowpeak=`samtools view -cL $CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt' $bowtie2_BAM_prefix'.rmdup.bam'`
			# new code - dump filtered reads first and then count the number of reads
			macs2_nreads_narrowpeak=`samtools view -L $CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt' $bowtie2_BAM_prefix'.rmdup.bam' | cut -f 1 | sort | uniq | wc -l`
			FRiP_narrowpeak=`bc <<< "scale=3; ($macs2_nreads_narrowpeak * 1.0) / $uniq_mapped_read"`

			# number of reads within MACS2 broad peaks (q threshold = 0.05)
			# old code - comment - sourya
			# -c (count option) returns incorrect numbers for paired end reads
			# macs2_nreads_broadpeak=`samtools view -cL $CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt' $bowtie2_BAM_prefix'.rmdup.bam'`
			# new code - dump filtered reads first and then count the number of reads
			macs2_nreads_broadpeak=`samtools view -L $CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt' $bowtie2_BAM_prefix'.rmdup.bam' | cut -f 1 | sort | uniq | wc -l`		
			FRiP_broadpeak=`bc <<< "scale=3; ($macs2_nreads_broadpeak * 1.0) / $uniq_mapped_read"`

			echo -e 'UniqMappedRead\tMappedReadNarrowpeak\tFRiPNarrowPeak\tMappedReadBroadpeak\tFRiPBroadPeak' > $FRiP_outfile
			echo -e '\n'$uniq_mapped_read'\t'$macs2_nreads_narrowpeak'\t'$FRiP_narrowpeak'\t'$macs2_nreads_broadpeak'\t'$FRiP_broadpeak >> $FRiP_outfile
		fi
	fi 	# end DEBUG_TXT condition

	# print the summary statistics of the aligned map file
	OutPeakStatFile=$CURR_MACS2_OUTDIR'/Peak_Statistics.txt'
	if [[ ! -f $OutPeakStatFile || $Overwrite == 1 ]]; then	
		# summary for narrow peaks
		npeakNarrowPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak' | wc -l`
		npeakQ1FiltNarrowPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak_Q0.05filt' | wc -l`
		npeakQ2FiltNarrowPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt' | wc -l`

		# summary for broad peaks
		npeakBroadPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak' | wc -l`
		npeakQ1FiltBroadPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.05filt' | wc -l`
		npeakQ2FiltBroadPeak=`cat $CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.01filt' | wc -l`

		echo -e 'TotNarrowPeak\tNarrowPeak_Q_0.05\tNarrowPeak_Q_0.01\tTotBroadPeak\tBroadPeak_Q_0.05\tBroadPeak_Q_0.01' > $OutPeakStatFile
		echo -e $npeakNarrowPeak'\t'$npeakQ1FiltNarrowPeak'\t'$npeakQ2FiltNarrowPeak'\t'$npeakBroadPeak'\t'$npeakQ1FiltBroadPeak'\t'$npeakQ2FiltBroadPeak >> $OutPeakStatFile
	fi

done 	# end MACS2COMMANDTYPE loop


#================================
# run ataqv QC pipeline using the aligned BAM file (with blacklisted regions)
# reference genome TSS file
# blacklist genome file (if provided)
# and the broadpeak file 
#================================
ATACQVOutDir=$OutDir'/QC_ataqv_ParkerLab'
mkdir -p $ATACQVOutDir

# we consider broad peaks generated from 'ExtSize' peak option
CURR_MACS2_OUTDIR=$OutDir'/MACS2_Ext_Tag'

if [[ ! -z $BlackListFile && -f $BlackListFile ]]; then	
	# prefix string of output json file name which will be the output of ataqv pipeline
	outJsonFilenamePrefix=$ATACQVOutDir'/'$PREFIX'_ataqv_with_blacklist_MACS2_Ext_Tag'
else
	# prefix string of output json file name which will be the output of ataqv pipeline
	outJsonFilenamePrefix=$ATACQVOutDir'/'$PREFIX'_ataqv_without_blacklist_MACS2_Ext_Tag'
fi
if [[ $nctrl -gt 0 ]]; then
	outJsonFilenamePrefix=$outJsonFilenamePrefix'_with_Control'
	CURR_MACS2_OUTDIR=$CURR_MACS2_OUTDIR'_with_Control/'
else
	outJsonFilenamePrefix=$outJsonFilenamePrefix'_No_Control'
	CURR_MACS2_OUTDIR=$CURR_MACS2_OUTDIR'_No_Control/'
fi

for FDRpct in 5 1; do
	inpBroadPeakFile=$CURR_MACS2_OUTDIR$PREFIX'.macs2_Broad_peaks.broadPeak_Q0.0'$FDRpct'filt'
	outJsonFilename=$outJsonFilenamePrefix'_Q0.0'$FDRpct'filt.ataqv.json.gz'	
	if [[ ! -z $BlackListFile && -f $BlackListFile ]]; then	
		# call ataqv pipeline with blacklisted genome region
		$ATAQVExec --peak-file $inpBroadPeakFile --name $PREFIX --metrics-file $outJsonFilename --excluded-region-file $BlackListFile --tss-file ${RefChrAnnotFile} --ignore-read-groups human $bowtie2_BAM_prefix'.rmdup.bam' > $ATACQVOutDir'/'$PREFIX'.ataqv.out'
	else
		# call ataqv pipeline
		$ATAQVExec --peak-file $inpBroadPeakFile --name $PREFIX --metrics-file $outJsonFilename --tss-file ${RefChrAnnotFile} --ignore-read-groups human $bowtie2_BAM_prefix'.rmdup.bam' > $ATACQVOutDir'/'$PREFIX'.ataqv.out'
	fi
done


#================================
# analyzing TF footprints
#================================
# for FDRpct in 5 1; do
for FDRpct in 5; do
	if [[ $nctrl -gt 0 ]]; then	
		inpPeakFile=$OutDir'/MACS2_Ext_Tag_with_Control/'$PREFIX'.macs2_peaks.narrowPeak_Q0.0'$FDRpct'filt'
		MotifOutDir=$OutDir'/Motif_MACS2_Ext_Tag_with_Control_narrowPeak_Q0.0'$FDRpct'filt'
	else
		inpPeakFile=$OutDir'/MACS2_Ext_Tag_No_Control/'$PREFIX'.macs2_peaks.narrowPeak_Q0.0'$FDRpct'filt'
		MotifOutDir=$OutDir'/Motif_MACS2_Ext_Tag_No_Control_narrowPeak_Q0.0'$FDRpct'filt'
	fi

	#for summitoffsetval in 500 200; do
	for summitoffsetval in 200; do
		
		motifPeakFile=$MotifOutDir'/Motif_Complete_Peaks_SummitOffset_'$summitoffsetval'/Peaks_Summit_Offset_'$summitoffsetval'bp.bed'

		# analyzing complete set of peaks, and extracting summits
		$RPackageExec ../Imp_Scripts/Motif_HOMER.R --MotifFindExec $HOMERMotifExec --RefGenome $RefGenome --PeakFile $inpPeakFile --OutDir $MotifOutDir --SummitOffset $summitoffsetval	
		
		# call HINT-ATAC specific motif footprinting
		$RPackageExec ../Imp_Scripts/Footprint_HINT_ATAC.R --MotifPeak $motifPeakFile --AllRead $ShiftedBAMFile --NFRRead $File_NFR --NFRANDNuclRead $File_Nucleosome_Merge --OutDir $MotifOutDir'/Motif_Complete_Peaks_SummitOffset_'$summitoffsetval'/Footprint_HINT_ATAC' --RefGenome $RefGenome --PE $paired_read				
		
		# motif finding by filtering peaks such that -log10(p-value) > 50
		pvalThr=50
		
		motifPeakFile=$MotifOutDir'/Motif_Peaks_PvalThr_'$pvalThr'_SummitOffset_'$summitoffsetval'/Filtered_Peaks_PvalThr_Summit_Offset_'$summitoffsetval'bp.bed'

		$RPackageExec ../Imp_Scripts/Motif_HOMER.R --MotifFindExec $HOMERMotifExec --RefGenome $RefGenome --PeakFile $inpPeakFile --OutDir $MotifOutDir --PValThr $pvalThr --SummitOffset $summitoffsetval

		# call HINT-ATAC specific motif footprinting
		$RPackageExec ../Imp_Scripts/Footprint_HINT_ATAC.R --MotifPeak $motifPeakFile --AllRead $ShiftedBAMFile --NFRRead $File_NFR --NFRANDNuclRead $File_Nucleosome_Merge --OutDir $MotifOutDir'/Motif_Peaks_PvalThr_'$pvalThr'_SummitOffset_'$summitoffsetval'/Footprint_HINT_ATAC' --RefGenome $RefGenome --PE $paired_read

	done
done

#================================
# plot the TSS enrichment of peaks and its surrounding regions
# using deeptools utility
# currently commented - sourya
#================================
TSSEnrichmentOutDir=$OutDir'/TSS_Enrichment_Peaks'
mkdir -p $TSSEnrichmentOutDir

for offsetval in 5000 1000; do
# for offsetval in 1000; do
	for FDRpct in 5 1; do
		if [[ $nctrl -gt 0 ]]; then	
			inpPeakFile=$OutDir'/MACS2_Ext_Tag_with_Control/'$PREFIX'.macs2_peaks.narrowPeak_Q0.0'$FDRpct'filt'
			CurrOutDir=$TSSEnrichmentOutDir'/MACS2_Ext_Tag_with_Control/macs2_narrowPeak_Q0.0'$FDRpct'filt_Offset_'$offsetval
		else
			inpPeakFile=$OutDir'/MACS2_Ext_Tag_No_Control/'$PREFIX'.macs2_peaks.narrowPeak_Q0.0'$FDRpct'filt'
			CurrOutDir=$TSSEnrichmentOutDir'/MACS2_Ext_Tag_No_Control/macs2_narrowPeak_Q0.0'$FDRpct'filt_Offset_'$offsetval
		fi
		# check if TSS file is provided as an input
		# Note: it should be separate from the reference gene annotation file
		# first three columns should have the TSS information
		if [[ ! -z $TSSFile && -f $TSSFile ]]; then	
			$RPackageExec ../Imp_Scripts/Peak_Enrichment.R --BigWigFile $BigWig_outdir1'/'$PREFIX'_NormCov.bw' --Label $PREFIX --DeepToolsDir $DeepToolsDir --TSSFile $TSSFile --PeakFile $inpPeakFile --OutDir $CurrOutDir --Offset $offsetval
		else 
			$RPackageExec ../Imp_Scripts/Peak_Enrichment.R --BigWigFile $BigWig_outdir1'/'$PREFIX'_NormCov.bw' --Label $PREFIX --DeepToolsDir $DeepToolsDir --PeakFile $inpPeakFile --OutDir $CurrOutDir --Offset $offsetval
		fi
	done
done







# # #================================
# # # here analyze the MACS2 output 
# # # from the default command
# # # and the extsize command
# # # compute peak overlap statistics
# # #================================

# # PeakOverlapOutDir=$OutDir'/Peak_Overlap_Statistics'
# # mkdir -p $PeakOverlapOutDir

# # if [ ! -f $PeakOverlapOutDir'/Peak_Overlap.txt' ]; then
# # 	/home/sourya/R-3.4.3/bin/Rscript $PeakOverlapCode --FileList $MACS2_outdir_default$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt':$MACS2_outdir_ext$PREFIX'.macs2_peaks.narrowPeak_Q0.01filt' --Labels 'MACS2_Default':'MACS2_Ext' --OutPrefix $PeakOverlapOutDir'/Peak_Ov' --Dump
# # fi

# # #=============================
# # # now we check for individual peaks (detected from MACS2)
# # # the number of reads mapped into this peak
# # # also we measure the peak length
# # #=============================
# # python ../src/peak_distribution.py -I $MACS2_outdir2 $PREFIX'.macs2_peaks.broadPeak' -R $bowtie2_BAM_prefix'.rmdup.MAPQ'$MAPQ_THR'.bam'
# # python ../src/peak_distribution.py -I $MACS2_outdir1 $PREFIX'.macs2_peaks.narrowPeak' -R $bowtie2_BAM_prefix'.rmdup.MAPQ'$MAPQ_THR'.bam'

# #==============
# # remove temporary files
# if [[ -f $curr_tagalign_file ]]; then
# 	rm $curr_tagalign_file
# fi

# if [[ -f $temp_NRF_PBC_file ]]; then
# 	rm $temp_NRF_PBC_file
# fi

# if [[ -f $Shifted_TagAlign_File ]]; then
# 	rm $Shifted_TagAlign_File
# fi

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
