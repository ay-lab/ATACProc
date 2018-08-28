#!/bin/bash

#=================================
# this script encapsulates IDR analysis between two replicates into a single script
# First, the input BAM files are analyzed to check their read similarity
# and the BAM file with higher read is subsampled
# then these modified BAM files are used for peak calling and subsequent IDR analysis
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: ./IDR_SubSampleBAM.sh [-h] [-A BamFile1] [-B BamFile2] [-d OutDir] [-n 'IDR_test'] [-P PathIDRCode] [-c 25000] [-C control.bam]
Example:

Options:    

  -- required:
	-A  BamFile1       	 First BAM file
	-B  BamFile2       	 Second BAM file
	-d  OutDir 		 	 Output directory (absolute path preferred) which will store the IDR results.
	-P 	PathIDRCode		 Path of the IDRCode package (Kundaje et. al. after its installation)
	-n 	PREFIX 			 Prefix of output files. Default 'IDR_ATAC'.
	-c  CountPeak		 No of peaks in both replicates that will be compared. Default 25000.
	-C  CONTROLBAM		 Control file (in eiher .BAM or tagalign file in .gz format)
EOF
}

# default values of peaks that need to be retained
CountPeak=25000

# executable containing the tag align shift code
TagAlignExec='../bin/TagAlign.sh'

# default control bam file
CONTROLBAM=""

# IDR analysis code using a pair of peak files
IDR_code='./IDRAnalysis.sh'

# default prefix string
PREFIX='IDR_ATAC'

# executable of sambamba 
# for subsampling of the bam file, samtools has a bug
# so using this package
sambamba_exec=`which sambamba`

# Sourya - Note the processing of input file argument since it can be more than one file
# Note the change of notations
while getopts "A:B:n:d:c:C:P:" opt;
do
	case "$opt" in
		A) BamFile1=$OPTARG;;
		B) BamFile2=$OPTARG;;
		n) PREFIX=$OPTARG;;
		d) OutDir=$OPTARG;;
		c) CountPeak=$OPTARG;;
		C) CONTROLBAM=$OPTARG;;
		P) IDRCodeDir=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $BamFile1 ]]; then
	echo 'User did not provide the first BAM file - exit for the moment !!'
	exit 1
fi

if [[ -z $BamFile2 ]]; then
	echo 'User did not provide the second BAM file - exit for the moment !!'
	exit 1
fi

if [[ -z $IDRCodeDir ]]; then
	echo 'User did not provide the path of IDRCode package (Kundaje et. al.) - exit for the moment !!'
	exit 1
fi

if [[ -z $OutDir ]]; then
	echo 'User did not provide output directory for storing the results - exit for the moment !!'
	exit 1
fi

# create the output directory
mkdir -p $OutDir

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

# count of READS for two BAM files
readcount1=`samtools view $BamFile1 | wc -l`
readcount2=`samtools view $BamFile2 | wc -l`

TagAlignFile1=$OutDir'/temp_1_tagalign.gz'
TagAlignFile2=$OutDir'/temp_2_tagalign.gz'

if [[ $readcount1 -gt $readcount2 ]]; then
	# the first BAM file needs to be subsampled, followed by their conversion to TAG Align format
	if [ ! -f $OutDir'/temp_1.bam' ]; then
		# fraction of subsampling
		# Note: we did not use the expr operator - simple expr does not work for float
		f=$(echo "scale=2;$readcount2/$readcount1" | bc)
		# we use sambamba for the subsampling
		# and use 8 threads for faster operation
		$sambamba_exec view -h -t 8 -s $f -f bam $BamFile1 -o $OutDir'/temp_1.bam'
	fi
	# conversion of the TAG Align format
	if [ ! -f $TagAlignFile1 ]; then
		$TagAlignExec -I $OutDir'/temp_1.bam' -N 0 -O $TagAlignFile1
	fi
	if [ ! -f $TagAlignFile2 ]; then
		$TagAlignExec -I $BamFile2 -N 0 -O $TagAlignFile2
	fi
else
	# the second BAM file needs to be subsampled, followed by their conversion to TAG Align format
	if [ ! -f $OutDir'/temp_2.bam' ]; then
		# fraction of subsampling
		# Note: we did not use the expr operator - simple expr does not work for float
		f=$(echo "scale=2;$readcount1/$readcount2" | bc)
		# we use sambamba for the subsampling
		# and use 8 threads for faster operation
		$sambamba_exec view -h -t 8 -s $f -f bam $BamFile2 -o $OutDir'/temp_2.bam'
	fi
	if [ ! -f $TagAlignFile2 ]; then
		$TagAlignExec -I $OutDir'/temp_2.bam' -N 0 -O $TagAlignFile2
	fi
	if [ ! -f $TagAlignFile1 ]; then
		$TagAlignExec -I $BamFile1 -N 0 -O $TagAlignFile1
	fi
fi

#==============================================
# calling the MACS2 using the generated tag align file
#==============================================

# first we have to fix the output folders containing the MACS2 output for both the samples
# the output folder name is like MACS2_0/1_C 
# (where 0/1 indicates first or second sample)
# _C is optional and included only when control bam file is provided as input

MACS2_outdir1=$OutDir'/MACS2_0'
MACS2_outdir2=$OutDir'/MACS2_1'
if [[ ! -z $CONTROLBAM ]]; then
	MACS2_outdir1=$MACS2_outdir1'_C'
	MACS2_outdir2=$MACS2_outdir2'_C'
fi
MACS2_outdir1=$MACS2_outdir1'/'
MACS2_outdir2=$MACS2_outdir2'/'
mkdir -p $MACS2_outdir1
mkdir -p $MACS2_outdir2


# first file - MACS2
MACS2PeakOutFile1=$MACS2_outdir1$PREFIX'.macs2_peaks.narrowPeak'
if [ ! -f $MACS2PeakOutFile1 ]; then
	MACS2_cmd='macs2 callpeak -t '$TagAlignFile1' -f BED -n '$PREFIX'.macs2 --nomodel --nolambda --shift -100 --extsize 200 --outdir '$MACS2_outdir1
	if [[ ! -z $CONTROLBAM ]]; then
		# include the control file also
		MACS2_cmd=$MACS2_cmd' -c '$CONTROLBAM
	fi
	# execute the command
	$MACS2_cmd
fi

# second file - MACS2
MACS2PeakOutFile2=$MACS2_outdir2$PREFIX'.macs2_peaks.narrowPeak'
if [ ! -f $MACS2PeakOutFile2 ]; then
	MACS2_cmd='macs2 callpeak -t '$TagAlignFile2' -f BED -n '$PREFIX'.macs2 --nomodel --nolambda --shift -100 --extsize 200 --outdir '$MACS2_outdir2
	if [[ ! -z $CONTROLBAM ]]; then
		# include the control file also
		MACS2_cmd=$MACS2_cmd' -c '$CONTROLBAM
	fi
	# execute the command
	$MACS2_cmd
fi

#====================================
# now call the IDR analysis using the generated peak files
#====================================
# we have to fix the output directory where the results of IDR will be stored
# depending on the presence of control parameters
# the folders will vary
# the folders have the following format: C(0/1) depending on the input options

IDR_OutDir=$OutDir'/'
if [[ ! -z $CONTROLBAM ]]; then
	IDR_OutDir=$IDR_OutDir'C1'
else
	IDR_OutDir=$IDR_OutDir'C0'
fi
IDR_OutDir=$IDR_OutDir'_Peak'$CountPeak'/'
mkdir -p $IDR_OutDir

$IDR_code -a $MACS2PeakOutFile1 -b $MACS2PeakOutFile2 -P $IDRCodeDir -d $IDR_OutDir -n $PREFIX -c $CountPeak

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
