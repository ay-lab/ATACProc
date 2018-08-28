#!/bin/bash 

#=================================
# this script is used to perform IDR analysis on a given ATAC seq replicates
# it uses the idrcode package provided by Anshul Kundaje et. al.
#=================================
# developed by - Sourya Bhattacharyya
# date: 11th july 2017
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

Options:    

  -- required:
	-a  FILE1        	 First file containing peak information (in either narrowpeak format or narrowpeak.gz format)
	-b  FILE2        	 Second file containing peak information (in either narrowpeak format or narrowpeak.gz format)
	-d  OutDir 		 	 Output directory containing the IDR results
	-P 	PathIDRCode		 Path of the IDRCode package (Kundaje et. al. after its installation)
	-n 	PREFIX 			 Prefix of output file
	-c 	SampledPeakCount Number of peaks which will be sampled from the input peak files  (default 25000)
EOF
}

# # name of the folder containing IDR results
# IDR_OutFold='IDR_Overlap0_PVal'

# executable (R code) of the batch consistency analysis
# IDRCodeDir='/home/sourya/packages/idrCode/'
exec1='batch-consistency-analysis.r'

# default values
PREFIX='IDR_ATAC'

OutDir=`pwd`

# Number of peaks sampled from the original peak detection output
SampledPeakCount=25000

while getopts "a:b:d:n:c:P:" opt;
do
	case "$opt" in
		a) FILE1=$OPTARG;;
		b) FILE2=$OPTARG;;
		d) OutDir=$OPTARG;;
		n) PREFIX=$OPTARG;;
		c) SampledPeakCount=$OPTARG;;
		P) IDRCodeDir=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $FILE1 ]]; then
	echo 'User should provide two input peak files (in a bed file or in gzipped bed file) !!'
	exit 1
else
	echo 'Input peak file 1: '$FILE1
fi

if [[ -z $FILE2 ]]; then
	echo 'User should provide two input peak files (in a bed file or in gzipped bed file) !!'
	exit 1
else
	echo 'Input peak file 2: '$FILE2
fi

if [[ -z $IDRCodeDir ]]; then
	echo 'User did not provide the path of IDRCode package (Kundaje et. al.) - exit for the moment !!'
	exit 1
fi

echo $OutDir
mkdir -p $OutDir

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

# log the replicates
# if [ ! -f $OutDir'/ReplicaNames.txt' ]; then
	echo -e "IDR Analysis of the following two peak files: \n File 1: ${FILE1} \n File 2: ${FILE2}" > $OutDir'/ReplicaNames.txt'
# fi

# check the extension of both input files
filebase1=$(basename "$FILE1")
filebase2=$(basename "$FILE2")

#===========================================
# extract first 25K (or the count provided) significant peaks from the given peak files and store them
# this will be actually used for IDR analysis
# such significance is decided by the 8th field (P value)
# as instructed in the ENCODE IDR documentation

if [[ $filebase1 =~ \.gz$ ]]; then
	# first file is a gzipped file
	#ConvFile1=${FILE1%.gz}'_first_'$SampledPeakCount'.gz'
	ConvFile1=${FILE1%.gz}'_first_'$SampledPeakCount
	echo 'Subsampled peak file corresponding to the first input: '$ConvFile1
	#if [ ! -f $ConvFile1 ]; then
		zcat $FILE1 | sort -k8,8nr > $OutDir'/temp1.txt'
		#head -n $SampledPeakCount $OutDir'/temp1.txt' | gzip -c > $ConvFile1
		head -n $SampledPeakCount $OutDir'/temp1.txt' > $ConvFile1
		rm $OutDir'/temp1.txt'
	#fi
else
	#ConvFile1=${FILE1}'_first_'$SampledPeakCount'.gz'
	ConvFile1=${FILE1}'_first_'$SampledPeakCount
	echo 'Subsampled peak file corresponding to the first input: '$ConvFile1
	#if [ ! -f $ConvFile1 ]; then
		cat $FILE1 | sort -k8,8nr > $OutDir'/temp1.txt'
		#head -n $SampledPeakCount $OutDir'/temp1.txt' | gzip -c > $ConvFile1
		head -n $SampledPeakCount $OutDir'/temp1.txt' > $ConvFile1
		rm $OutDir'/temp1.txt'
	#fi
fi

if [[ $filebase2 =~ \.gz$ ]]; then
	# first file is a gzipped file
	#ConvFile2=${FILE2%.gz}'_first_'$SampledPeakCount'.gz'
	ConvFile2=${FILE2%.gz}'_first_'$SampledPeakCount
	echo 'Subsampled peak file corresponding to the second input: '$ConvFile2
	#if [ ! -f $ConvFile2 ]; then
		zcat $FILE2 | sort -k8,8nr > $OutDir'/temp2.txt'
		#head -n $SampledPeakCount $OutDir'/temp2.txt' | gzip -c > $ConvFile2
		head -n $SampledPeakCount $OutDir'/temp2.txt' > $ConvFile2
		rm $OutDir'/temp2.txt'
	#fi
else
	#ConvFile2=${FILE2}'_first_'$SampledPeakCount'.gz'
	ConvFile2=${FILE2}'_first_'$SampledPeakCount
	echo 'Subsampled peak file corresponding to the second input: '$ConvFile2
	#if [ ! -f $ConvFile2 ]; then
		cat $FILE2 | sort -k8,8nr > $OutDir'/temp2.txt'
		#head -n $SampledPeakCount $OutDir'/temp2.txt' | gzip -c > $ConvFile2
		head -n $SampledPeakCount $OutDir'/temp2.txt' > $ConvFile2
		rm $OutDir'/temp2.txt'
	#fi
fi

#===========================================

# we employ p value as the measure for rank determination of the peaks
# We note that only narrow peaks are analyzed for the significance test - so the 6th argument is F (no broadpeak)
# we also note that the criteria of peak overlap is set as 1 bp. 
# So the 5th argument is placed as 0 - if it is 0.5, 50% overlap criteria is imposed

# this output directory also notes the settings used for this IDR
CurrOutDir=$OutDir	#'/'	#$IDR_OutFold
mkdir -p $CurrOutDir

# the prefix also contains the output directory where all results will be stored
CurrOutPrefix=$CurrOutDir'/'$PREFIX

# parameter description is provided in the ENCODE IDR documentation
if [ ! -f $CurrOutPrefix'-overlapped-peaks.txt' ]; then
	# command for batch consistency analysis
	# Note: the input to this program should be uncompressed narrow peak file

	# first unzip the files
	#gunzip $ConvFile1
	#file1=${ConvFile1%.gz}
	#gunzip $ConvFile2
	#file2=${ConvFile2%.gz}

	# then call the batch consistency command
	cd $IDRCodeDir
	#Rscript $exec1 $file1 $file2 -1 $CurrOutPrefix 0 F p.value
	Rscript $exec1 $ConvFile1 $ConvFile2 -1 $CurrOutPrefix 0 F p.value
	cd -

	# now re-zip the peak files
	#gzip $file1
	#gzip $file2 

fi

#======================
# add - sourya
# here we call a custom R function
# which plots IDR scatter analysis between this pair of samples

Rscript IDRScatterPlot.r $IDRCodeDir $ConvFile1 $ConvFile2 $IDRCodeDir'/genome_table.txt' $CurrOutPrefix

#----------------------------------
# after generating the IDR statistics for this pair of replicates, 
# now quantify the similarities

# number of peaks in the input peak files
# and also the number of overlapped peaks
#npeak1=`zcat $ConvFile1 | wc -l`
#npeak2=`zcat $ConvFile2 | wc -l`
npeak1=`cat $ConvFile1 | wc -l`
npeak2=`cat $ConvFile2 | wc -l`

Rscript IDRSummary.r $CurrOutPrefix'-overlapped-peaks.txt' $npeak1 $npeak2

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
