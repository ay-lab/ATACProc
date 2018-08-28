#!/bin/bash 

#=================================
# this script encapsulates IDR analysis between different replicates into a single script
# It has input peak files (2 or more) 
# The script calls IDRAnalysis.sh for pairwise analysis
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: ./IDRMain.sh [-h] [-I peakfile1.narrowpeak] [-I peakfile2.narrowpeak] [-P PathIDRCode] [-d OutDir] [-n PREFIXSTR]
Example:
	./IDRMain.sh -I peak1.narrowPeak -I peak2.narrowPeak -I peak3.narrowPeak -P /home/sourya/packages/idrCode/ -d /home/sourya/OutDir_IDR -n 'IDR_test'

Options:    

  -- required:
	-I  InpFile        	 A list of input peak files (obtained from MACS2 - in .narrowPeak or .narrowPeak.gz format). 
						 At least two peak files are required.
	-P 	PathIDRCode		 Path of the IDRCode package (Kundaje et. al. after its installation)
	-d  OutDir 		 	 Output directory (absolute path preferred) which will store the IDR results.
	-n 	PREFIX 			 Prefix of output files. Default 'IDR_ATAC'.
EOF
}

# default variables and values
IDR_code='./IDRAnalysis.sh'
PREFIX='IDR_ATAC'

# code containing the IDR + consistency plot
# dir2='/home/sourya/packages/idrCode/'
exec2='batch-consistency-plot.r'

# Sourya - Note the processing of input file argument since it can be more than one file
# Note the change of notations
while getopts "I:n:d:P:" opt;
do
	case "$opt" in
		I) InpFile+=($OPTARG);;
		n) PREFIX=$OPTARG;;
		d) OutDir=$OPTARG;;
		P) dir2=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $InpFile ]]; then
	echo 'User did not provide any input peak file - exit for the moment !!'
	exit 1
fi

if [[ -z $dir2 ]]; then
	echo 'User did not provide the path of IDRCode package (Kundaje et. al.) - exit for the moment !!'
	exit 1
fi

if [[ -z $OutDir ]]; then
	echo 'User did not provide output directory for storing the results - exit for the moment !!'
	exit 1
fi

# number of input files provided
nsample=${#InpFile[@]}
echo 'Number of input files : '$nsample

if [ $nsample -lt 2 ]; then
	echo 'User needs to provide at least two peak files for comparison - exit for the moment !!'
	exit 1
fi

# generate the output directory
mkdir -p $OutDir

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

#=================================
# batch replicate analysis
#=================================
# if [ ! -f $OutDir'/Replicate_Names.txt' ]; then
	echo 'Analyzing the '$nsample' Number of replicates --- ' > $OutDir'/Replicate_Names.txt'
	for (( i=0; i<${nsample}; i++ ))
	do
		echo 'Sample '$i' is : '${InpFile[$i]} >> $OutDir'/Replicate_Names.txt'
	done
# fi

#===================
# add - sourya
# here we analyze individual peak files (input)
# and accordingly assign the number of common peaks to be considered
PeakStatFile=$OutDir'/Input_Peak_Statistics.txt'

echo 'Summarizing the peak count statistics for individual input files: ' > $PeakStatFile

# first get the minimum no of peaks across all the samples
for (( i=0; i<${nsample}; i++ ))
do
	peakfile=${InpFile[$i]}
	pc=`cat $peakfile | wc -l`
	echo "Analyzing the peak file: $peakfile " >> $PeakStatFile
	echo "Peak count: $pc " >> $PeakStatFile
	if [ $i == 0 ]; then
		minpc=$pc
	else
		if [ $minpc > $pc ]; then
			minpc=$pc
		fi
	fi
done

# assign the minimum number of peaks for consideration
if [[ $minpc -gt 200000 ]]; then
	CountPeak=150000
elif [[ $minpc -gt 150000 ]]; then
	CountPeak=100000
elif [[ $minpc -gt 100000 ]]; then
	CountPeak=75000
elif [[ $minpc -gt 75000 ]]; then
	CountPeak=50000
else
	CountPeak=25000
fi

echo "Value of CountPeak (number of common peaks to be analyzed for all replicates): $CountPeak " >> $PeakStatFile
#===================

# loop for pairwise execution of samples
for (( i=0; i<${nsample}-1; i++ ))
do
	for (( j=$i+1; j<${nsample}; j++ ))
	do
		# pair of samples
		sample1=${InpFile[$i]}
		sample2=${InpFile[$j]}
		# execute the sample pairs
		# Note the output directory name - it is the sample directory plus the pairwise comparison
		$IDR_code -a $sample1 -b $sample2 -P $dir2 -d $OutDir'/'$i'_and_'$j -n $PREFIX -c $CountPeak
	done
done

#=================================
# batch consistency plots
#=================================
if [ ! -f $OutDir'/IDR_Batch_Plot-plot.pdf' ]; then

	# the pattern of input prefix present in every replicate 
	#inppfx=$IDR_OutFold'/'$PREFIX
	inppfx=$PREFIX

	# no of pairs of samples
	x=$nsample
	y=`expr $nsample - 1`
	z=`expr $x \* $y`
	npairs=`expr $z / 2`
	echo 'npairs: '$npairs

	# output command for IDR plot
	cmd='Rscript '$exec2' '$npairs' '$OutDir'/IDR_Batch_Plot'
	for (( i=0; i<${nsample}-1; i++ ))
	do
		for (( j=$i+1; j<${nsample}; j++ ))
		do
			cmd=$cmd' '$OutDir'/'$i'_and_'$j'/'$inppfx
		done
	done
	echo 'cmd: '$cmd

	# execute the command
	# first go to the directory containing the R code of the IDR
	cd $dir2
	$cmd
	cd -

	# now convert the generated postscript plot file to a pdf file
	ps2pdf $OutDir'/IDR_Batch_Plot-plot.ps' $OutDir'/IDR_Batch_Plot-plot.pdf' 

fi
   
#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------


