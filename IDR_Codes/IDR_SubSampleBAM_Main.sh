#!/bin/bash

#=================================
# this script encapsulates IDR analysis between different replicates into a single script
# provided that the inputs are in BAM format
# and they need to be resampled + peak calling + IDR
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: ./IDR_SubSampleBAM_Main.sh [-h] [-I inpfile1.bam] [-I inpfile2.bam] [-d OutDir] [-P PathIDRCode] [-n 'IDR_test'] [-c 25000] [-C control.bam]

Options:    

  -- required:
	-I  InpFile       	 A list of input bam files. At least two bam files are required.
	-d  OutDir 		 	 Output directory (absolute path preferred) which will store the IDR results.
	-P 	PathIDRCode		 Path of the IDRCode package (Kundaje et. al. after its installation)
	-n 	PREFIX 			 Prefix of output files. Default 'IDR_ATAC'.
	-c  CountPeak		 No of peaks in both replicates that will be compared. Default 25000.
	-C  CONTROLBAM		 Control file (in eiher .BAM or tagalign file in .gz format)
EOF
}

# default values of peaks that need to be retained
CountPeak=25000

# code containing the IDR + consistency plot
# dir2='/home/sourya/packages/idrCode/'
exec2='batch-consistency-plot.r'

# default control bam file
CONTROLBAM=""

# default prefix string
PREFIX='IDR_ATAC'

# Sourya - Note the processing of input file argument since it can be more than one file
# Note the change of notations
while getopts "I:n:d:c:C:P:" opt;
do
	case "$opt" in
		I) InpFile+=($OPTARG);;
		n) PREFIX=$OPTARG;;
		d) OutDir=$OPTARG;;
		c) CountPeak=$OPTARG;;
		C) CONTROLBAM=$OPTARG;;
		P) dir2=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $InpFile ]]; then
	echo 'User did not provide any input BAM file - exit for the moment !!'
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
if [ ! -f $OutDir'/Replicate_Names.txt' ]; then
	echo 'Analyzing the '$nsample' Number of replicates --- ' > $OutDir'/Replicate_Names.txt'
	for (( i=0; i<${nsample}; i++ ))
	do
		echo 'Sample '$i' is : '${InpFile[$i]} >> $OutDir'/Replicate_Names.txt'
	done
fi

# loop for pairwise execution of samples
for (( i=0; i<${nsample}-1; i++ ))
do
	for (( j=$i+1; j<${nsample}; j++ ))
	do
		# pair of samples
		sample1=${InpFile[$i]}
		sample2=${InpFile[$j]}
		# execute the sample pairs
		# Note the output directory name
		if [[ ! -z $CONTROLBAM ]]; then
			./IDR_SubSampleBAM.sh -A $sample1 -B $sample2 -d $OutDir'/'$i'_and_'$j -P $dir2 -n $PREFIX -c $CountPeak -C $CONTROLBAM
		else
			./IDR_SubSampleBAM.sh -A $sample1 -B $sample2 -d $OutDir'/'$i'_and_'$j -P $dir2 -n $PREFIX -c $CountPeak
		fi
	done
done

#=================================
# batch consistency plots
#=================================

# the pattern of input prefix present in every replicate 
# depends on the control sample and tagmentation option

if [[ ! -z $CONTROLBAM ]]; then
	inppfx='C1'
else
	inppfx='C0'
fi
inppfx=$inppfx'_Peak'$CountPeak'/'$PREFIX

# basic plotting file name format 
# without the extension '-plot.pdf'
plotfilename='IDR_Batch_Plot'
if [[ ! -z $CONTROLBAM ]]; then
	plotfilename=$plotfilename'_C1'
else
	plotfilename=$plotfilename'_C0'
fi

#if [ ! -f $OutDir'/'$plotfilename'-plot.pdf' ]; then

	# no of pairs of samples
	x=$nsample
	y=`expr $nsample - 1`
	z=`expr $x \* $y`
	npairs=`expr $z / 2`
	echo 'npairs: '$npairs

	# output command for IDR plot
	cmd='Rscript '$exec2' '$npairs' '$OutDir'/'$plotfilename
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
	ps2pdf $OutDir'/'$plotfilename'-plot.ps' $OutDir'/'$plotfilename'-plot.pdf' 

#fi

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------
