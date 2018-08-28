#!/bin/bash

#=================================
# this program creates a tag align file (in .gz format) from one or more input aligned (in .bam / .gz format) files
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: 

1) ./TagAlign.sh [-h] [-I inpfile1] [-N 0] [-f 4] [-r 5] [-O Outfile]			For ATAC seq or ChIPMentation data
2) ./TagAlign.sh [-h] [-I inpfile1] [-N 1] [-O Outfile]							For standard ChIP seq data

Example to process multiple input files:
./TagAlign.sh [-h] [-I inpfile1] [-I inpfile2] [-N 1] [-O Outfile]					
	Here all the input files will be processed separately, and their outputs will be combined in a single file

Options:    

  -- required:
	-I  InpFile          Input files (aligned in .bam format) or already in gzipped bed format (.gz). 
						 User can provide multiple input files together.
	-N  NoShift			 It is a binary variable. If 1, the aligned files strand information are not altered. 
						 For standard ChIP seq data, this should be set as 1. 
						 For a ChIPMentation or ATAC seq data, this should be set as 0, since in these cases the tagaligned files are formed 
						 by shifting forward and reverse strands to cover the length of transposon. Default 1.
	-f  fwdshift		 If NoShift=0, this value signifies the amount of shift a forward strand will 
						 require to cover the length of transposon. Default 4.
	-r  revshift         If NoShift=0, this value signifies the amount of shift a reverse strand 
						 will require to cover the length of transposon. Default 5.
	-O  OutFile 		 Output tagalign file (in .gz format) combining the input files
	-q  MAPQ_THR		 Quality threshold that will be applied on the given input BAM file (default 30)
  
EOF
}

# threshold of mapq quality
MAPQ_THR=30

# default configurations
NoShift=1
fwdshift=4
revshift=5

# Sourya - Note the processing of input file argument since it can be more than one file
# Note the change of notations

while getopts "I:N:f:r:O:q:" opt;
do
	case "$opt" in
		I) InpFile+=($OPTARG);;
		N) NoShift=$OPTARG;;
		f) fwdshift=$OPTARG;;
		r) revshift=$OPTARG;;
		O) OutFile=$OPTARG;;
		q) MAPQ_THR=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

echo 'Within utility function TagAlign'

# # this line should be added when processing a list of inputs using the same command line option
# shift $(( OPTIND - 1 ))

# number of input files provided
ninp=${#InpFile[@]}
echo 'Number of input files : '$ninp

if [[ $ninp == 0 ]]; then
	echo 'User should provide one or more aligned (.bam) or existing tagalign (.gz) files to combine them - exit for the moment !!'
	exit 1
fi

# echo 'List of input files: '$InpFile

if [[ -z $OutFile ]]; then
	echo 'User did not provide the output file name - exit for the moment !!'
	exit 1
fi

#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

# also check the extension of input file
filebase1=$(basename "${InpFile[0]}")
if [[ $filebase1 =~ \.bam$ ]]; then
	bamext=1
	echo 'Input files are provided in .bam format'
else
	if [[ $filebase1 =~ \.gz$ ]]; then
		bamext=0
		echo 'Input files are already in .gz format'
	else
		echo 'User should provide either one or more aligned (.bam) files or previously generated TagAlign files in .gz format !! Exit '
		exit 1
	fi
fi

# similarly check the extension of output file and if required, append the gzipped extension
filebase2=$(basename "$OutFile")
if [[ $filebase2 =~ \.gz$ ]]; then
	echo 'User has correctly provided gzipped outfile name'
else
	echo 'Appending gzipped extension in the output file name'
	OutFile=$OutFile'.gz'
fi

# check the number of input files and proceed accordingly
if [ $ninp == 1 ]; then
	# Only one input file is provided
	echo 'Converting the file: '${InpFile[0]}
	if [ $bamext == 1 ]; then
		# here one input bam file is provided
		# so convert the bam file according to the shifting / non-shifting criteria
		if [ $NoShift == 0 ]; then
			samtools view -b -F 1548 -q $MAPQ_THR ${InpFile[0]} | bamToBed -i stdin | awk -v f=$fwdshift -v r=$revshift 'BEGIN {OFS = "\t"} ; function pos(x){return ((x < 0.0) ? 0 : x)}  {if ($6 == "+") print $1, $2 + f, $3 + f, $4, $5, $6; else print $1, pos($2 - r), pos($3 - r), $4, $5, $6}' | gzip -c > $OutFile
		else
			samtools view -b -F 1548 -q $MAPQ_THR ${InpFile[0]} | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $OutFile
		fi
	else
		# already a gzipped file is provided
		# we can just copy the file
		cp ${InpFile[0]} $OutFile
	fi
else
	# Multiple input files are provided
	if [ $bamext == 1 ]; then
		# input files are provided in bam format
		# we have to convert them individually, and then combine them
		convfilelist=''
		for (( i=0; i<${ninp}; i++ ));
		do
			# convert the current file into a temporary output file
			echo 'Converting the file: '${InpFile[i]}
			curroutfile='temp_'$i'.gz'
			if [ $NoShift == 0 ]; then
				samtools view -b -F 1548 -q $MAPQ_THR ${InpFile[i]} | bamToBed -i stdin | awk -v f=$fwdshift -v r=$revshift 'BEGIN {OFS = "\t"} ; function pos(x){return ((x < 0.0) ? 0 : x)}  {if ($6 == "+") print $1, $2 + f, $3 + f, $4, $5, $6; else print $1, pos($2 - r), pos($3 - r), $4, $5, $6}' | gzip -c > $curroutfile
			else
				samtools view -b -F 1548 -q $MAPQ_THR ${InpFile[i]} | bamToBed -i stdin | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $curroutfile
			fi
			# also update the command for combining these generated files
			convfilelist=$convfilelist' '$curroutfile
		done
		zcat $convfilelist | gzip -c > $OutFile
		# remove temporary files
		for (( i=0; i<${ninp}; i++ ));
		do
			rm 'temp_'$i'.gz'
		done
	else
		# input files are already in gzipped format
		# we can just combine them
		convfilelist=''
		for val in "${InpFile[@]}"; do
			convfilelist=$convfilelist' '$val
		done
		zcat $convfilelist | gzip -c > $OutFile
	fi
fi


#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------

