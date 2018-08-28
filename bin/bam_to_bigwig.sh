#!/bin/bash

#========================================
# sample script for converting input bam file to bigwig format

# author: Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#========================================

# usage info
usage(){
cat << EOF

usage: 
./bam_to_bigwig.sh [-h] [-I InpFile] [-g refgenome] [-d OutDir]
Example:
./bam_to_bigwig.sh -I Inp.bam -g 'hg19' -d '/home/sample_ATAC'

Options:
  -- required:
	-I  InpFile          Input BAM file.
	-g  refgenome   	 Reference genome for chromosome size etc.
	-d  OutDir 			 Output directory which will contain all the bigwig file and other data
	-n  OutFilePrefix	 If specified, output bigwig file name will be 'OutFilePrefix.bw' under the directory 'OutDir'
EOF
}

# default output directory
OutDir=`pwd`'/'

# initialization of the prefix string
OutFilePrefix=""

while getopts "I:g:d:n:" opt;
do
	case "$opt" in
		I) InpFile=$OPTARG;;
		g) refgenome=$OPTARG;;
		d) OutDir=$OPTARG;;
		n) OutFilePrefix=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

if [[ -z $InpFile ]]; then
	echo 'No input BAM file is provided - exit !!'
	exit 1
fi

if [[ -z $refgenome ]]; then
	echo 'No reference genome is provided - exit !!'
	exit 1
fi

filebase1=$(basename "${InpFile}")
if [[ $filebase1 =~ \.bam$ ]]; then
	echo 'Input files are provided in .bam format'
else
	echo 'Input file is not in BAM format - exit !!'
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

if [ ! -f $refgenome'.chrom.sizes' ]; then
	# this utility program from UCSC, fetches the chromosome size of the target genome
	# and stores that in the specified text file
	echo 'Getting the chromosome size'
	fetchChromSizes $refgenome > $refgenome'.chrom.sizes'
fi

# convert the bam file to a bedgraph file
# ensure that the bedgraph file contains only valid chromosomes
# if [ ! -f $OutDir'/Inp.bedGraph' ]; then
	genomeCoverageBed -bga -ibam $InpFile -g $refgenome'.chrom.sizes' | awk '( $1 ~ /^chr([1-9]|2[0-2]|1[0-9]|X|M|Y)$/ )' - > $OutDir'/Inp.bedGraph'
# fi

# sort the generated bedgraph file using the utility of UCSC genome browser
# if [ ! -f $OutDir'/Inp.Sorted.bedGraph' ]; then
	bedSort $OutDir'/Inp.bedGraph' $OutDir'/Inp.Sorted.bedGraph'
# fi

# from the bedgraph file, generate the BigWig file
# using an utility of UCSC genome browser
if [[ -z $OutFilePrefix ]]; then
	outbigwigfile=$OutDir'/Inp_BigWig.bw'
else
	outbigwigfile=$OutDir'/'$OutFilePrefix'.bw'
fi 

# if [ ! -f $outbigwigfile ]; then
	bedGraphToBigWig $OutDir'/Inp.Sorted.bedGraph' $refgenome'.chrom.sizes' $outbigwigfile
# fi

#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------

