#!/bin/bash -ex
#PBS -l nodes=1:ppn=4
#PBS -l mem=10GB
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

#=================================
# this program denotes a sample pipeline for ATAC-seq data
# applicable only a single fastq or alignment file is provided
#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================

# usage info
usage(){
cat << EOF

usage: 

Options:    

  -- required:
	-I  InpFile   		 Input alignment file (Bowtie2 aligned)
	-n  PREFIX           Prefix of output files.
	-d  OutDir 			 Set the output directory which will contain all the results
 	-w 	BigWigGenome	 The reference genome which is used to convert BAM file to a BigWig file. (such as 'hg19', 'mm9', etc.)

EOF
}

while getopts "n:I:d:w:" opt;
do
	case "$opt" in
		n) PREFIX=$OPTARG;;
		I) InpFile=$OPTARG;;
		d) OutDir=$OPTARG;;
		w) BigWigGenome=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done

# executable to convert the sorted bam file to the bigwig format
BigWigCreateExec='/home/sourya/proj/utils/bam_to_bigwig.sh'

#======================
# convert the alignment file to the bigwig data format
# for track visualization
#======================
BigWigoutdir=$OutDir'/Out_BigWig'
mkdir -p $BigWigoutdir

# we use sorted (before duplicate removal) bam file
$BigWigCreateExec -I $InpFile -g $BigWigGenome -d $BigWigoutdir -n $PREFIX

