#!/bin/bash

#=================================
# sample script for IDR execution
# where peaks generated from multiple ChIP-seq replicates are provided as input
#=================================

# main executable of IDR script
# when peak files are used as input
IDRScript='./IDR_Codes/IDRMain.sh'

# main executable of IDR script
# when BAM files are used as input
IDRScriptBAM='./IDR_Codes/IDR_SubSampleBAM_Main.sh'

#******************************
# path containing the IDRCode package by Anshul Kundaje et. al.
# user should replace this path with their custom installation directory
IDRCodePackage='/home/sourya/packages/idrCode/'
#******************************


#====================
# IDR testing script 1
# examining IDR between two peak files
# top 25K common peaks between two samples are experimented
#====================

SampleBaseDir='/home/sourya/test1/'

$IDRScript -I $SampleBaseDir'Sample1/MACS2_Default_Tag_No_Control/Sample1.macs2_peaks.narrowPeak_Q0.01filt' -I $SampleBaseDir'Sample2/MACS2_Default_Tag_No_Control/Sample2.macs2_peaks.narrowPeak_Q0.01filt' -d $SampleBaseDir'/Sample_IDR_Peaks' -P $IDRCodePackage


#====================
# IDR testing script 2
# examining IDR between two BAM files
# first these BAM files are subsampled
# and their peaks are estimated using MACS2
# top 25K common peaks between two samples are experimented
# no control BAM file is provided
# user may specify one or more control BAM files using -C option
# like -C control1.bam -C control2.bam etc.
#====================

SampleBaseDir='/home/sourya/test2/'

$IDRScriptBAM -I $SampleBaseDir'Sample1/Alignment_MAPQ30/Sample1.align.sort.MAPQ30.bam' -I $SampleBaseDir'Sample2/Alignment_MAPQ30/Sample2.align.sort.MAPQ30.bam' -d $SampleBaseDir'/Sample_IDR_BAMFiles' -P $IDRCodePackage -c 25000
