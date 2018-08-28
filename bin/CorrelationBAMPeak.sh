#!/bin/bash 

#=================================
# this program is a supporting script for the ATAC seq pipeline

# inputs:
# 1) a set of input bam files (requires sorted bam files, and possibly indexed as well)
# 2) a set of peak files
# 3) max no of peaks to be considered

# union of the given peak files is used to compute the coverage with respect to individual input bam files
# with respect to a minimum coverage (num of reads) threshold
# and a given threshold of the max  no of peaks to be considered
# bam files are subsampled to cover only the subset of union peak set
# correlation between these subsampled bam files 

#=================================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================================


# usage info
usage(){
cat << EOF

Options:    

  -- required:
  	-B 	BAM 			 One or more input bam files (sorted)
  	-P 	NarrowPeak 		 One or more input narrow peak files (corresponding to the input bam files)
  						 and in the same order as the input bam files
	-L  Labels 			 One or more strings (labels) corresponding to input bam files.
 	-D 	OutDir 			 Output directory storing the correlation results
 	-r 	ReadCount		 A threshold (integer) of the number of reads (coverage) that each peak should 
 					 	 minimally cover. Default = 5 (according to the Greenleaf 2018 paper)
 	-c 	PeakCount		 Number of peaks to be randomly from the union set of peaks. 
 						 Default = 50000 (according to the Greenleaf 2018 paper)
	-O 	Overwrite		 this boolean option signifies whether existing output files would 
						 be overwritten (1) or not (0).
						 Default = 0

EOF
}

# default minimum coverage threshold for each peak
ReadCountThr=5

# default threshold of the number of peaks
PeakCountThr=50000

# default output directory
OutDir=`pwd`'/'

# this boolean option signifies whether existing output 
# files would be overwritten (1) or not (0).
# Default = 0
Overwrite=0

while getopts "B:P:D:r:c:O:L:" opt;
do
	case "$opt" in
		B) BAMFILES+=($OPTARG);;	# one or more bam input files can be provided
		P) PEAKFILES+=($OPTARG);;	# one or more peak input files can be provided
		L) Labels+=($OPTARG);;		# labels corresponding to individual input bam files
		D) OutDir=$OPTARG;;
		r) ReadCountThr=$OPTARG;;
		c) PeakCountThr=$OPTARG;;
		O) Overwrite=$OPTARG;;
		\?) usage
			echo "error: unrecognized option -$OPTARG";
			exit 1
			;;
	esac
done


#----------------------------------
# important - sourya
# change the current directory as the dir containing this executable
# since other source files relative to the current directory needs to be called
current_dir=$(pwd)
script_dir=$(dirname $0)
cd $script_dir
#----------------------------------

nbamfiles=${#BAMFILES[@]}
echo 'number of bam files provided: '$nbamfiles

npeakfiles=${#PEAKFILES[@]}
echo 'number of peak files provided: '$npeakfiles

nlabels=${#Labels[@]}
echo 'number of labels provided: '$npeakfiles

# if [[ $nbamfiles != $npeakfiles ]]; then
# 	echo "Number of input bam files and the number of peak files do not match - return !!!"
# 	exit 1
# fi

# check if the input bam files are all indexed
# otherwise index the bam files
listbamfiles=''
for (( i=0; i<${nbamfiles}; i++ ));
do 
	currbamfile=${BAMFILES[i]}
	if [[ $i == 0 ]]; then
		listbamfiles=$currbamfile
	else
		listbamfiles=$listbamfiles' '$currbamfile
	fi
	echo 'processing the bam file index: '$i'  name: '$currbamfile
	if [ ! -f $currbamfile'.bai' ]; then
		samtools index $currbamfile
	fi
done
echo 'listbamfiles: '$listbamfiles

# list of labels
if [[ $nlabels == $nbamfiles ]]; then
	listlabels=''
	# also required a colon separated list 
	listlabelsRscript=''
	for (( i=0; i<${nlabels}; i++ ));
	do
		if [[ $i == 0 ]]; then
			listlabels=${Labels[i]}
			listlabelsRscript=${Labels[i]}
		else
			listlabels=$listlabels' '${Labels[i]}
			listlabelsRscript=$listlabelsRscript':'${Labels[i]}
		fi
	done
	echo 'listlabels: '$listlabels
	echo 'listlabelsRscript: '$listlabelsRscript
fi

# list of peak files (if provided)
if [[ $npeakfiles == $nbamfiles ]]; then
	listpeakfiles=''
	for (( i=0; i<${npeakfiles}; i++ ));
	do 
		echo 'processing the peak file index: '$i'  name: '${PEAKFILES[i]}
		if [[ $i == 0 ]]; then
			listpeakfiles=${PEAKFILES[i]}
		else
			listpeakfiles=$listpeakfiles' '${PEAKFILES[i]}
		fi
	done
	echo 'listpeakfiles: '$listpeakfiles
fi

#=============================
# two cases:
# 1) when peak files are provided, and the subsampled union of peaks are used for correlation
# 2) or, when peaks are not provided, and whole bam files are used for correlation
#==============================
if [[ $npeakfiles == $nbamfiles ]]; then

	# union of the input peak files
	UnionPeakFile=$OutDir'/Union_Peaks_Original.bed'

	echo '***** before computing '$UnionPeakFile'  *****'

	if [[ ! -f $UnionPeakFile || $Overwrite == 1 ]]; then
		cat $listpeakfiles | cut -f1-3 | sort -k1,1 -k2,2n | mergeBed -i stdin > $UnionPeakFile
	fi

	echo '***** after computing '$UnionPeakFile'  *****'	

	# now perform the coverage of this union peaks
	# with respect to input bam files
	coverageoutfile=$OutDir'/Union_Peaks_Original_CoverageVal.bed'

	echo '***** before computing '$coverageoutfile'  *****'

	if [[ ! -f $coverageoutfile || $Overwrite == 1 ]]; then
		# earlier command - sourya
		# bedtools multicov -bams ${listbamfiles} -bed ${UnionPeakFile} > ${coverageoutfile}
		
		# modified command - sourya
		# we found that supplying all of the bam files together in the "multicov" function
		# results errors, probably due to mismatching headers in different bam files
		# so we supply one bam file at a time,
		# compute the coverage with respect to individual bam files
		# and sequentially merge all the information
		for (( i=0; i<${nbamfiles}; i++ ));
		do
			# temp output file
			tempout=$OutDir'/temp_out_union_peaks_coverage.bed'
			bedtools multicov -bams ${BAMFILES[i]} -bed ${UnionPeakFile} > $tempout
			if [ $i == 0 ]; then
				# first iteration - rename the temporary output file to the 
				# final output file
				mv $tempout $coverageoutfile
			else
				# subsequent iteration
				# use bedtools map function to merge the existing contents 
				# of "coverageoutfile" with the new "tempout" contents
				tempout2=$OutDir'/temp_out_union_peaks_coverage2.bed'
				bedtools map -c 4 -o mean -null '0'	-a $coverageoutfile -b $tempout > $tempout2
				# remove the old instance of coverage output file
				# and use the newly constructed file
				rm $coverageoutfile
				mv $tempout2 $coverageoutfile
			fi
		done
		# remove the temporary files
		if [ -f $tempout ]; then
			rm $tempout
		fi
		if [ -f $tempout2 ]; then
			rm $tempout2
		fi
	fi

	echo '***** after computing '$coverageoutfile'  *****'

	coverageoutfileThr=$OutDir'/Union_Peaks_CoverageVal_MinReadThr.bed'

	echo '***** before computing '$coverageoutfileThr'  *****'

	# select only those peaks which have coverage value >= ReadCountThr
	# for all the bam files considered
	if [[ ! -f $coverageoutfileThr || $Overwrite == 1 ]]; then
		awk -v T="$ReadCountThr" -v N="$npeakfiles" '{f=0; for (i=0;i<N;i++) {if ($(NF-i)<T) {f=1}}; {if (f==0) {print $0}}}' $coverageoutfile > $coverageoutfileThr
	fi
 
	echo '***** after computing '$coverageoutfileThr'  *****'

	# subset of the peaks (randomly selected)
	# such that the total number of peaks = PeakCountThr
	# check if the coverge thresholded peaks have higher number of peaks
	# than the mentioned threshold
	coverageoutfileThrSubSample=$OutDir'/Union_Peaks_CoverageVal_MinReadThr_Subsampled_'$PeakCountThr'.bed'

	echo '***** before computing '$coverageoutfileThrSubSample'  *****'

	if [[ ! -f $coverageoutfileThrSubSample || $Overwrite == 1 ]]; then
		npeakAboveThr=`cat $coverageoutfileThr | wc -l`
		echo 'npeakAboveThr: '$npeakAboveThr
		if [[ $npeakAboveThr -gt $PeakCountThr ]]; then
			echo 'Above the mentioned peak count threshold - random subset'
			shuf $coverageoutfileThr | head -n $PeakCountThr | cut -f1-3 | sort -k1,1 -k2,2n > $coverageoutfileThrSubSample
		else
			cat $coverageoutfileThr | cut -f1-3 | sort -k1,1 -k2,2n > $coverageoutfileThrSubSample 
		fi
	fi

	echo '***** after computing '$coverageoutfileThrSubSample'  *****'

	# dumping intermediate results
	# minimum mapping quality is maintained at 30
	OutDumpFile=$OutDir'/results_SubsampledPeak.npz'

	echo '***** before computing '$OutDumpFile'  *****'

	# --minMappingQuality 30 

	if [[ $nlabels == $nbamfiles ]]; then
		multiBamSummary BED-file --BED $coverageoutfileThrSubSample --bamfiles $listbamfiles --labels $listlabels -out $OutDumpFile
	else
		multiBamSummary BED-file --BED $coverageoutfileThrSubSample --bamfiles $listbamfiles --smartLabels -out $OutDumpFile 
	fi

	echo '***** after computing '$OutDumpFile'  *****'

	# using spearman correlation
	
	# OutPlotFile=$OutDir'/Correlation_Spearman_SubsampledPeak.pdf'
	# OutMatrixFile=$OutDir'/Correlation_Spearman_SubsampledPeak.matrix'
	# --plotFile $OutPlotFile --corMethod spearman  --outFileCorMatrix $OutMatrixFile

	# --skipZeros --removeOutliers
	
	OutPlotFileHeatMap=$OutDir'/Correlation_Spearman_SubsampledPeak_Heatmap.pdf'
	OutPlotFileScatter=$OutDir'/Correlation_Spearman_SubsampledPeak_Scatterplot.pdf'
	if [[ ! -f $OutPlotFileHeatMap || ! -f $OutPlotFileScatter || $Overwrite == 1 ]]; then
		plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileHeatMap --corMethod spearman --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix $OutDir'/Correlation_Spearman_SubsampledPeak_Corr.mat' 
		plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileScatter --corMethod spearman --whatToPlot scatterplot 
	fi

	# # using pearson correlation
	# commented - check https://www.biostars.org/p/195328/
	
	# # OutPlotFile=$OutDir'/Correlation_Pearson_SubsampledPeak.pdf'
	# # OutMatrixFile=$OutDir'/Correlation_Pearson_SubsampledPeak.matrix'
	# # --plotFile $OutPlotFile --corMethod pearson  --outFileCorMatrix $OutMatrixFile 

	# # --skipZeros --removeOutliers

	# OutPlotFileHeatMap=$OutDir'/Correlation_Pearson_SubsampledPeak_Heatmap.pdf'
	# OutPlotFileScatter=$OutDir'/Correlation_Pearson_SubsampledPeak_Scatterplot.pdf'

	# if [[ ! -f $OutPlotFileHeatMap || ! -f $OutPlotFileScatter || $Overwrite == 1 ]]; then
	# 	plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileHeatMap --corMethod pearson --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix $OutDir'/Correlation_Pearson_SubsampledPeak_Corr.mat' 
	# 	plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileScatter --corMethod pearson --whatToPlot scatterplot 
	# fi

	#========================
	# now extract the peak intervals (after subsampling)
	# and get the peak intensities for individual input peak files
	OnlyPeakSubsample=$OutDir'/Peak_PVal_Subsampled_'$PeakCountThr'.bed'
	cat $coverageoutfileThrSubSample | cut -f1-3 > $OnlyPeakSubsample

	# merge with individual peak input files
	for (( i=0; i<${npeakfiles}; i++ ));
	do
		tempfile=$OutDir'/Peak_PVal_temp.bed'
		# 8th field in the peak file contains -log10(p) score
		bedtools map -c 8 -o mean -null '0'	-a $OnlyPeakSubsample -b ${PEAKFILES[i]} > $tempfile
		rm $OnlyPeakSubsample
		mv $tempfile $OnlyPeakSubsample
	done

	# now call a R script which would plot the correlation
	# for these peaks
	if [[ $nlabels == $nbamfiles ]]; then
		Rscript CorrelationPeakPlot.r --InpPeakFile $OnlyPeakSubsample --OutDir $OutDir --InpLabels $listlabelsRscript
	else
		Rscript CorrelationPeakPlot.r --InpPeakFile $OnlyPeakSubsample --OutDir $OutDir
	fi
	#========================

else

	# here no peak files are provided
	# so simple correlation using the whole bam files is required

	# dumping intermediate results
	# minimum mapping quality is maintained at 30
	# --minMappingQuality 30 

	OutDumpFile=$OutDir'/results.npz'
	if [[ $nlabels == $nbamfiles ]]; then
		multiBamSummary bins --bamfiles $listbamfiles --labels $listlabels -out $OutDumpFile
	else
		multiBamSummary bins --bamfiles $listbamfiles --smartLabels -out $OutDumpFile 
	fi

	# using spearman correlation
	
	# OutPlotFile=$OutDir'/Correlation_Spearman_SubsampledPeak.pdf'
	# OutMatrixFile=$OutDir'/Correlation_Spearman_SubsampledPeak.matrix'
	# --plotFile $OutPlotFile --corMethod spearman  --outFileCorMatrix $OutMatrixFile

	# --skipZeros --removeOutliers

	OutPlotFileHeatMap=$OutDir'/Spearman_Heatmap.pdf'
	OutPlotFileScatter=$OutDir'/Spearman_Scatterplot.pdf'

	if [[ ! -f $OutPlotFileHeatMap || ! -f $OutPlotFileScatter || $Overwrite == 1 ]]; then
		plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileHeatMap --corMethod spearman --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix $OutDir'/Correlation_Spearman_Corr.mat' 
		plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileScatter --skipZeros --removeOutliers --corMethod spearman --whatToPlot scatterplot 
	fi



	# # using pearson correlation
	# commented - check https://www.biostars.org/p/195328/

	# # OutPlotFile=$OutDir'/Correlation_Pearson_SubsampledPeak.pdf'
	# # OutMatrixFile=$OutDir'/Correlation_Pearson_SubsampledPeak.matrix'
	# # --plotFile $OutPlotFile --corMethod pearson  --outFileCorMatrix $OutMatrixFile 

	# OutPlotFileHeatMap=$OutDir'/Pearson_Heatmap.pdf'
	# OutPlotFileScatter=$OutDir'/Pearson_Scatterplot.pdf'

	# if [[ ! -f $OutPlotFileHeatMap || ! -f $OutPlotFileScatter || $Overwrite == 1 ]]; then
	# 	plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileHeatMap --corMethod pearson --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix $OutDir'/Correlation_Pearson_Corr.mat' 
	# 	plotCorrelation --corData $OutDumpFile --plotFile $OutPlotFileScatter --skipZeros --removeOutliers --corMethod pearson --whatToPlot scatterplot 
	# fi



fi


#----------------------------------
# important - sourya
# now restore the original directory
cd $current_dir
#----------------------------------

