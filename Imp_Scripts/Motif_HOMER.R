#!/usr/bin/env Rscript

#=========================
# call motifs from ATAC-seq peaks using HOMER
#=========================

library(optparse)
library(data.table)

options(scipen = 999)
options(datatable.fread.datatable=FALSE)

#=================================
# function to create peak summit information
#=================================
GeneratePeakSummitFile <- function(PeakSummitFile, PeakData, offset=500) {
	
	if (ncol(PeakData) > 9) {
		# use the relative peak summit information (10th field) and generate an offset of 
		outDF <- cbind.data.frame(PeakData[,1], (PeakData[,2] + PeakData[,10] - offset), (PeakData[,2] + PeakData[,10] + offset))
	} else {
		# use the midpoint of the peaks as the summit
		outDF <- cbind.data.frame(PeakData[,1], (as.integer((PeakData[,2] + PeakData[,3])/2) - offset), (as.integer((PeakData[,2] + PeakData[,3])/2) + offset))
	}
	write.table(outDF, PeakSummitFile, row.names=F, col.names=F, sep="\t", quote=F, append=F)

}	# end function

#===========================================================
option_list = list(

	make_option(c("--MotifFindExec"), type="character", default=NULL, help="HOMER motif finding executable"),
	make_option(c("--RefGenome"), type="character", default=NULL, help="Reference genome name."),
	make_option(c("--PeakFile"), type="character", default=NULL, help="ATAC-seq Peak file."),
	make_option(c("--PValThr"), type="numeric", default=0, help="Threshold of -log10(p-value) above which peaks will be considered. Default = 0, means no Threshold is imposed."),
	make_option(c("--QValThr"), type="numeric", default=0, help="Threshold of -log10(q-value) above which peaks will be considered. Default = 0, means no threshold is imposed."),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory."),
	make_option(c("--SizeVal"), type="integer", action="store", default=200, help="Size argument of HOMER motif finding. Default = 200"),
	make_option(c("--SummitOffset"), type="integer", action="store", default=500, help="Offset around the peak summit position to be considered for motif finding. Default = 500")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

system(paste("mkdir -p", opt$OutDir))

PValThr <- as.numeric(opt$PValThr)
QValThr <- as.numeric(opt$QValThr)
if (QValThr > 0) {
	PValThr <- 0
}

if ((PValThr == 0) & (QValThr == 0)) {
	# CurrOutDir <- paste0(opt$OutDir, '/Motif_Complete_Peaks_Size_', opt$SizeVal, '_SummitOffset_', opt$SummitOffset)
	CurrOutDir <- paste0(opt$OutDir, '/Motif_Complete_Peaks_SummitOffset_', opt$SummitOffset)
} else if (QValThr > 0) {
	# CurrOutDir <- paste0(opt$OutDir, '/Motif_Peaks_QvalThr_', QValThr, '_Size_', opt$SizeVal, '_SummitOffset_', opt$SummitOffset)
	CurrOutDir <- paste0(opt$OutDir, '/Motif_Peaks_QvalThr_', QValThr, '_SummitOffset_', opt$SummitOffset)
} else {
	# CurrOutDir <- paste0(opt$OutDir, '/Motif_Peaks_PvalThr_', PValThr, '_Size_', opt$SizeVal, '_SummitOffset_', opt$SummitOffset)
	CurrOutDir <- paste0(opt$OutDir, '/Motif_Peaks_PvalThr_', PValThr, '_SummitOffset_', opt$SummitOffset)
}
system(paste("mkdir -p", CurrOutDir))

# read the complete peak data
PeakData <- data.table::fread(opt$PeakFile)

# filter peaks if there is any p-value or q-value specific threshold is provided
# then call the motif finding routine
if (PValThr > 0) {	
	PeakData_Filt <- PeakData[which(PeakData[, 8] > PValThr), ]
	if (nrow(PeakData_Filt) > 0) {		
		# write the filtered peaks
		FiltPeakFileName <- paste0(CurrOutDir, '/Filtered_Peaks_PvalThr.bed')
		write.table(PeakData_Filt, FiltPeakFileName, row.names=F, col.names=F, sep="\t", quote=F, append=F)
		# extract the peak summits and +/- opt$SummitOffset bp from the summits
		FiltPeakFileNameSummit <- paste0(CurrOutDir, '/Filtered_Peaks_PvalThr_Summit_Offset_', opt$SummitOffset, 'bp.bed')
		GeneratePeakSummitFile(FiltPeakFileNameSummit, PeakData_Filt, offset=opt$SummitOffset)
		# now call motif using these summit information
		# currently commented - sourya
		# system(paste(opt$MotifFindExec, FiltPeakFileNameSummit, opt$RefGenome, CurrOutDir, " -size ", opt$SizeVal, " -mask"))
	}		
} else if (QValThr > 0) {	
	PeakData_Filt <- PeakData[which(PeakData[, 9] > QValThr), ]
	if (nrow(PeakData_Filt) > 0) {
		# write the filtered peaks
		FiltPeakFileName <- paste0(CurrOutDir, '/Filtered_Peaks_QvalThr.bed')
		write.table(PeakData_Filt, FiltPeakFileName, row.names=F, col.names=F, sep="\t", quote=F, append=F)
		# extract the peak summits and +/- opt$SummitOffset bp from the summits
		FiltPeakFileNameSummit <- paste0(CurrOutDir, '/Filtered_Peaks_QvalThr_Summit_Offset_', opt$SummitOffset, 'bp.bed')
		GeneratePeakSummitFile(FiltPeakFileNameSummit, PeakData_Filt, offset=opt$SummitOffset)
		# now call motif using these summit information
		# currently commented - sourya
		# system(paste(opt$MotifFindExec, FiltPeakFileNameSummit, opt$RefGenome, CurrOutDir, " -size ", opt$SizeVal, " -mask"))
	}
} else {
	# extract the peak summits and +/- opt$SummitOffset bp from the summits
	FiltPeakFileNameSummit <- paste0(CurrOutDir, '/Peaks_Summit_Offset_', opt$SummitOffset, 'bp.bed')	
	GeneratePeakSummitFile(FiltPeakFileNameSummit, PeakData, offset=opt$SummitOffset)

	# now call motif using these summit information
	# currently commented - sourya
	# system(paste(opt$MotifFindExec, FiltPeakFileNameSummit, opt$RefGenome, CurrOutDir, " -size ", opt$SizeVal, " -mask"))
}




