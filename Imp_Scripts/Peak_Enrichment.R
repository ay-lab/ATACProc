#!/usr/bin/env Rscript

#=========================
# analyze ATAC-seq peaks and plot the enrichment for peaks and surrounding regions
# also analyze separately for promoter peaks and enhancer peaks
# using reference TSS information
#=========================

suppressMessages(library(GenomicRanges))
library(optparse)
library(data.table)

options(scipen = 999)
options(datatable.fread.datatable=FALSE)

#=================================
# function to create peak summit information
#=================================
GeneratePeakSummitFile <- function(PeakSummitFile, PeakData) {
	
	if (ncol(PeakData) > 9) {
		# use the relative peak summit information (10th field) and generate an offset of 
		outDF <- cbind.data.frame(PeakData[,1], (PeakData[,2] + PeakData[,10] - 5), (PeakData[,2] + PeakData[,10] + 5))
	} else {
		# use the midpoint of the peaks as the summit
		outDF <- cbind.data.frame(PeakData[,1], (as.integer((PeakData[,2] + PeakData[,3])/2) - 5), (as.integer((PeakData[,2] + PeakData[,3])/2) + 5))
	}
	write.table(outDF, PeakSummitFile, row.names=F, col.names=F, sep="\t", quote=F, append=F)

}	# end function


#=================================
# function to compute overlap of 1D bins
#=================================
Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {

	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(ov1[,1])
		ov_idx_file2 <- unique(ov1[,2])		
	} else {
		ov_idx_file1 <- ov1[,1]
		ov_idx_file2 <- ov1[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}

#=================================
# function to plot the heatmap using deeptools
#=================================
PlotHeatMap <- function(CurrOutDir, DeepToolsDir, outmatfile, Label) {

	# then use this matrix to plot profile
	outprofileplotfile <- paste0(CurrOutDir, '/out_mat_profile_plot.pdf')
	outprofileplotfile1 <- paste0(CurrOutDir, '/out_mat_profile_plot_1.pdf')
	outprofileplotfile2 <- paste0(CurrOutDir, '/out_mat_heatmap_plot.pdf')
	outprofileplotfile3 <- paste0(CurrOutDir, '/out_mat_heatmap_plot_1.pdf')

	system(paste0(DeepToolsDir, "/plotProfile --matrixFile ", outmatfile, " --outFileName ", outprofileplotfile, " --plotHeight 7 --plotWidth 10 --samplesLabel ", Label, " --plotTitle ATACPeakTSSEnrichment --plotFileFormat pdf --yMin 0.5 --yMax 40 --colors red yellow blue"))

	system(paste0(DeepToolsDir, "/plotProfile --matrixFile ", outmatfile, " --outFileName ", outprofileplotfile1, " --plotHeight 7 --plotWidth 10 --samplesLabel ", Label, " --plotTitle ATACPeakTSSEnrichment --plotFileFormat pdf --yMin 0.5 --yMax 5 --colors red yellow blue"))

	system(paste0(DeepToolsDir, "/plotHeatmap --matrixFile ", outmatfile, " --outFileName ", outprofileplotfile2, " --heatmapHeight 10 --heatmapWidth 8 --samplesLabel ", Label, " --plotTitle ATACPeakTSSEnrichment --plotFileFormat pdf --yMin 0.5 --yMax 40 --zMin 0 --zMax 50"))

	system(paste0(DeepToolsDir, "/plotHeatmap --matrixFile ", outmatfile, " --outFileName ", outprofileplotfile3, " --heatmapHeight 10 --heatmapWidth 8 --samplesLabel ", Label, " --plotTitle ATACPeakTSSEnrichment --plotFileFormat pdf --yMin 0.5 --yMax 5 --zMin 0 --zMax 50"))

}	# end function


#===========================================================
option_list = list(

	make_option(c("--BigWigFile"), type="character", default=NULL, help="BigWig file of ATAC-seq data."),
	make_option(c("--Label"), type="character", default=NULL, help="Label or sample name of ATAC-seq data."),
	make_option(c("--DeepToolsDir"), type="character", default=NULL, help="Deeptools executable directory."),
	make_option(c("--TSSFile"), type="character", default=NULL, help="File containing reference genome TSS information."),
	make_option(c("--PeakFile"), type="character", default=NULL, help="File containing ATAC-seq peak information."),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory."),
	make_option(c("--Offset"), type="integer", action="store", default=5000, help="Offset with respect to summit (in bp) to compute enrichment. Default = 5000 means 5 Kb around peak summits would be used for enrichment.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

system(paste("mkdir -p", opt$OutDir))

# read the input peaks
PeakData <- data.table::fread(opt$PeakFile)

# extract the peak summits of the complete peak file
CurrOutDir <- paste0(opt$OutDir, '/Complete_Peaks')
system(paste("mkdir -p", CurrOutDir))

PeakSummitFile <- paste0(CurrOutDir, '/Peak_Summits.bed')
if (file.exists(PeakSummitFile) == FALSE) {
	GeneratePeakSummitFile(PeakSummitFile, PeakData)
}

# now apply deeptools utility to compute enrichment
outmatfile <- paste0(CurrOutDir, '/deeptools_out_mat_TSS.gz')
if (file.exists(outmatfile) == FALSE) {
	system(paste0(opt$DeepToolsDir, "/computeMatrix reference-point -R ", PeakSummitFile, " -S ", opt$BigWigFile, " -a ", opt$Offset, " -b ", opt$Offset, " --skipZeros --outFileName ", outmatfile))
}
PlotHeatMap(CurrOutDir, opt$DeepToolsDir, outmatfile, opt$Label)

# if TSS information is also provided, find the enrichment of promoter and 
# enhancer peaks separately
if (!is.null(opt$TSSFile)) {
	TSSData <- data.table::fread(opt$TSSFile)
	# 2.5 Kb overlap on both side of TSS data
	ov <- Overlap1D(PeakData[,1:3], cbind.data.frame(TSSData[,1:2],TSSData[,2]), boundary=0, offset=2500, uniqov=TRUE)
	PromPeakData <- PeakData[ov$A_AND_B, ]
	EnhPeakData <- PeakData[ov$A_MINUS_B, ]

	# process the promoter peaks
	if (nrow(PromPeakData) > 0) {
		CurrOutDir <- paste0(opt$OutDir, '/Promoter_Peaks')
		system(paste("mkdir -p", CurrOutDir))
		PeakSummitFile <- paste0(CurrOutDir, '/Peak_Summits.bed')
		if (file.exists(PeakSummitFile) == FALSE) {
			GeneratePeakSummitFile(PeakSummitFile, PromPeakData)
		}
		# now apply deeptools utility to compute enrichment
		outmatfile <- paste0(CurrOutDir, '/deeptools_out_mat_TSS.gz')
		if (file.exists(outmatfile) == FALSE) {
			system(paste0(opt$DeepToolsDir, "/computeMatrix reference-point -R ", PeakSummitFile, " -S ", opt$BigWigFile, " -a ", opt$Offset, " -b ", opt$Offset, " --skipZeros --outFileName ", outmatfile))
		}
		PlotHeatMap(CurrOutDir, opt$DeepToolsDir, outmatfile, opt$Label)
	}

	# process the enhancer peaks
	if (nrow(EnhPeakData) > 0) {
		CurrOutDir <- paste0(opt$OutDir, '/Enhancer_Peaks')
		system(paste("mkdir -p", CurrOutDir))
		PeakSummitFile <- paste0(CurrOutDir, '/Peak_Summits.bed')
		if (file.exists(PeakSummitFile) == FALSE) {
			GeneratePeakSummitFile(PeakSummitFile, EnhPeakData)
		}
		# now apply deeptools utility to compute enrichment
		outmatfile <- paste0(CurrOutDir, '/deeptools_out_mat_TSS.gz')
		if (file.exists(outmatfile) == FALSE) {
			system(paste0(opt$DeepToolsDir, "/computeMatrix reference-point -R ", PeakSummitFile, " -S ", opt$BigWigFile, " -a ", opt$Offset, " -b ", opt$Offset, " --skipZeros --outFileName ", outmatfile))
		}
		PlotHeatMap(CurrOutDir, opt$DeepToolsDir, outmatfile, opt$Label)
	}
}



