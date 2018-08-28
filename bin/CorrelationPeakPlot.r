#!/usr/bin/env Rscript

#===========================================================
# R script for plotting the pairwise correlation of the peak intensity values
# for a given pair of input peaks
# Input: one file containing the peak intervals (first three columns)
# and the peak intensity values in the subsequent columns
# another input is the labels of the given samples

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
# February 26, 2018
#===========================================================

# used for parsing the command line arguments
library(optparse)

# library(ggplot2)

# plot dimension values
PlotWidth <- 10
PlotHeight <- 7

# font size used in texts
FontSize <- 20

# different colors used in heatmap
ColorVec <- c('blue', 'cyan', 'green', 'yellow', 'orange', 'red')

option_list = list(
	make_option(c("--InpPeakFile"), type="character", default=NULL, help="Input file containing peak locations and the peak intensity values for all the candidate input samples"),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory containing the results"),	
	make_option(c("--InpLabels"), type="character", default=NULL, help="Comma or colon separated list of labels associated with individual samples (default %default)")

); 

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

# read the input peak file
InpData <- read.table(opt$InpPeakFile, header=F)

# number of samples is the number of columns of the input file
# minus the first three columns
NumSample <- ncol(InpData) - 3

# read the labels of the input samples 
# if not provided, assign numeric labels 1 to NumSample
if (is.null(opt$InpLabels)) {
	InpLabelList <- as.character(seq(1, NumSample))
} else {
	InpLabelList <- as.character(unlist(strsplit(opt$InpLabels,"[,:]")))
}

if (is.null(opt$OutDir)) {
	OutDir <- getwd()
} else {
	OutDir <- opt$OutDir
}

cat(sprintf("\n\n *** NumSample: %s ", NumSample))
cat(sprintf("\n\n *** InpLabelList: %s ", InpLabelList))

TextFile <- paste0(OutDir, '/Correlation_Peak_Spearman.txt')
con <- file(TextFile, "w") 

# pairwise processing of the input samples
for (i in (1:(NumSample-1))) {
	for (j in ((i+1):NumSample)) {

		XAxisData <- InpData[, 3+i]
		YAxisData <- InpData[, 3+j]
		AbsDiffVec <- abs(XAxisData - YAxisData)
		MinDiff <- min(AbsDiffVec)
		MaxDiff <- max(AbsDiffVec)
		ColorVal_CurrData <- ceiling(((AbsDiffVec - MinDiff) * length(ColorVec)) / ((MaxDiff - MinDiff) * 1.0))

		# plot the peak correlation
		plotfile1 <- paste0(OutDir, '/Correlation_Peak_', InpLabelList[i], '_', InpLabelList[j], '.pdf')	
		pdf(plotfile1, width=PlotWidth, height=PlotHeight)
		plot(XAxisData, YAxisData, cex=0.25, col=ColorVal_CurrData, xlab=paste0("Peak_Log10P_", InpLabelList[i]), ylab=paste0("Peak_Log10P_", InpLabelList[j]))
		title("Correlation between peak intensity")
		dev.off()

		# also print the correlation 
		Corr_val <- cor(XAxisData, YAxisData, method="spearman")
		outtext <- paste0("\n\n First peak file label : ", InpLabelList[i], "\n\n Second peak file label : ", InpLabelList[j], "\n\n Correlation value: ", Corr_val)
		writeLines(outtext, con=con, sep="\n")
	}
}

close(con)




