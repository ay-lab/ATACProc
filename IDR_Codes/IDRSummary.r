#!/usr/bin/env Rscript

#===========================================================
# R script for summarizing the results of IDR analysis between different sample replicates

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript result_summary.r $inpfile
#===========================================================

args <- commandArgs(TRUE)

# file containing the overlapped peak information
CommonPeakFile <- args[1]
inpdir <- dirname(CommonPeakFile)

npeak1 <- as.integer(args[2])
npeak2 <- as.integer(args[3])

# print(sprintf("\n CommonPeakFile: %s ", CommonPeakFile))
# print(sprintf("\n npeak1: %s ", npeak1))
# print(sprintf("\n npeak2: %s ", npeak2))

# information of the common peak
# Note: the file contains a header line
CommonPeakInfo <- read.table(CommonPeakFile, header=TRUE)

# number of overlapped peaks (considering all the IDR values)
ncommonpeak <- length(CommonPeakInfo[,1])
fracpeak1 <- (ncommonpeak * 1.0) / npeak1
fracpeak2 <- (ncommonpeak * 1.0) / npeak2

# print(sprintf("\n ncommonpeak: %s ", ncommonpeak))
# print(sprintf("\n fracpeak1: %s ", fracpeak1))
# print(sprintf("\n fracpeak2: %s ", fracpeak2))

# find the rows where IDR is lower than a specified threshold
# we employ three different thresholds:
# 1) 0.01, 2) 0.05, and 3) 0.1
# Note:  Threshold of 0.01 (newly added) is recommended in the ENCODE  

NumIDRPass0 <- length(which(CommonPeakInfo[,10] <= 0.01))
FracIDRPass0 <- (NumIDRPass0 * 1.0) / ncommonpeak
NumIDRPass1 <- length(which(CommonPeakInfo[,10] <= 0.05))
FracIDRPass1 <- (NumIDRPass1 * 1.0) / ncommonpeak
NumIDRPass2 <- length(which(CommonPeakInfo[,10] <= 0.1))
FracIDRPass2 <- (NumIDRPass2 * 1.0) / ncommonpeak

# print(sprintf("\n NumIDRPass0: %s ", NumIDRPass0))
# print(sprintf("\n FracIDRPass0: %s ", FracIDRPass0))
# print(sprintf("\n NumIDRPass1: %s ", NumIDRPass1))
# print(sprintf("\n FracIDRPass1: %s ", FracIDRPass1))
# print(sprintf("\n NumIDRPass2: %s ", NumIDRPass2))
# print(sprintf("\n FracIDRPass2: %s ", FracIDRPass2))

# divide the input overlapped peak files into two different structures
# corresponding to the peak information of two different inputs
# the seq() function also includes the row number for every interaction
# this row number serves as the id of peaks
PeakInfoInput1 <- cbind(seq(1:ncommonpeak), CommonPeakInfo[,1:4])
PeakInfoInput2 <- cbind(seq(1:ncommonpeak), CommonPeakInfo[,5:8])

# sort the data according to the significance value (last column of both the data)
# decreasing order is employed
PeakInfoInput1_Sort <- PeakInfoInput1[ order(-PeakInfoInput1[,5]),]
PeakInfoInput2_Sort <- PeakInfoInput2[ order(-PeakInfoInput2[,5]),]

# we check the cumulative percent of samples in both peak sets
# and find out the overlap of peaks
fraction_overlap <- c()

for (x in seq(0, 1, 0.1)) {
	if ((x != 0) && (x != 1)) {
		# number of elements of both peak lists
		nsample <- as.integer(ncommonpeak * x)
		# common elements in both peak lists
		# the common factor is the first column: peak id
		OverlapSet <- PeakInfoInput1_Sort[1:nsample, 1] %in% PeakInfoInput2_Sort[1:nsample, 1]
		ncommon <- length(OverlapSet[OverlapSet==TRUE])
		frac_common <- (ncommon * 1.0 / nsample)
		fraction_overlap <- c(fraction_overlap, frac_common)

		# we also note down two different fraction overlap statistics
		# corresponding to 10\%, 20% and 50% strongest peaks
		if (x == 0.1) {
			frac_overlap_10Pct = frac_common
		}
		if (x == 0.2) {
			frac_overlap_20Pct = frac_common
		}
		if (x == 0.5) {
			frac_overlap_50Pct = frac_common
		}

		# print(sprintf("\n Percentile value: %s ", x))
		# print(sprintf("\n nsample: %s ", nsample))
		# print(sprintf("\n ncommon: %s ", ncommon))
		# print(sprintf("\n frac_common: %s ", frac_common))
	}
}

# print(sprintf("\n Mean of fraction overlap: %s ", mean(fraction_overlap)))

# # we check the percent of samples in both peak sets
# # and find out the overlap of peaks
# # for individual 10% bins

# nbins <- 5	#10
# fraction_overlap2 <- c()
# sampleperbin <- as.integer(ncommonpeak / nbins)

# for (b in (1:nbins)) {
# 	if (b == 1) {
# 		si <- 1
# 		ei <- si + sampleperbin - 1
# 	} else {
# 		si <- ei + 1
# 		if (b == nbins) {
# 			ei <- ncommonpeak
# 		} else {
# 			ei <- si + sampleperbin - 1
# 		}
# 	}
# 	OverlapSet <- PeakInfoInput1_Sort[si:ei, 1] %in% PeakInfoInput2_Sort[si:ei, 1]
# 	ncommon <- length(OverlapSet[OverlapSet==TRUE])
# 	frac_common <- (ncommon * 1.0 / (ei-si+1))
# 	fraction_overlap2 <- c(fraction_overlap2, frac_common)
# 	print(sprintf("\n si: %s ", si))
# 	print(sprintf("\n ei: %s ", ei))
# 	print(sprintf("\n ncommon: %s ", ncommon))
# 	print(sprintf("\n frac_common: %s ", frac_common))
# }

# print(sprintf("\n Mean of fraction overlap2: %s ", mean(fraction_overlap2)))


# write the results in a text file
OutFilename <- paste0(inpdir, '/Stat.tab')

fp <- file(OutFilename, open="w")
write(paste0('NPeak1', '\t', 'NPeak2', '\t', 'CommonPeak', '\t', 'FracPeak1', '\t', 'FracPeak2', '\t', 'IDR_0.01_Peak', '\t', 'Frac_IDR_0.01_Peak', '\t', 'IDR_0.05_Peak', '\t', 'Frac_IDR_0.05_Peak', '\t', 'IDR_0.1_Peak', '\t', 'Frac_IDR_0.1_Peak', '\t', 'MeanOverlap', '\t', 'Overlap10', '\t', 'Overlap20', '\t', 'Overlap50'), file=fp, append=T)
write(paste(npeak1, '\t', npeak2, '\t', ncommonpeak, '\t', fracpeak1, '\t', fracpeak2, '\t', NumIDRPass0, '\t', FracIDRPass0, '\t', NumIDRPass1, '\t', FracIDRPass1, '\t', NumIDRPass2, '\t', FracIDRPass2, '\t', mean(fraction_overlap), '\t', frac_overlap_10Pct, '\t', frac_overlap_20Pct, '\t', frac_overlap_50Pct), file=fp, append=T)
close(fp)



