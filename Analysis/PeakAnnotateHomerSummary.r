#!/usr/bin/env Rscript

#=============================
# this code plots the percentage of genomic annotations
# for each peak file generated from the MACS2 command
# of a ChIP-seq pipeline
# input is a HOMER annotation of the peak file
# corresponding to different genomic segments
#=============================

args <- commandArgs(TRUE)
if(length(args)<1)	{
	q("no")
}

HomerAnnotFile <- args[1]
HomerPeakTSSDistFile <- args[2]

InpDir <- dirname(HomerAnnotFile)

#====================================
# first process the peak annotations
# produced by HOMER
#====================================

# first initialize different annotation categories
npeak_3UTR <- 0
npeak_TTS <- 0
npeak_Exon <- 0
npeak_Intron <- 0
npeak_Intergenic <- 0
npeak_Promoter <- 0
npeak_5UTR <- 0

# open the file and read line by line
# extract performance values
finp <- file(HomerAnnotFile, "r")
lineset <- readLines(finp)
for (i in 1:length(lineset)) {
	curr_line <- trimws(lineset[i], which = "both")
	curr_line_split <- strsplit(curr_line,"\t")[[1]]
	# cat(sprintf("\n curr_line : %s ", curr_line_split))
	# cat(sprintf("\n elem 1: %s elem 2: %s elem 3: %s elem 4: %s ", curr_line_split[1], curr_line_split[2], curr_line_split[3], curr_line_split[4]))
	if (regexpr("3UTR", curr_line) > 0) {
		npeak_3UTR <- as.numeric(curr_line_split[2])
	} else if (regexpr("TTS", curr_line) > 0) {
		npeak_TTS <- as.numeric(curr_line_split[2])
	} else if (regexpr("Exon", curr_line) > 0) {
		npeak_Exon <- as.numeric(curr_line_split[2])
	} else if (regexpr("Intron", curr_line) > 0) {
		npeak_Intron <- as.numeric(curr_line_split[2])
	} else if (regexpr("Intergenic", curr_line) > 0) {
		npeak_Intergenic <- as.numeric(curr_line_split[2])
	} else if (regexpr("Promoter", curr_line) > 0) {
		npeak_Promoter <- as.numeric(curr_line_split[2])
	} else if (regexpr("5UTR", curr_line) > 0) {
		npeak_5UTR <- as.numeric(curr_line_split[2])
	} 
}

# close the input file
close(finp)

cat(sprintf("\n npeak_3UTR : %s ", npeak_3UTR))
cat(sprintf("\n npeak_TTS : %s ", npeak_TTS))
cat(sprintf("\n npeak_Exon : %s ", npeak_Exon))
cat(sprintf("\n npeak_Intron : %s ", npeak_Intron))
cat(sprintf("\n npeak_Intergenic : %s ", npeak_Intergenic))
cat(sprintf("\n npeak_Promoter : %s ", npeak_Promoter))
cat(sprintf("\n npeak_5UTR : %s ", npeak_5UTR))

# now prepare a vector of the above categories (provided non zero instances)
# to create a pie chart
slices <- c()
lbls <- c()
if (npeak_3UTR > 0) {
	slices <- c(slices, npeak_3UTR)
	lbls <- c(lbls, "3UTR")
}
if (npeak_TTS > 0) {
	slices <- c(slices, npeak_TTS)
	lbls <- c(lbls, "TTS")
}
if (npeak_Exon > 0) {
	slices <- c(slices, npeak_Exon)
	lbls <- c(lbls, "Exon")
}
if (npeak_Intron > 0) {
	slices <- c(slices, npeak_Intron)
	lbls <- c(lbls, "Intron")
}
if (npeak_Intergenic > 0) {
	slices <- c(slices, npeak_Intergenic)
	lbls <- c(lbls, "Intergenic")
}
if (npeak_Promoter > 0) {
	slices <- c(slices, npeak_Promoter)
	lbls <- c(lbls, "Promoter")
}
if (npeak_5UTR > 0) {
	slices <- c(slices, npeak_5UTR)
	lbls <- c(lbls, "5UTR")
}

# convert the vector to include the percentage values as well
# for displaying in the pie chart
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 

OutPlotFile <- paste0(InpDir, "/Pie_Chart_Peak_Annotation.pdf")
pdf(OutPlotFile, width=6, height=4)
pie(slices, labels = lbls, col=rainbow(length(lbls)), main="Pie Chart of peak annotation", radius = 1, cex = 0.5)
dev.off()

#====================================
# then process the distance from TSS sites (nearest)
# for individual peaks
# the histogram data is already provided
#====================================
# first remove the header line from the input file
tempfile <- paste0(InpDir, '/temp_TSSDistFile.bed')
system(paste("awk \'NR>1\'", HomerPeakTSSDistFile, ">", tempfile))

# now process the temporary file
PeakTSSData <- read.table(tempfile, header=T)
OutPlotFile <- paste0(InpDir, "/Peak_TSS_Distance.pdf")
pdf(OutPlotFile, width=6, height=4)
plot(PeakTSSData[,1], PeakTSSData[,2], cex=0.5, col="red", xlab="Distance from TSS", ylab="ChIP fragment depth (per bp per peak)")
title("Peak distribution near TSS sites")
dev.off()

# then remove the temporary file
system(paste("rm", tempfile))

