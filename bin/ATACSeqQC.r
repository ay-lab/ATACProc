#!/usr/bin/env Rscript

#==================================
# used for ATAC seq quality analysis
# adapted from the link:
# https://www.bioconductor.org/packages/3.7/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html

# following installation is to be done in R console:

# library(BiocInstaller)
# biocLite(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
#            "BSgenome.Hsapiens.UCSC.hg19", "TxDb.Hsapiens.UCSC.hg19.knownGene",
#            "phastCons100way.UCSC.hg19", 
#            "BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene", 
#            "phastCons60way.UCSC.mm10"))

# running command line arguments:
# 1) input bam file name
# 2) output directory where the results would be stored
# 3) boolean flag whether paired end (1) or single end (0) read is input
# 4) reference genome name - such as 'hg19' (default), 'mm9', 'mm10', etc.

#==================================
# author: Sourya Bhattacharyya
# Vijay-AY lab
#==================================

#============
# install package

# suppressPackageStartupMessages({
#   library(ATACseqQC)
#   library(ChIPpeakAnno)
#   library(BSgenome.Hsapiens.UCSC.hg19)
#   library(BSgenome.Hsapiens.UCSC.hg38)
#   library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#   library(phastCons100way.UCSC.hg19)
#   library(phastCons100way.UCSC.hg38)
#   library(MotifDb)
#   library(GenomicRanges)
#   library(GenomicAlignments)
# })

if(!require(GenomicAlignments)){
	library(BiocInstaller)
    biocLite(c("GenomicAlignments"))
}
library(GenomicAlignments)

if(!require(GenomicRanges)){
	library(BiocInstaller)
    biocLite(c("GenomicRanges"))
}
library(GenomicRanges)

if(!require(MotifDb)){
	library(BiocInstaller)
    biocLite(c("MotifDb"))
}
library(MotifDb)

if(!require(ChIPpeakAnno)){
	library(BiocInstaller)
    biocLite(c("ChIPpeakAnno"))
}
library(ChIPpeakAnno)

if(!require(ATACseqQC)){
	library(BiocInstaller)
    biocLite(c("ATACseqQC"))    
}
library(ATACseqQC)

if(!require(BSgenome.Hsapiens.UCSC.hg19)){
	library(BiocInstaller)
    biocLite(c("BSgenome.Hsapiens.UCSC.hg19"))
}
library(BSgenome.Hsapiens.UCSC.hg19)

if(!require(TxDb.Hsapiens.UCSC.hg19.knownGene)){
	library(BiocInstaller)
    biocLite(c("TxDb.Hsapiens.UCSC.hg19.knownGene"))
}
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

if(!require(phastCons100way.UCSC.hg19)){
	library(BiocInstaller)
    biocLite(c("phastCons100way.UCSC.hg19"))
}
library(phastCons100way.UCSC.hg19)

if(!require(BSgenome.Hsapiens.UCSC.hg38)){
	library(BiocInstaller)
    biocLite(c("BSgenome.Hsapiens.UCSC.hg38"))
}
library(BSgenome.Hsapiens.UCSC.hg38)

if(!require(TxDb.Hsapiens.UCSC.hg38.knownGene)){
	library(BiocInstaller)
    biocLite(c("TxDb.Hsapiens.UCSC.hg38.knownGene"))
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if(!require(phastCons100way.UCSC.hg38)){
	library(BiocInstaller)
    biocLite(c("phastCons100way.UCSC.hg38"))
}
library(phastCons100way.UCSC.hg38)


# if(!require(BSgenome.Mmusculus.UCSC.mm9)){
# 	library(BiocInstaller)
#     biocLite(c("BSgenome.Mmusculus.UCSC.mm9"))
#     library(BSgenome.Mmusculus.UCSC.mm9)
# }

# if(!require(TxDb.Mmusculus.UCSC.mm9.knownGene)){
# 	library(BiocInstaller)
#     biocLite(c("TxDb.Mmusculus.UCSC.mm9.knownGene"))
#     library(TxDb.Mmusculus.UCSC.mm9.knownGene)
# }

# if(!require(phastCons60way.UCSC.mm9)){
# 	library(BiocInstaller)
#     biocLite(c("phastCons60way.UCSC.mm9"))
#     library(phastCons60way.UCSC.mm9)
# }

# if(!require(BSgenome.Mmusculus.UCSC.mm10)){
# 	library(BiocInstaller)
#     biocLite(c("BSgenome.Mmusculus.UCSC.mm10"))
#     library(BSgenome.Mmusculus.UCSC.mm10)
# }

# if(!require(TxDb.Mmusculus.UCSC.mm10.knownGene)){
# 	library(BiocInstaller)
#     biocLite(c("TxDb.Mmusculus.UCSC.mm10.knownGene"))
#     library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# }

# if(!require(phastCons60way.UCSC.mm10)){
# 	library(BiocInstaller)
#     biocLite(c("phastCons60way.UCSC.mm10"))
#     library(phastCons60way.UCSC.mm10)
# }

# end install package
#============

knitr::opts_chunk$set(warning=FALSE, message=FALSE)

args <- commandArgs(TRUE)

cat(sprintf("length args: ", length(args)))

if (length(args) == 0) {
	cat(sprintf("No input argument is provided -- return "))
	return()
}

# input bam file name
# convert input file to the absolute path
inpbamfilename <- normalizePath(args[1])

# output directory containing results
if (length(args) > 1) {
	outdir <- args[2]
} else {
	outdir <- dirname(inpbamfilename)
}
# convert to absolute path
outdir <- normalizePath(outdir)
system(paste("mkdir -p", outdir))

# option (boolean) if the input bam file is of paired end 
# default 0
if (length(args) > 2) {
	paired_end_input <- as.integer(args[3])
} else {
	paired_end_input <- 0
}

# reference genome 
if (length(args) > 3) {
	refgenome <- args[4]
} else {
	refgenome <- 'hg19'	
}

cat(sprintf("\n inpbamfilename: %s ", inpbamfilename))
cat(sprintf("\n outdir: %s ", outdir))
cat(sprintf("\n refgenome: %s ", refgenome))
cat(sprintf("\n paired_end_input: %s ", paired_end_input))

# store the current directory
currdir <- getwd()

# go to the directory "outdir"
setwd(outdir)

# write the QC measures in a text file
bam_QC_textfile <- paste0(outdir, '/Summary_QC.txt')

# input the bamFile from the ATACseqQC package 
# bamfile <- system.file("extdata", inpbamfilename, package="ATACseqQC", mustWork=TRUE)
bamfile <- inpbamfilename
bamfile.labels <- gsub(".bam", "", basename(bamfile))

# following code is commented for the moment - sourya
# it is working fine
if (1) {
	# Check alignment metrics and mapping quality
	fp_out <- file(bam_QC_textfile, "w")
	outtext <- paste0("\n *** Quality control measures corresponding to the alignment file ", inpbamfilename, " is  **** \n\n ")
	writeLines(outtext, con=fp_out, sep="\n")
	close(fp_out)

	# quality control summary
	capture.output(bamQC(bamfile, outPath=NULL), file=bam_QC_textfile, append=TRUE)

	cat(sprintf("\n\n *** Computed the sumary statistics of the input ATAC seq file ****  \n\n"))
}

# following code is commented for the moment - sourya
# the function is not found
if (0) {
	# Estimate the library complexity
	estimateLibComplexity(readsDupFreq(bamfile))
	# # plot the library complexity
	# plotfile <- paste0(outdir, '/Library_Complexity.pdf')
	# pdf(plotfile, width=10, height=7)
	# plot(lib_compl)
	# title("Library complexity - ATAC seq")
	# dev.off()	
}

#==================================
# following codes work only when the input is paired end
#==================================
if (paired_end_input == 1) {
	
	# generate fragement size distribution
	# works only when the input is paired end
	fragSize <- fragSizeDist(bamfile, bamfile.labels)
	cat(sprintf("\n\n *** Computed the fragment size distribution (paired end read) of the input ATAC seq file ****  \n\n"))

	#==========================
	# Nucleosome positioning
	# Adjust the read start sites
	#==========================
	# bamfile tags to be read in
	tags <- c('AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YS', 'YT')

	## files will be output into outPath
	outPath <- paste0(outdir, '/splited')
	system(paste("mkdir -p", outPath))

	## shift the coordinates of 5'ends of alignments in the bam file

	# comment - sourya - only for chromosome 1
	# seqlev <- "chr1" ## subsample data for quick run
	# which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")

	# considering all the chromosomes
	seqlev <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")
	which <- as(seqinfo(Hsapiens)[seqlev], "GRanges")
	
	gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE)

	cat(sprintf("\n\n *** Performed function - readBamFile ****  \n\n"))

	gal1 <- shiftGAlignmentsList(gal)
	shiftedBamfile <- file.path(outPath, "shifted.bam")
	export(gal1, shiftedBamfile)

	cat(sprintf("\n\n *** Exported the shifted bam file for the input ATAC seq file ****  \n\n"))

	## --------------------------------------------------------------------------
	txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

	# # comment - sourya
	# ## run program for chromosome 1 only
	# txs <- txs[seqnames(txs) %in% "chr1"]

	genome <- Hsapiens
	## split the reads into NucleosomeFree, mononucleosome, 
	## dinucleosome and trinucleosome.
	objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, conservation=phastCons100way.UCSC.hg19)

	cat(sprintf("\n\n *** Performed function - splitGAlignmentsByCut ****  \n\n"))

	## --------------------------------------------------------------------------
	null <- writeListOfGAlignments(objs, outPath)
	## list the files generated by splitBam.
	dir(outPath)

	cat(sprintf("\n\n *** Performed function - writeListOfGAlignments ****  \n\n"))

	## ----eval=FALSE------------------------------------------------------------
	#  objs <- splitBam(bamfile, tags=tags, outPath=outPath,
	#                   txs=txs, genome=genome,
	#                   conservation=phastCons100way.UCSC.hg19)

	## ----fig.height=4, fig.width=4---------------------------------------------
	bamfiles <- file.path(outPath, c("NucleosomeFree.bam", "mononucleosome.bam", "dinucleosome.bam", "trinucleosome.bam"))

	## Plot the cumulative percentage of tag allocation in nucleosome-free 
	## and mononucleosome bam files.

	# modified - sourya
	# earlier - using only chromosome 1
	# cumulativePercentage(bamfiles[1:2], as(seqinfo(Hsapiens)["chr1"], "GRanges"))
	# now using all the candidate chromosomes
	cumulativePercentage(bamfiles[1:2], as(seqinfo(Hsapiens)[seqlev], "GRanges"))

	cat(sprintf("\n\n *** Performed function - cumulativePercentage ****  \n\n"))

	## ----fig.height=8, fig.width=4---------------------------------------------
	TSS <- promoters(txs, upstream=0, downstream=1)
	TSS <- unique(TSS)
	## estimate the library size for normalization
	(librarySize <- estLibSize(bamfiles))
	## calculate the signals around TSSs.
	NTILE <- 101
	dws <- ups <- 1010
	sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")], TSS=TSS, librarySize=librarySize, seqlev=seqlev, TSS.filter=0.5, n.tile = NTILE, upstream = ups, downstream = dws)
	## log2 transformed signals
	sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
	#plot heatmap
	featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws), zeroAt=.5, n.tile=NTILE)

	cat(sprintf("\n\n *** Performed function - featureAlignedHeatmap ****  \n\n"))

	## ----fig.show="hide"-------------------------------------------------------
	## get signals normalized for nucleosome-free and nucleosome-bound regions.
	out <- featureAlignedDistribution(sigs, reCenterPeaks(TSS, width=ups+dws), zeroAt=.5, n.tile=NTILE, type="l", ylab="Averaged coverage")

	cat(sprintf("\n\n *** Performed function - featureAlignedDistribution ****  \n\n"))

	## --------------------------------------------------------------------------
	## rescale the nucleosome-free and nucleosome signals to 0~1
	range01 <- function(x){(x-min(x))/(max(x)-min(x))}
	out <- apply(out, 2, range01)
	matplot(out, type="l", xaxt="n", xlab="Position (bp)", ylab="Fraction of signal")
	axis(1, at=seq(0, 100, by=10)+1, labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
	abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")

	## --------------------------------------------------------------------------
	## foot prints
	CTCF <- query(MotifDb, c("CTCF"))
	CTCF <- as.list(CTCF)
	print(CTCF[[1]], digits=2)
	sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]], genome=genome, min.score="90%", seqlev=seqlev, upstream=100, downstream=100)

	## ----fig.height=6, fig.width=6---------------------------------------------
	featureAlignedHeatmap(sigs$signal, feature.gr=reCenterPeaks(sigs$bindingSites, width=200+width(sigs$bindingSites[1])), annoMcols="score", sortBy="score", n.tile=ncol(sigs$signal[[1]]))

	cat(sprintf("\n\n *** Performed function - featureAlignedHeatmap - after motif analysis ****  \n\n"))

	spearman_corr <- sigs$spearman.correlation

	fp_out <- file(bam_QC_textfile, "a")
	outtext <- paste0("\n\n *** Spearman correlation of the signal:  ", spearman_corr, "  **** \n\n ")
	writeLines(outtext, con=fp_out, sep="\n")
	close(fp_out)

	# rename the file (plots) 
	system(paste("mv Rplots.pdf All_Plots.pdf"))

} 	# end check paired end read

# return to the original directory
setwd(currdir)

