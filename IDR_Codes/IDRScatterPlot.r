#!/usr/bin/env Rscript

#===========================================================
# R script for scatter plot between a pair of peak files

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript result_summary.r $inpfile
#===========================================================

args <- commandArgs(TRUE)

# directory containing IDR code package
IDRCodeDir <- args[1]	#"/home/sourya/packages/idrCode/"

# the pair of peak outputs for comparison
peakfile1 <- args[2]
peakfile2 <- args[3]
# genome table.txt file provided in the IDR code package
genometablefile <- args[4]
# output prefix including the directory path
# and the prefix string of the output plot file names
curroutprefix <- args[5]

# system path includes the path of IDR code 
source(paste0(IDRCodeDir, "functions-all-clayton-12-13.r"))

chr.size <- read.table(genometablefile)

half.width <- NULL
overlap.ratio <- 0
is.broadpeak <- F
sig.value <- "p.value"

# width and height values employed in these plots
plotwidth <- 8
plotheight <- 6

rep1 <- process.narrowpeak(paste(peakfile1, sep=""), chr.size, 
	half.width=half.width, summit="offset", broadpeak=is.broadpeak)

rep2 <- process.narrowpeak(paste(peakfile2, sep=""), chr.size, 
	half.width=half.width, summit="offset", broadpeak=is.broadpeak)

uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, 
	sig.value1=sig.value, sig.value2=sig.value, overlap.ratio=overlap.ratio)

em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)
idr.local <- 1-em.output$em.fit$e.z
IDR <- c()
o <- order(idr.local)
IDR[o] <- cumsum(idr.local[o])/c(1:length(o))

idr_output <- data.frame(chr1=em.output$data.pruned$sample1[, "chr"], start1=em.output$data.pruned$sample1[, "start.ori"], stop1=em.output$data.pruned$sample1[, "stop.ori"], sig.value1=em.output$data.pruned$sample1[, "sig.value"],  chr2=em.output$data.pruned$sample2[, "chr"], start2=em.output$data.pruned$sample2[, "start.ori"], stop2=em.output$data.pruned$sample2[, "stop.ori"], sig.value2=em.output$data.pruned$sample2[, "sig.value"], idr.local=1-em.output$em.fit$e.z, IDR=IDR)

# this idr_output is already placed in the file "idr_overlapped_peaks.txt" 

filtered_peaks <- idr_output[idr_output[,10]<=0.01,]
dim(filtered_peaks) # get the number of peaks

ez.list <- get.ez.tt.all(em.output, uri.output$data12.enrich$merge1, uri.output$data12.enrich$merge2)

par(mar=c(5,5,0,0.5), mfrow = c(1,3), oma=c(5,0,2,0))

idr_output$col[idr_output[,10]<=0.01]="black"

idr_output$col[idr_output[,10]>=0.01]="red"

# first graph
pdf(paste0(curroutprefix,'_Signal_Replicates.pdf'), width=plotwidth, height=plotheight)
plot(log(idr_output[,4]),log(idr_output[,8]),col=idr_output[,11], pch=19, cex = 0.05, xlab="log(signal) Rep1", ylab="log(signal) Rep2")
legend("topleft", c("IDR=>0.01","IDR<=0.01"), col=c("red","black"), pch=19, bty="n", lty=c(1,1), lwd=c(2,2))
dev.off()

# second graph
pdf(paste0(curroutprefix,'_Peak_Rank_Replicates.pdf'), width=plotwidth, height=plotheight)
plot(rank(-idr_output[,4]),rank(-idr_output[,8]),col=idr_output[,11], pch=19, cex = 0.05, xlab="Peak rank Rep1", ylab="Peak rank Rep2")
legend("topleft", c("IDR=>0.01","IDR<=0.01"), col=c("red","black"), pch=19, bty="n", lty=c(1,1), lwd=c(1,1))
dev.off()

# third graph
pdf(paste0(curroutprefix,'_SignificantPeaks_vs_IDR.pdf'), width=plotwidth, height=plotheight)
plot(ez.list$IDR, ylab="IDR", xlab="num of significant peaks")
dev.off()





