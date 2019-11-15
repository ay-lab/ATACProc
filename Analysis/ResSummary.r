#!/usr/bin/env Rscript

#==================================
# used to print the summary result statistics 
# from a collection of ATAC seq samples
# this script is to be called on the top most directory structure 
# containing all the ATAC seq sample folders 

# author: Sourya Bhattacharyya
# Vijay-AY lab
#==================================

library(optparse)
library(ggplot2)
library(plotly)

#==============
# function to plot scatter using plotly package
#==============
PlotScatter_Data <- function(InpDF, ylabel, plotfile) {
	colnames(InpDF) <- c('X', 'Y')
	currplot <- plotly::plot_ly(InpDF, x= ~X, y= ~Y, name=ylabel, type="scatter", mode="markers", marker=list(size=10, color= "blue")) %>% layout(xaxis = list(title = 'Samples', zeroline = FALSE, showticklabels = FALSE), yaxis = list(title = ylabel, zeroline = FALSE))
	htmlwidgets::saveWidget(currplot, plotfile)
}

#===========================================================
option_list = list(
	make_option(c("--BaseDir"), type="character", default=NULL, help="Base directory containing all the ATAC-seq samples. Mandatory parameter."),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory to contain the summary results. If empty, current directory is used.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$BaseDir)) {
	print_help(opt_parser)
	stop("ERROR !!!!!!! Base output directory is not provided - check the option --BaseDir \n", call.=FALSE)
}

if (is.null(opt$OutDir)) {
	OutDir <- getwd()
} else {
	OutDir <- opt$OutDir
	system(paste("mkdir -p", OutDir))
}

# check if the last character of the output base directory is '/'
# unless append that character
baseresdir <- opt$BaseDir
if (substr(baseresdir,nchar(baseresdir),nchar(baseresdir)) != '/') {
	baseresdir <- paste0(baseresdir, '/')
}

# template name of the directories containing MACS2 results
# either default parameters are used
# or extsize based parameters are employed
MACS2_def_dir <- 'MACS2_Default_Tag'
MACS2_ext_dir <- 'MACS2_Ext_Tag'
# MACS2_noduprem_ext_dir <- 'MACS2_NoDupRem_Align_Ext_Tag'

# template name of the folders (peak outputs) depending on whether control (input) are used for peak detection
Ctrl_0_Fold <- '_No_Control'
Ctrl_1_Fold <- '_with_Control'

# file formats of NRF, read statistics, FRiP and Peak statistics
# which are present for every sample
ReadStatFileNameFmt <- 'Read_Count_Stat.txt'
NRFfilenamefmt <- 'out_NRF'
FRiPFileNameFmt <- 'out_FRiP.txt'
PeakCountFileFmt <- 'Peak_Statistics.txt'

# output text file to contain the summary results
outtextfile <- paste0(baseresdir, 'Results_All_Samples_Summary.txt')

#=================================
file_process <- FALSE

# process individual directories under the main results directory
dir.list <- list.dirs(path = baseresdir, full.names = FALSE, recursive = FALSE)
for (dr in dir.list) {

	# cat(sprintf("\n Examining directory: %s \n", dr))	

	#==================
	# following file stores the count of reads throughout various stages of filtering
	#==================
	ReadStatFile <- paste0(baseresdir, dr, "/", ReadStatFileNameFmt)
	if (file.exists(ReadStatFile) && (file.access(ReadStatFile, 4) == 0)) {
		cat(sprintf("\n Found the file: %s \n", ReadStatFile))

		# search for a file with name *_R1*.fastq.gz in the current file
		filenames <- Sys.glob(paste0(baseresdir, dr, "/*_R1*.fastq.gz"))
		currSampleName <- substr(basename(filenames[1]), start=1, stop=regexpr("_R1", basename(filenames[1]))-1)

		x <- readLines(ReadStatFile)
		lastline <- strsplit(x[length(x)], "\t")[[1]]

		# the line should have 8 or 7 fields
		# 8 fields if fastq file is used in the pipeline
		# 7 fields if already aligned file is used in the pipeline
		nfields <- length(lastline)
		cat(sprintf("\n No of fields in the read statistics file: %s ", nfields))

	 	TotRead <- as.integer(lastline[nfields-6])	#lastline[2]
	 	MappableRead <- as.integer(lastline[nfields-5])	#lastline[3]
	 	Frac_Mappable_Read <- ((MappableRead * 1.0) / TotRead)
	 	Frac_Unmappable_Read <- (((TotRead - MappableRead) * 1.0) / TotRead)
	 	Read_remain_after_RandomDel <- as.integer(lastline[nfields-4])	#lastline[4]
	 	Frac_reads_remain_after_RandomDel <- ((Read_remain_after_RandomDel * 1.0) / TotRead)
	 	Frac_reads_deleted_random <- (((MappableRead - Read_remain_after_RandomDel) * 1.0) / TotRead)
	 	Read_remain_after_Mitochondrial_Read_Del <- as.integer(lastline[nfields-3])	#lastline[5]
	 	Frac_reads_remain_after_MtReadDel <- ((Read_remain_after_Mitochondrial_Read_Del * 1.0) / TotRead)
	 	Frac_reads_deleted_MtRead <- (((Read_remain_after_RandomDel - Read_remain_after_Mitochondrial_Read_Del) * 1.0) / TotRead)
	 	UniqMappedRead <- as.integer(lastline[nfields-2])	#lastline[6]
	 	Frac_reads_unique_mapped <- ((UniqMappedRead * 1.0) / TotRead)
	 	Frac_reads_del_multimap <- (((Read_remain_after_Mitochondrial_Read_Del - UniqMappedRead) * 1.0) / TotRead)
	 	ReadQualThr <- as.integer(lastline[nfields-1])	#lastline[7]
 		Frac_reads_remain_QualThr <- ((ReadQualThr * 1.0) / TotRead)
 		Frac_reads_del_QualThr <- (((UniqMappedRead - ReadQualThr) * 1.0) / TotRead)
		Dupl_Rem_Read <- lastline[nfields]	#lastline[9]
		Frac_reads_remain_Dupl <- ((Dupl_Rem_Read * 1.0) / TotRead)
		Frac_reads_del_Dupl <- (((ReadQualThr - Dupl_Rem_Read) * 1.0) / TotRead)

		# append the entries in the final vector
		CurrOutVec <- c(basename(dr), currSampleName, TotRead, MappableRead, Frac_Mappable_Read, Frac_Unmappable_Read, Read_remain_after_RandomDel, Frac_reads_remain_after_RandomDel, Frac_reads_deleted_random, Read_remain_after_Mitochondrial_Read_Del, Frac_reads_remain_after_MtReadDel, Frac_reads_deleted_MtRead, UniqMappedRead, Frac_reads_unique_mapped, Frac_reads_del_multimap, ReadQualThr, Frac_reads_remain_QualThr, Frac_reads_del_QualThr, Dupl_Rem_Read, Frac_reads_remain_Dupl, Frac_reads_del_Dupl)

		#==================
		# the following file stores the NRF / library complexity value
		#==================
		filenames <- Sys.glob(paste0(baseresdir, dr, "/*", NRFfilenamefmt, "*.txt"))
		if (length(filenames) > 0) {
			NRF_textfile <- filenames[1]
			if (file.exists(NRF_textfile) && (file.access(NRF_textfile, 4) == 0)){
			 	x <- readLines(NRF_textfile)
			 	# the 2nd line in string splitted structure
			 	lastline <- strsplit(x[length(x)], "\t")[[1]]
			 	UniqMappedPos <- lastline[2]
			  	NRF_val <- lastline[3]
			  	M1 <- lastline[4]
			  	M2 <- lastline[5]
			  	PBC1 <- lastline[6]
			  	PBC2 <- lastline[7]
			  	# adjust the output vector
			  	CurrOutVec <- c(CurrOutVec, UniqMappedPos, NRF_val, M1, M2, PBC1, PBC2)
		 	} else {
		 		CurrOutVec <- c(CurrOutVec, rep('NA', 6))
		 	}
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 6))
	 	}

 		#==================
	 	# check the peak directories and find corresponding statistics
	 	#==================

		#==================
		# FRiP and Peak count measures - default peak calling - no control
		#==================	
		FRiP_textfile_def_noctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_0_Fold, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_def_noctrl) && (file.access(FRiP_textfile_def_noctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_def_noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeak_def_noctrl <- lastline[2]
		  	FRiP_def_noctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, MappedReadPeak_def_noctrl, FRiP_def_noctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 2))
		}
	
		PeakCount_TextFile_def_noctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_0_Fold, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_def_noctrl) && (file.access(PeakCount_TextFile_def_noctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_def_noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeak_def_noctrl <- lastline[1]
		  	TotPeak_Q_Five_Pct_def_noctrl <- lastline[2]
		  	TotPeak_Q_One_Pct_def_noctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, TotPeak_def_noctrl, TotPeak_Q_Five_Pct_def_noctrl, TotPeak_Q_One_Pct_def_noctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 3))
		}

		#==================
		# FRiP and Peak count measures - Ext peak calling - no control
		#==================	
		FRiP_textfile_ext_noctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_0_Fold, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_ext_noctrl) && (file.access(FRiP_textfile_ext_noctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_ext_noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeak_ext_noctrl <- lastline[2]
		  	FRiP_ext_noctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, MappedReadPeak_ext_noctrl, FRiP_ext_noctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 2))
		}

		PeakCount_TextFile_ext_noctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_0_Fold, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_ext_noctrl) && (file.access(PeakCount_TextFile_ext_noctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_ext_noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeak_ext_noctrl <- lastline[1]
		  	TotPeak_Q_Five_Pct_ext_noctrl <- lastline[2]
		  	TotPeak_Q_One_Pct_ext_noctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, TotPeak_ext_noctrl, TotPeak_Q_Five_Pct_ext_noctrl, TotPeak_Q_One_Pct_ext_noctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 3))
		}

		#==================
		# FRiP and Peak count measures - default peak calling - with control
		#==================	
		FRiP_textfile_def_ctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_1_Fold, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_def_ctrl) && (file.access(FRiP_textfile_def_ctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_def_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeak_def_ctrl <- lastline[2]
		  	FRiP_def_ctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, MappedReadPeak_def_ctrl, FRiP_def_ctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 2))
		}

		PeakCount_TextFile_def_ctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_1_Fold, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_def_ctrl) && (file.access(PeakCount_TextFile_def_ctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_def_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeak_def_ctrl <- lastline[1]
		  	TotPeak_Q_Five_Pct_def_ctrl <- lastline[2]
		  	TotPeak_Q_One_Pct_def_ctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, TotPeak_def_ctrl, TotPeak_Q_Five_Pct_def_ctrl, TotPeak_Q_One_Pct_def_ctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 3))
		}

		#==================
		# FRiP and Peak count measures - Ext peak calling - with control
		#==================	
		FRiP_textfile_ext_ctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_1_Fold, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_ext_ctrl) && (file.access(FRiP_textfile_ext_ctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_ext_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeak_ext_ctrl <- lastline[2]
		  	FRiP_ext_ctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, MappedReadPeak_ext_ctrl, FRiP_ext_ctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 2))
		}

		PeakCount_TextFile_ext_ctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_1_Fold, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_ext_ctrl) && (file.access(PeakCount_TextFile_ext_ctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_ext_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeak_ext_ctrl <- lastline[1]
		  	TotPeak_Q_Five_Pct_ext_ctrl <- lastline[2]
		  	TotPeak_Q_One_Pct_ext_ctrl <- lastline[3]
		  	# adjust the output vector
		  	CurrOutVec <- c(CurrOutVec, TotPeak_ext_ctrl, TotPeak_Q_Five_Pct_ext_ctrl, TotPeak_Q_One_Pct_ext_ctrl)
		} else {
			CurrOutVec <- c(CurrOutVec, rep('NA', 3))
		}

		# now convert the current vector in a data frame
		CurrDF <- data.frame(CurrOutVec, nrow=1, ncol=length(CurrOutVec))
		colnames(CurrDF) <- c('Dir', 'SampleName', 'Total_Read', 'Number_of_Mappable_Reads', 'Fraction_of_Mappable_Reads', 'Fraction_of_Unmappable_Reads', 'Reads_after_random_chromsome_deletion', 'Fraction_reads_remain_after_random_chromsome_deletion', 'Fraction_reads_in_random_chromsome', 'Reads_excluding_Mitochondrial_Reads', 'Fraction_reads_excluding_MtRead', 'Fraction_mitochondrial_reads', 'UniqMappedRead', 'Fraction_reads_unique_mapped', 'Fraction_reads_multimap', 'Reads_remain_after_QualThr', 'Fraction_reads_remain_QualThr', 'Frac_reads_low_qual', 'Reads_after_dupl_remove', 'Fraction_de-duplicated_reads', 'Fraction_duplicate_reads', 'UniqMappedPos', 'NRF', 'M1', 'M2', 'PBC1', 'PBC2', 'MappedReadPeak_Def_noctrl(Q<0.05)', 'FRiP_Def_NoCtrl(Q<0.05)', 'nPeak_Def_NoCtrl', 'nPeak_Def_NoCtrl(Q<0.05)', 'nPeak_Def_NoCtrl(Q<0.01)', 'MapReadPeak_Ext_NoCtrl(Q<0.05)', 'FRiP_Ext_NoCtrl(Q<0.05)', 'nPeak_Ext_NoCtrl', 'nPeak_Ext_NoCtrl(Q<0.05)', 'nPeak_Ext_NoCtrl(Q<0.01)', 'MapReadPeak_Def_Ctrl(Q<0.05)',  'FRiP_Def_Ctrl(Q<0.05)', 'nPeak_Def_Ctrl', 'nPeak_Def_Ctrl(Q<0.05)', 'nPeak_Def_Ctrl(Q<0.01)', 'MapReadPeak_Ext_Ctrl(Q<0.05)', 'FRiP_Ext_Ctrl(Q<0.05)', 'nPeak_Ext_Ctrl', 'nPeak_Ext_Ctrl(Q<0.05)', 'nPeak_Ext_Ctrl(Q<0.01)')

	 	if (file_process == FALSE) {
	 		FinalDF <- CurrDF
	 		file_process <- TRUE
	 	} else {
	 		FinalDF <- rbind.data.frame(FinalDF, CurrDF)
	 	}
	 	
	}	# end processing current sample condition 

}	# end directory traverse

# now remove one or more columns of this data frame, if they are all 'NA'
NA_ColList <- c()
for (i in (1:ncol(FinalDF))) {
	idx <- which(FinalDF[, i] == 'NA')
	if (length(idx) == nrow(FinalDF)) {
		# every entry of this column is NA. So discard this column
		NA_ColList <- c(NA_ColList, i)
	}
}
# if (length(NA_ColList) > 0) {
# 	FinalDF_Modified <- FinalDF[-c(NA_ColList)]
# 	cat(sprintf("\n *** Dropped one or more columns since all entries were NA - before dropping : number of columns : %s after dropping columns : %s ", ncol(FinalDF), ncol(FinalDF_Modified)))
# 	write.table(FinalDF_Modified, outtextfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
# } else {
# 	write.table(FinalDF, outtextfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
# }

#============================
# a few summary statements for the output text file
#============================

CommentsFile <- paste0(OutDir, '/Field_Description.txt')

# open the output text file
con <- file(CommentsFile, "a")

outtext <- paste0("\n\n\n\n\n *** Important parameters ***** \n\n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Total_Read: number of reads in individual fastq file(s)")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Number_of_Mappable_Reads, Fraction_of_Mappable_Reads, and Fraction_of_Unmappable_Reads: number (and fraction) of reads mappable and unmappable to the reference genome. May not be uniquely mappable reads.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Reads_after_random_chromsome_deletion, Fraction_reads_remain_after_random_chromsome_deletion, and Fraction_reads_in_random_chromsome: number (and fraction) of reads remaining (and deleted) after deleting reads from random chromosomes such as chr1_*, chr2_*, chrUN, ....")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Reads_excluding_Mitochondrial_Reads, Fraction_reads_excluding_MtRead and Fraction_mitochondrial_reads: number (and fraction) of reads remaining (and deleted) after removing the mitochondrial reads")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  UniqMappedRead, Fraction_reads_unique_mapped, and Fraction_reads_multimap: number (and fraction) of reads uniquely mapped (and multimapped) to the reference genome.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Reads_remain_after_QualThr, Fraction_reads_remain_QualThr and Frac_reads_low_qual: number (and fraction) of reads remaining (and deleted) after removing low quality reads (MAPQ threshold)")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Reads_not_in_blackList_genome, Fraction_reads_not_in_blackList_genome, and Fraction_reads_in_blackList_genome: number (and fraction) of reads not in (and in) blacklist segments.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  Reads_after_dupl_remove, Fraction_de-duplicated_reads and Fraction_duplicate_reads: number (and fraction) of reads remaining (and deleted) after removing duplicate reads")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  UniqMapPos: number of distinct genome position where at least one read maps uniquely.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  NRF (Non redundant fraction): number of distinct genome positions for uniquely mapped reads / number of uniquely mapped reads")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  M1: number of genomic locations where exactly one read maps uniquely.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  M2: number of genomic locations where exactly two reads map uniquely.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  PBC1: M1 / UniqMapPos ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  PBC2: M1 / M2")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n\n\n\n  MACS2 outputs corresponding to peaks with default MACS2 command and no control ----- input  missing values are replaced by NA \n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n MappedReadPeak_Def_noctrl: mapped reads in peaks ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n FRiP_Def_NoCtrl: MappedReadPeak_Def_noctrl / UniqMappedRead ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_NoCtrl: number of peaks (determined by p value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_NoCtrl(Q<0.05): number of peaks (determined by q value threshold of 0.05) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_NoCtrl(Q<0.01): number of peaks (determined by q value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n\n\n\n  MACS2 outputs corresponding to peaks with --extsize option (recommended in existing ATAC seq pipeline) and no control input ----- missing values are replaced by NA \n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n MapReadPeak_Ext_NoCtrl: mapped reads in peaks ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n FRiP_Ext_NoCtrl: MapReadPeak_Ext_NoCtrl / UniqMapRead ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_NoCtrl: number of peaks (determined by p value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_NoCtrl(Q<0.05): number of peaks (determined by q value threshold of 0.05) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_NoCtrl(Q<0.01): number of peaks (determined by q value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n\n\n\n  MACS2 outputs corresponding to peaks with default MACS2 command ----- but here control input is present ----- missing values are replaced by NA \n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n MapReadPeak_Def_Ctrl: mapped reads in peaks ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n FRiP_Def_Ctrl: MapReadPeak_Def_Ctrl / UniqMapRead ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_Ctrl: number of peaks (determined by p value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_Ctrl(Q<0.05): number of peaks (determined by q value threshold of 0.05) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Def_Ctrl(Q<0.01): number of peaks (determined by q value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n\n\n\n  MACS2 outputs corresponding to peaks with --extsize option (recommended in existing ATAC seq pipeline) ----- here control input is provided ----- missing values are replaced by NA \n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n MapReadPeak_Ext_Ctrl: mapped reads in peaks ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n FRiP_Ext_Ctrl: MapReadPeak_Ext_Ctrl / UniqMapRead ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_Ctrl: number of peaks (determined by p value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_Ctrl(Q<0.05): number of peaks (determined by q value threshold of 0.05) ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n nPeak_Ext_Ctrl(Q<0.01): number of peaks (determined by q value threshold of 0.01) ")
writeLines(outtext, con=con, sep="\n")


# close output summary text file
close(con)

#===================
# now read the summary file once more, and plot different statistics
#===================
FinalDF <- read.table(outtextfile, header=T, sep="\t", stringsAsFactors=F)

# plot total number of reads for each sample
plotfile <- paste0(OutDir, '/TotalReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 3])
PlotScatter_Data(plotdf, "Total Reads", plotfile)

# plot fraction of mappable reads for each sample
plotfile <- paste0(OutDir, '/Fraction_MappableReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 5])
PlotScatter_Data(plotdf, "Fraction Mappable Reads", plotfile)

# plot fraction of mitochondrial reads for each sample
plotfile <- paste0(OutDir, '/Fraction_MitochondrialReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 12])
PlotScatter_Data(plotdf, "Fraction mitochondrial Reads", plotfile)

# plot fraction of uniquely mapped reads for each sample
plotfile <- paste0(OutDir, '/Fraction_UniqueMappReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 14])
PlotScatter_Data(plotdf, "Fraction unique mapped Reads", plotfile)

# plot fraction of low quality reads for each sample
plotfile <- paste0(OutDir, '/Fraction_LowQualReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 18])
PlotScatter_Data(plotdf, "Fraction low quality Reads", plotfile)

# plot fraction of duplicate reads for each sample
plotfile <- paste0(OutDir, '/Fraction_DuplicateReadCount_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 21])
PlotScatter_Data(plotdf, "Fraction duplicates", plotfile)

# plot NRF for each sample
plotfile <- paste0(OutDir, '/NRF_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 23])
PlotScatter_Data(plotdf, "NRF", plotfile)

# plot M1 for each sample
plotfile <- paste0(OutDir, '/M1_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 24])
PlotScatter_Data(plotdf, "M1", plotfile)

# plot M2 for each sample
plotfile <- paste0(OutDir, '/M2_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 25])
PlotScatter_Data(plotdf, "M2", plotfile)

# plot PBC1 for each sample
plotfile <- paste0(OutDir, '/PBC1_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 26])
PlotScatter_Data(plotdf, "PBC1", plotfile)

# plot PBC2 for each sample
plotfile <- paste0(OutDir, '/PBC2_Distribution.html')
plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, 27])
PlotScatter_Data(plotdf, "PBC2", plotfile)

# plot FRiP for each sample - no control, default MACS2 peaks
# provided the column is not filled with NA
colno <- 29
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/FRiP_Def_NoCtrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "FRiP_Def_NoCtrl", plotfile)
}

# plot number of peaks for each sample - FDR = 0.05 
# no control, MACS2 default peaks
# provided the column is not filled with NA
colno <- 31
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/NumPeak_Def_NoCtrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "NumPeak_Def_NoCtrl", plotfile)
}

# plot FRiP for each sample - no control, MACS2 Extsize peaks
# provided the column is not filled with NA
colno <- 34
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/FRiP_Ext_NoCtrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "FRiP_Ext_NoCtrl", plotfile)
}

# plot number of peaks for each sample - FDR = 0.05 
# no control, MACS2 Extsize peaks
# provided the column is not filled with NA
colno <- 36
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/NumPeak_Ext_NoCtrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "NumPeak_Ext_NoCtrl", plotfile)
}

# plot FRiP for each sample - with control, MACS2 default peaks
# provided the column is not filled with NA
colno <- 39
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/FRiP_Def_Ctrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "FRiP_Def_Ctrl", plotfile)
}

# plot number of peaks for each sample - FDR = 0.05 
# no control, MACS2 Extsize peaks
# provided the column is not filled with NA
colno <- 41
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/NumPeak_Def_Ctrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "NumPeak_Def_Ctrl", plotfile)
}

# plot FRiP for each sample - with control, MACS2 Extsize peaks
# provided the column is not filled with NA
colno <- 44
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/FRiP_Ext_Ctrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "FRiP_Ext_Ctrl", plotfile)
}

# plot number of peaks for each sample - FDR = 0.05 
# with control, MACS2 Extsize peaks
# provided the column is not filled with NA
colno <- 46
if ((colno %in% NA_ColList) == FALSE) {
	plotfile <- paste0(OutDir, '/NumPeak_Ext_Ctrl_Distribution.html')
	plotdf <- data.frame(X=FinalDF[, 2], Y=FinalDF[, colno])
	PlotScatter_Data(plotdf, "NumPeak_Ext_Ctrl", plotfile)
}

