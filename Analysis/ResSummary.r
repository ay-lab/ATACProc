#!/usr/bin/env Rscript

#==================================
# used to print the summary result statistics 
# from a collection of ATAC seq samples
#==================================

# author: Sourya Bhattacharyya
# Vijay-AY lab

#=========================================
# execution command:
# Rscript summary_res_script.r $basedir $OutFile
# parameters:
# 1) $basedir: directory under which results of all samples are present 
# 2) $OutFile: output file (preferably in excel format) which stores the summary results
# 		if no such file is provided, a file "Results_All_Samples_Summary.xls" is created 
# 		and placed under the directory $basedir
#=========================================

args <- commandArgs(TRUE)
if(length(args)<1)	{
	q("no")
}

# directory under which all the results have been stored
# all the subdirectories indicate different samples
# corresponding to the experimental condition mentioned 
# in the main directory
baseresdir <- args[1]

# template name of the directories containing MACS2 results
# either default parameters are used
# or extsize based parameters are employed
MACS2_def_dir <- 'MACS2_Default_Tag'
MACS2_ext_dir <- 'MACS2_Ext_Tag'
MACS2_noduprem_ext_dir <- 'MACS2_NoDupRem_Align_Ext_Tag'

# template name of the folders (peak outputs) depending on whether control (input) are used for peak detection
Ctrl_0_Fold <- '_No_Control'
Ctrl_1_Fold <- '_with_Control'

# file formats of NRF, read statistics, FRiP and Peak statistics
# which are present for every sample
ReadStatFileNameFmt <- 'Read_Count_Stat.txt'
NRFfilenamefmt <- 'out_NRF'
FRiPFileNameFmt <- 'out_FRiP.txt'
PeakCountFileFmt <- 'Peak_Statistics.txt'

# check if the last character of the output base directory is '/'
# unless append that character
if (substr(baseresdir,nchar(baseresdir),nchar(baseresdir)) != '/') {
	baseresdir <- paste0(baseresdir, '/')
}

# output excel file (consolidated summary of different performance measures)
if (length(args) > 1) {
	outtextfile <- args[2]
} else {
	outtextfile <- paste0(baseresdir, 'Results_All_Samples_Summary.xls')
}

# open the output text file
con <- file(outtextfile, "w")

# header string of the output summary file
cat(paste(c("Dir", "TotRead", "MappableRead", "%MappRead", "RandomDelRead", "%RandomDelRead", "DelMitchndRead", "%DelMitch", "UniqMapRead", "%UniqMapp", "QualThr", "%QualThr", "rmDupRead", "%rmDupRead", "UniqMapPos", "NRF", "M1", "M2", "PBC1", "PBC2", "MapReadPeak_Def_NoCtrl", "FRiP_Def_NoCtrl", "nPeak_Def_NoCtrl", "nPeak_Def_NoCtrl(Q<0.05)", "nPeak_Def_NoCtrl(Q<0.01)", "MapReadPeak_Ext_NoCtrl", "FRiP_Ext_NoCtrl", "nPeak_Ext_NoCtrl", "nPeak_Ext_NoCtrl(Q<0.05)", "nPeak_Ext_NoCtrl(Q<0.01)", "MapReadPeak_Def_Ctrl", "FRiP_Def_Ctrl", "nPeak_Def_Ctrl", "nPeak_Def_Ctrl(Q<0.05)", "nPeak_Def_Ctrl(Q<0.01)", "MapReadPeak_Ext_Ctrl", "FRiP_Ext_Ctrl", "nPeak_Ext_Ctrl", "nPeak_Ext_Ctrl(Q<0.05)", "nPeak_Ext_Ctrl(Q<0.01)"), collapse='\t'), file=con)

#=================================
# process individual directories under the main results directory
dir.list <- list.dirs(path = baseresdir, full.names = FALSE, recursive = FALSE)
for (dr in dir.list) {

	# cat(sprintf("\n Examining directory: %s \n", dr))

	#==================
	# following file stores the count of reads throughout various stages of filtering
	#==================
	ReadStatFile <- paste0(baseresdir, dr, "/", ReadStatFileNameFmt)
	if (file.exists(ReadStatFile) && (file.access(ReadStatFile, 4) == 0)){

		cat(sprintf("\n Found the file: %s \n", ReadStatFile))

		x <- readLines(ReadStatFile)
		lastline <- strsplit(x[length(x)], "\t")[[1]]

		# the line can have 8 or 7 fields
		# 8 fields if fastq file is used in the pipeline
		# 7 fields if already aligned file is used in the pipeline
		nfields <- length(lastline)
		cat(sprintf("\n No of fields in the read statistics file: %s ", nfields))

	 	TotRead <- lastline[nfields-6]	#lastline[2]
	 	MappableRead <- lastline[nfields-5]	#lastline[3]
	 	Percent_Mappable_Read <- as.double(MappableRead) / as.double(TotRead)
	 	RandomDelRead <- lastline[nfields-4]
	 	Percent_RandomDelRead <- as.double(RandomDelRead) / as.double(TotRead)
	 	Del_Mitochondrial_Read <- lastline[nfields-3]	#lastline[4]
	 	Percent_Remain_Mit_Read <- as.double(Del_Mitochondrial_Read) / as.double(TotRead)
	 	UniqMappedRead <- lastline[nfields-2]	#lastline[5]
	 	Percent_UniqMappedRead <- as.double(UniqMappedRead) / as.double(TotRead)
	 	ReadQualThr <- lastline[nfields-1]	#lastline[6]
	 	Percent_QualThr <- as.double(ReadQualThr) / as.double(TotRead)
	 	Dupl_Rem_Read <- lastline[nfields]	#lastline[7]
	 	Percent_DuplRead <- as.double(Dupl_Rem_Read) / as.double(TotRead)	  			
	} else {
		TotRead <- 'NA'
		MappableRead <- 'NA'
		Percent_Mappable_Read <- 'NA'
		RandomDelRead <- 'NA'
		Percent_RandomDelRead <- 'NA'
		Del_Mitochondrial_Read <- 'NA'
		Percent_Remain_Mit_Read <- 'NA'
		UniqMappedRead <- 'NA'
		Percent_UniqMappedRead <- 'NA'
		ReadQualThr <- 'NA'
		Percent_QualThr <- 'NA'
		Dupl_Rem_Read <- 'NA'
		Percent_DuplRead <- 'NA'
	}

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
	 	} else {
		 	UniqMappedPos <- 'NA'	
		 	NRF_val <- 'NA'
		 	M1 <- 'NA'
		 	M2 <- 'NA'
		 	PBC1 <- 'NA'
		 	PBC2 <- 'NA'
	 	}
	} else {
	 	UniqMappedPos <- 'NA'		
	 	NRF_val <- 'NA'
	 	M1 <- 'NA'
	 	M2 <- 'NA'
	 	PBC1 <- 'NA'
	 	PBC2 <- 'NA'
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
	} else {
		MappedReadPeak_def_noctrl <- 'NA'
		FRiP_def_noctrl <- 'NA'
	}

	PeakCount_TextFile_def_noctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_0_Fold, "/", PeakCountFileFmt)
	if (file.exists(PeakCount_TextFile_def_noctrl) && (file.access(PeakCount_TextFile_def_noctrl, 4) == 0)){
	 	x <- readLines(PeakCount_TextFile_def_noctrl)
	 	lastline <- strsplit(x[length(x)], "\t")[[1]]
	 	TotPeak_def_noctrl <- lastline[1]
	  	TotPeak_Q_Five_Pct_def_noctrl <- lastline[2]
	  	TotPeak_Q_One_Pct_def_noctrl <- lastline[3]
	} else {
		TotPeak_def_noctrl <- 'NA'
		TotPeak_Q_Five_Pct_def_noctrl <- 'NA'
		TotPeak_Q_One_Pct_def_noctrl <- 'NA'
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
	} else {
		MappedReadPeak_ext_noctrl <- 'NA'
		FRiP_ext_noctrl <- 'NA'
	}

	PeakCount_TextFile_ext_noctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_0_Fold, "/", PeakCountFileFmt)
	if (file.exists(PeakCount_TextFile_ext_noctrl) && (file.access(PeakCount_TextFile_ext_noctrl, 4) == 0)){
	 	x <- readLines(PeakCount_TextFile_ext_noctrl)
	 	lastline <- strsplit(x[length(x)], "\t")[[1]]
	 	TotPeak_ext_noctrl <- lastline[1]
	  	TotPeak_Q_Five_Pct_ext_noctrl <- lastline[2]
	  	TotPeak_Q_One_Pct_ext_noctrl <- lastline[3]
	} else {
		TotPeak_ext_noctrl <- 'NA'
		TotPeak_Q_Five_Pct_ext_noctrl <- 'NA'
		TotPeak_Q_One_Pct_ext_noctrl <- 'NA'
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
	} else {
		MappedReadPeak_def_ctrl <- 'NA'
		FRiP_def_ctrl <- 'NA'
	}

	PeakCount_TextFile_def_ctrl <- paste0(baseresdir, dr, "/", MACS2_def_dir, Ctrl_1_Fold, "/", PeakCountFileFmt)
	if (file.exists(PeakCount_TextFile_def_ctrl) && (file.access(PeakCount_TextFile_def_ctrl, 4) == 0)){
	 	x <- readLines(PeakCount_TextFile_def_ctrl)
	 	lastline <- strsplit(x[length(x)], "\t")[[1]]
	 	TotPeak_def_ctrl <- lastline[1]
	  	TotPeak_Q_Five_Pct_def_ctrl <- lastline[2]
	  	TotPeak_Q_One_Pct_def_ctrl <- lastline[3]
	} else {
		TotPeak_def_ctrl <- 'NA'
		TotPeak_Q_Five_Pct_def_ctrl <- 'NA'
		TotPeak_Q_One_Pct_def_ctrl <- 'NA'
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
	} else {
		MappedReadPeak_ext_ctrl <- 'NA'
		FRiP_ext_ctrl <- 'NA'
	}

	PeakCount_TextFile_ext_ctrl <- paste0(baseresdir, dr, "/", MACS2_ext_dir, Ctrl_1_Fold, "/", PeakCountFileFmt)
	if (file.exists(PeakCount_TextFile_ext_ctrl) && (file.access(PeakCount_TextFile_ext_ctrl, 4) == 0)){
	 	x <- readLines(PeakCount_TextFile_ext_ctrl)
	 	lastline <- strsplit(x[length(x)], "\t")[[1]]
	 	TotPeak_ext_ctrl <- lastline[1]
	  	TotPeak_Q_Five_Pct_ext_ctrl <- lastline[2]
	  	TotPeak_Q_One_Pct_ext_ctrl <- lastline[3]
	} else {
		TotPeak_ext_ctrl <- 'NA'
		TotPeak_Q_Five_Pct_ext_ctrl <- 'NA'
		TotPeak_Q_One_Pct_ext_ctrl <- 'NA'
	}

	#==============================================
	# write the statistic in the results log file
	cat('\n', paste(c(basename(dr), TotRead, MappableRead, Percent_Mappable_Read, RandomDelRead, Percent_RandomDelRead, Del_Mitochondrial_Read, Percent_Remain_Mit_Read, UniqMappedRead, Percent_UniqMappedRead, ReadQualThr, Percent_QualThr, Dupl_Rem_Read, Percent_DuplRead, UniqMappedPos, NRF_val, M1, M2, PBC1, PBC2, MappedReadPeak_def_noctrl, FRiP_def_noctrl, TotPeak_def_noctrl, TotPeak_Q_Five_Pct_def_noctrl, TotPeak_Q_One_Pct_def_noctrl, MappedReadPeak_ext_noctrl, FRiP_ext_noctrl, TotPeak_ext_noctrl, TotPeak_Q_Five_Pct_ext_noctrl, TotPeak_Q_One_Pct_ext_noctrl, MappedReadPeak_def_ctrl, FRiP_def_ctrl, TotPeak_def_ctrl, TotPeak_Q_Five_Pct_def_ctrl, TotPeak_Q_One_Pct_def_ctrl, MappedReadPeak_ext_ctrl, FRiP_ext_ctrl, TotPeak_ext_ctrl, TotPeak_Q_Five_Pct_ext_ctrl, TotPeak_Q_One_Pct_ext_ctrl), collapse='\t'), file=con) 

}

#============================
# a few summary statements for the output text file
#============================
outtext <- paste0("\n\n\n\n\n *** Important parameters ***** \n\n\n")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  TotRead: number of reads in individual fastq file(s)")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  MappableRead and %MappRead: number (and percentage) of reads mappable to the reference genome. May not be uniquely mappable reads.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  RandomDelRead and %RandomDelRead: number (and percentage) of reads remaining after deleting random chromosomes such as chr1_*, chr2_*, chrUN, ....")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  DelMitchndRead and %DelMitch: number (and percentage) of reads remaining after deleting the mitochondrial reads")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  UniqMapRead and %UniqMapp: number (and percentage) of reads uniquely mapped to the reference genome.")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  QualThr and %QualThr: number (and percentage) of reads remaining after deleting the aligned reads below quality threshold (default 30)")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n  rmDupRead and %rmDupRead: number (and percentage) of reads remaining after deleting the deplicated reads")
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

outtext <- paste0("\n\n MapReadPeak_Def_NoCtrl: mapped reads in peaks ")
writeLines(outtext, con=con, sep="\n")

outtext <- paste0("\n\n FRiP_Def_NoCtrl: MapReadPeak_Def_NoCtrl / UniqMapRead ")
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





