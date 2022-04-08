#!/usr/bin/env Rscript

#==================================
# footprinting using HINT-ATAC package
# http://www.regulatory-genomics.org/hint/introduction/

# author: Sourya Bhattacharyya
# Vijay-AY lab

# check 
# http://www.regulatory-genomics.org/hint/tutorial/
#==================================

library(optparse)

#===========================================================
option_list = list(

	make_option(c("--AllRead"), type="character", default=NULL, help="Alignment file containing all reads."),
	make_option(c("--NFRRead"), type="character", default=NULL, help="Alignment file containing nucleosome free regions (NFR) reads."),
	make_option(c("--NFRANDNuclRead"), type="character", default=NULL, help="Alignment file containing nucleosome free regions (NFR) plus all nucleosome (1M, 2M, 3M) merged reads."),
	make_option(c("--RefGenome"), type="character", default=NULL, help="Reference genome name."),
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory to contain the motif."),
	make_option(c("--PE"), type="integer", action="store", default=0, help="If 1, indicates paired end input data. Default = 0"),
	make_option(c("--FP"), type="integer", action="store", default=1, help="Footorinting option. Value can be 1 (default), 2, or 3. (1): footoprint using the nucleosome free reads (NFR) will be computed. Default setting. Best for default ATAC-seq protocol (check Li et. al. Genome Biology 2019). 2: footoprint using the nucleosome free reads (NFR) and also the nucleosome containing reads (NFR + 1N + 2N + 3N ...) will be computed (two different footprint outputs - time consuming). Best for Omni-ATAC protocol (check Li et. al. Genome Biology 2019). (3): footoprint using NFR, NFR with nucleosome reads, and all reads will be computed (three different footprint outputs - highly time consuming). Default = 1"),
	make_option(c("--MotifPeak"), type="character", default=NULL, help="Peak or summit file which was used by HOMER to generate corresponding motifs. Mandatory parameter.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# create the output directory
system(paste("mkdir -p", opt$OutDir))

# prefix string of output file name
OUTPREFIX <- 'footprints_HINT_ATAC'


##===========
## processing all reads
## only if FP option > 2
##===========
if (opt$FP > 2) {
	if (!is.null(opt$AllRead)) {
		curroutdir <- paste0(opt$OutDir, '/all')
		system(paste("mkdir -p", curroutdir))
		if ((file.exists(paste0(curroutdir, '/', OUTPREFIX, '.bed')) == FALSE) | (file.exists(paste0(curroutdir, '/', OUTPREFIX, '.info')) == FALSE)) {
			if (opt$PE == 1) {
				cat(sprintf("\n start footprint HINT ATAC PE reads - all reads"))
				system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --paired-end --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$AllRead, opt$MotifPeak))
			} else {
				system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$AllRead, opt$MotifPeak))
			}
		}
		# now call motif matching for the obtained footprints
		# JASPAR database is used by default for motif finding
		# 10% random background region is tested - using the option --rand-proportion 10
		cat(sprintf("\n start motifanalysis of HINT ATAC - all reads"))
		motifoutdir <- paste0(curroutdir, '/motifanalysis_matching_out')
		system(paste("mkdir -p", motifoutdir))
		system(paste("rgt-motifanalysis matching --organism ", opt$RefGenome, " --rand-proportion 10 --input-files ", paste0(curroutdir, '/', OUTPREFIX, '.bed'), " --output-location ",  motifoutdir))
	}
}

##===========
## processing NFR and nucleosome reads (1N, 2N, ...)
## only if FP option > 1
##===========
if (opt$FP > 1) {
	if (!is.null(opt$NFRANDNuclRead)) {
		curroutdir <- paste0(opt$OutDir, '/NFRANDNucl')
		system(paste("mkdir -p", curroutdir))
		if ((file.exists(paste0(curroutdir, '/', OUTPREFIX, '.bed')) == FALSE) | (file.exists(paste0(curroutdir, '/', OUTPREFIX, '.info')) == FALSE)) {	
			if (opt$PE == 1) {
				cat(sprintf("\n start footprint HINT ATAC PE reads - nucleosome free and nucleosome reads"))
				system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --paired-end --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$NFRANDNuclRead, opt$MotifPeak))
			} else {
				system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$NFRANDNuclRead, opt$MotifPeak))
			}
		}
		# now call motif matching for the obtained footprints
		# JASPAR database is used by default for motif finding	
		# 10% random background region is tested - using the option --rand-proportion 10
		cat(sprintf("\n start motifanalysis of HINT ATAC - nucleosome free and nucleosome reads"))
		motifoutdir <- paste0(curroutdir, '/motifanalysis_matching_out')
		system(paste("mkdir -p", motifoutdir))	
		system(paste("rgt-motifanalysis matching --organism ", opt$RefGenome, " --rand-proportion 10 --input-files ", paste0(curroutdir, '/', OUTPREFIX, '.bed'), " --output-location ",  motifoutdir))
	}
}

##===========
## processing NFR reads
## default option
##===========
if (!is.null(opt$NFRRead)) {
	curroutdir <- paste0(opt$OutDir, '/NFR')
	system(paste("mkdir -p", curroutdir))
	if ((file.exists(paste0(curroutdir, '/', OUTPREFIX, '.bed')) == FALSE) | (file.exists(paste0(curroutdir, '/', OUTPREFIX, '.info')) == FALSE)) {
		if (opt$PE == 1) {
			cat(sprintf("\n start footprint HINT ATAC PE reads - nucleosome free reads"))
			system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --paired-end --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$NFRRead, opt$MotifPeak))
		} else {
			system(paste("rgt-hint footprinting --atac-seq --organism ", opt$RefGenome, " --output-location ", curroutdir, " --output-prefix ", OUTPREFIX, " ", opt$NFRRead, opt$MotifPeak))
		}
	}
	# now call motif matching for the obtained footprints
	# JASPAR database is used by default for motif finding	
	# 10% random background region is tested - using the option --rand-proportion 10
	cat(sprintf("\n start motifanalysis of HINT ATAC - nucleosome free reads"))
	motifoutdir <- paste0(curroutdir, '/motifanalysis_matching_out')
	system(paste("mkdir -p", motifoutdir))	
	system(paste("rgt-motifanalysis matching --organism ", opt$RefGenome, " --rand-proportion 10 --input-files ", paste0(curroutdir, '/', OUTPREFIX, '.bed'), " --output-location ",  motifoutdir))
}



