# ATACProc - a pipeline for processing ATAC-seq data

Devloper: Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Immunology, CA 92037, USA


#######################

ATACProc is a pipeline to analyze ATAC-seq data. Currently datasets involving one of the four reference genomes, namely hg19, hg38, mm9 and mm10 are supported. Important features of this pipeline are:

1) Supports single or paired-end fastq or BAM formatted data.

2) Generates alignment summary and QC statistics.

3) Peak calls using MACS2, for multiple FDR thresholds (0.01 and 0.05)

4) Generating raw and coverage normalized BigWig tracks for visualizing the data in UCSC genome browser.

5) Irreproducible Discovery Rate (IDR) analysis (https://github.com/nboley/idr) between a set of peak calls or even a set of input alignment (BAM) files (in which case, peaks are estimated first) corresponding to a set of biological or technical ATAC-seq replicates. 

6) **New in version 2.0:** Support discarding reads falling in blacklisted genomic regions

7) **New in version 2.0:** Support extracting nucleosome free reads (NFR), one or more nucleosome containing regions (denoted as +1M), for TF footprinting analysis.

8) **New in version 2.0:** Compatibility to the package ATAQV (https://github.com/ParkerLab/ataqv) for generating summary statistics across a set of samples.

#######################

Release notes
-----------------

**Version 2.1 - July 2020**

Minor change of picard duplicate removal syntax, according to the picard tool version 2.8.14 
We recommend using this (or later) versions

**Version 2.0 - November 2019**

1) Included TF footprinting, optional discarding of blacklisted genomic regions, motif analysis

2) Updated summary statistics incorporating support for ATAQV package (https://github.com/ParkerLab/ataqv)

3) Discarded R package ATACseqQC (https://bioconductor.org/packages/release/bioc/html/ATACseqQC.html) and corresponding operations, mainly due to its time complexity and reliability issues.


*Version 1.0 - July 2018:*

1) Released first version of ATAC-seq pipeline, supporting generation of QC metrics, peak calls, signal tracks for visualizing in UCSC genome browser. 

2) Also supports IDR between a set of peaks / alignments for a set of replicates.


Theory
----------

Papers / links for understanding ATAC-seq QCs:

1) https://github.com/crazyhottommy/ChIP-seq-analysis  (very useful; contains many papers 
and links for understanding ChIP-seq and ATAC-seq data)

2) https://www.encodeproject.org/data-standards/terms/#library

3) https://www.biostars.org/p/187204/

4) http://seqanswers.com/forums/archive/index.php/t-59219.html

5) https://github.com/kundajelab/atac_dnase_pipelines

6) https://github.com/ParkerLab/bioinf525#sifting

7) https://github.com/taoliu/MACS/issues/145

8) https://www.biostars.org/p/207318/

9) https://www.biostars.org/p/209592/

10) https://www.biostars.org/p/205576/


Understanding peak calling

1) https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137

Understanding TF footprinting

1) https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2

Understanding IDR analysis

1) https://github.com/nboley/idr



Installation
-------------

Following packages / libraries should be installed before running this pipeline:

1) Python 2.7 

2) R environment (we have used 3.4.3)

	User should also install the following R packages, by running the following command inside R prompt:

	install.packages(c(“optparse”, “ggplot2”, “data.table”, “plotly”))

	Also user needs to install the bioconductor package GenomicRanges <https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html>

3) Bowtie2 (we have used version 2.3.3.1) <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>

4) samtools (we have used version 1.6) <http://samtools.sourceforge.net/>

5) PICARD tools (we have used 2.8.14 version now; previously we were using version 2.7.1) <https://broadinstitute.github.io/picard/>

6) Utilities "bedGraphToBigWig", "bedSort", "bigBedToBed", "hubCheck" and "fetchChromSizes" - to be downloaded from UCSC repository <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>

7) deepTools (we have used version 2.0) <https://deeptools.readthedocs.io/en/develop/>

8) MACS2 (we have used version 2.1.1) https://github.com/taoliu/MACS

9) HOMER (we recommend using the latest version) http://homer.ucsd.edu/homer/

10) The package *ataqv* (https://github.com/ParkerLab/ataqv). User needs to download the GitHub release (.tar.gz) file in a convenient location, extract it, and provide corresponding path in a configuration file (mentioned below).

11) Regulatory genomics toolbox (https://www.regulatory-genomics.org/) 

	First user needs to install the module *RGT* using the following commands:

		pip install --user cython numpy scipy
		pip install --user RGT

	A folder *rgtdata* would be created inside the home directory. Next step is to configure that folder by typing the following commands:

		cd ~/rgtdata
		python setupGenomicData.py --hg19
		python setupGenomicData.py --hg38
		python setupGenomicData.py --mm9
		python setupGenomicData.py --mm10

		(Note: it is better to run the last four commands together in a qsub / cluster environment, otherwise it'll be time consuming).


	Then, user needs to set up the motif configuration data, via executing the following commands (preferable to run in qsub / cluster environment)

		cd ~/rgtdata
		python setupLogoData.py --all


**User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. Alternatively, installation PATHS for some of these packages are to be mentioned in a separate configuration file (described below)**

**Following packages / libraries are to be installed for executing IDR code**

9) sambamba (we have used version 0.6.7) <http://lomereiter.github.io/sambamba/>

10) IDRCode (https://drive.google.com/file/d/0B_ssVVyXv8ZSX3luT0xhV3ZQNWc/view?usp=sharing). User should unzip the archieve and store in convenient location. Path of this archieve is to be provided for executing IDR code.



Execution
----------

User should first clone this pipeline in a convenient location, using the following command: 

git clone https://github.com/ay-lab/ATACProc.git

A sample script "pipeline_exec.sh" contains basic execution commands, to invoke the main executable "pipeline.sh" (located inside the folder "bin"). The executable has the following command line options:

Options:

Mandatory parameters:

	-C  ConfigFile		    
         	Configuration file to be separately provided. Mandatory parameter. Current package includes four sample configuration files named "configfile_*" corresponding to the reference genomes hg19, hg38, mm9 and mm10. Detailed description of the entries in this configuration file are mentioned later.
	              
	-f  FASTQ1          
        	Read 1 (or forward strand) of paired-end sequencing data  [.fq|.gz|.bz2]. Or, even an aligned genome (.bam file; single or paired end alignment) can be provided.
	        
	-r  FASTQ2          
            R2 of pair-end sequencing data [.fq|.gz|.bz2]. If not provided, and the -f parameter is not a BAM file, the input is assumed to be single ended.

	-n  PREFIX           
            Prefix string of output files. For example, -n "TEST" means that the output filenames start with the string "TEST". Generally, sample names with run ID, lane information, etc. can be used as a prefix string.

	-g  BOWTIE2_GENOME   
            Bowtie2 indexed reference genome. Basically, the folder containing bwt2 indices (corresponding to the reference genome) are to be provided. Mandatory parameter if the user provides fastq files as input (-f and -r options). If user provides .bam files as an input (-f option) then this field is optional.

	-d  OutDir 			  
            Output directory to store the results for the current sample.

	-c  CONTROLBAM		 
         	Control file(s) used for peak calling using MACS2. One or more alignment files can be provided to be used as a control. It may not be specified at all, in which case MACS2 operates without any control. Control file can be either in *BAM* or in *tagalign.gz* format (the standalone script *bin/TagAlign.sh* in this repository converts BAM file to tagalign.gz format). For multiple control files, they all are required to be of the same format (i.e. either all BAM or all tagalign.gz). Example: -c control1.bam -c control2.bam puts two control files for using in MACS2.
		
	-w BigWigGenome	 
			Reference genome as a string. Allowed values are hg19 (default), hg38, mm9 and mm10. If -g option is enabled (i.e. the Bowtie2 index genome is provided), this field is optional. Otherwise, mandatory parameter.				
		
	-D  DEBUG_TXT		 
			Binary variable. If 1 (recommended), dumps QC statistics. For a set of samples, those QC statistics can be used later to profile QC variation among different samples.				
		
	-O 	Overwrite		 
			Binary variable. If 1, overwrites the existing files (if any). Default = 0.

	-F 	Footprint 	 	
			This flag specifies the footprinting option. Value can be 1 (default), 2, or 3
			1: footoprint using the nucleosome free reads (NFR) will be computed. 
			   Default setting. Best for default ATAC-seq protocol (check Li et. al. Genome Biology 2019)
			2: footoprint using the nucleosome free reads (NFR) and also the nucleosome containing reads (NFR + 1N + 2N + 3N ...) 
			   will be computed (two different footprint outputs - time consuming). 
			   Best for Omni-ATAC protocol (check Li et. al. Genome Biology 2019)
 			3: footoprint using NFR, NFR with nucleosome reads, and all reads will be computed 
			   (three different footprint outputs - highly time consuming).	
			   
Optional parameters:
	-q  MAPQ_THR		 
			Mapping quality threshold for bowtie2 alignment. Aligned reads with quality below this threshold are discarded. Default = 30. 
		 
	-t  NUMTHREADS              
			Number of sorting, Bowtie2 mapping THREADS [Default = 1]. If multiprocessing core is available, user should specify values > 1 such as 4 or 8, for faster execution of Bowtie2.
		
	-m  MAX_MEM          
			Set max memory used for PICARD duplication removal [Default = 8G].
		
	-a  ALIGNVALIDMAX	 
			Set the number of (max) valid alignments which will be searched [Default = 4] for Bowtie2.
		
	-l  MAXFRAGLEN 		 
			Set the maximum fragment length to be used for Bowtie2 alignment [Default = 2000]
			

Entries in the configuration file (first parameter)
---------------------------------------------------

The configuration file follows the format parameter=value

And is to be filled with the following entries:

	picardexec=
		Path of Picard tool executable
		Example: /home/sourya/packages/picard-tools/picard-tools-2.7.1/picard.jar

	HOMERPath=
		Path of HOMER (after installation)
		Example: /home/sourya/packages/HOMER/bin/

	DeepToolsDir=
		Path of deepTools executable
		Example: /home/sourya/packages/deepTools/deepTools2.0/bin/	

	NarrowPeakASFile=
		file (SQL) required to convert the narrowPeak file to the bigBed format
		Download the file from this link (and save):
		https://genome-source.gi.ucsc.edu/gitlist/kent.git/blob/master/src/hg/lib/encode/narrowPeak.as
		Specify the location of this downloaded file:
		Example: /home/sourya/genomes/chrsize/narrowPeak.as

	BigNarrowPeakASFile=
		file (SQL) required to convert the bignarrowPeak file to the bigBed format
		Download the file from this link (and save):
		https://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.as
		Specify the location of this downloaded file:
		Example: /home/sourya/genomes/chrsize/bigNarrowPeak.as
		
	BroadPeakASFile=
		file (SQL) required to convert the broadPeak file to the bigBed format
		Download the file from this link (and save):
		https://genome-source.gi.ucsc.edu/gitlist/kent.git/blob/master/src/hg/lib/encode/broadPeak.as
		Specify the location of this downloaded file:
		Example: /home/sourya/genomes/chrsize/broadPeak.as
		
	RefChrSizeFile=
		files containing chromosome size information
		two column file storing the size of individual chromosomes
		Downloaded from the link (depends on the reference Chromosome employed):
		For example, the hg38.chrom.sizes file for the hg38 database is located at 
		http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes.
		Alternatively, Use the "fetchChromSizes" script from the UCSC repository 
		to get the appropriate chromosome size file.
		Specify the location of this downloaded file:
		Example: /home/sourya/genomes/chrsize/hg38.chrom.sizes
		
	RefChrFastaFile=
		Fasta file of the reference Chromosome. 
		Can be downloaded from the link: http://hgdownload.cse.ucsc.edu/downloads.html
		Example: /home/sourya/genomes/Complete_Genome/hg38/hg38.fa
		
	RefChrAnnotFile=
		file containing reference genome specific annotation (.gtf format). 
		To be downloaded from the following links:
		hg38: ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/
		hg19: ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
		mm9: ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/
		mm10: ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/
		Example: /home/sourya/genomes/Annotation/hg38/hg38.gtf

	BlackListFile=
		file containing blacklisted regions corresponding to the reference genome. 
		To be downloaded from the link: https://github.com/Boyle-Lab/Blacklist/tree/master/lists (v2)
		File can be gzipped or normal text format.
		*Note: This parameter is optional.*
		Example: /home/sourya/genomes/BlackListed_Regions/hg38-blacklist.v2.bed

	ATAQVPath=
		Path of ataqv package (https://github.com/ParkerLab/ataqv) executable. 
		User needs to download the GitHub release (.tar.gz) file, extract it, and provide the ataqv executable path here.
		Example: /home/sourya/packages/ataqv/ataqv-1.0.0/bin/ataqv

	TSSFile=
		File containing TSS information for the reference genome. Obtained using the gene annotation (GTF) file.
		Example: /home/sourya/genomes/Annotation/hg38/hg38_TSS.gtf

	
	The last parameter, *TSSFile*, needs a special mention. User can apply the following awk script to the reference genome annotation file (indicated in the parameter *RefChrAnnotFile*) to produce a file with TSS information.

	Assuming user has downloaded the reference genome specific gene annotation file using one of the ftp links provided above, when the reference genome is either hg19, hg38 or mm10, user can apply the following awk script to obtain a TSS file (input_TSS.gtf) from the gene annotation file (input.gtf) (Note: it is always best to check the .gtf file format) :

		awk -F'[\t]' '{if ((substr($1,1,1)!="#") && ($3=="transcript")) {if ($7=="+") {print "chr"$1"\t"$4"\t"$4"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9} else {print "chr"$1"\t"$5"\t"$5"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9}}}' input.gtf > input_TSS.gtf

	When the reference genome is mm9, user can apply the following script (it is best to check the .gtf file format):

		awk -F'[\t]' '{if ((substr($1,1,1)!="#") && ($3=="exon")) {if ($7=="+") {print "chr"$1"\t"$4"\t"$4"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9} else {print "chr"$1"\t"$5"\t"$5"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9}}}' mm9.gtf > mm9_TSS.gtf

Describing output of ATAC-seq pipeline
-----------------------------------------

Within the folder *OutDir* (specified by the configuration option -d) following files (f) and folders (F) exist:

	F1: Alignment_MAPQ${MAPQ_THR}

		f1-1: Bowtie2_Init_Align.sam
			Initial alignment by Bowtie2 (if fastq files are provided as the input.)
		f1-2: UniqMappedRead.bam
			Uniquely mapped reads.
		f1-3: Bowtie2_del_Random.bam
			Alignment after excluding reads from chromosomes other than autosomal chromosomes, chrX, and chrM.
		f1-4: Bowtie2_del_Mitch.bam: 
			Alignment after excluding reads from chrM.
		f1-5: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.bam
			Sorted, and MAPQ thresholded alignment.
		f1-6: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.rmdup.bam
			De-duplicated alignment (used for subsequent operations)
		f1-7: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.picard_metrics.txt
			PICARD metrics log file corresponding to the duplicate removal operation.
		f1-8: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}_TN5_Shift.bam
			**New in version 2.0:** De-duplicated reads with shifted forward (+4bp) and reverse strands (-5bp) by Tn5 transposase. Used to extract the nucleosome free and nucleosome containing regions.
		f1-9: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}_TN5_Shift.bed
			**New in version 2.0:** Bed converted f7, used for MACS2 peak calling.
		f1-10: NucleosomeFree.bam
			**New in version 2.0:** Alignment with nucleosome free regions (NFR)
		f1-11: mononucleosome.bam
			**New in version 2.0:** Alignment with mononucleosome fragments
		f1-12: dinucleosome.bam
			**New in version 2.0:** Alignment with dinucleosome fragments
		f1-13: trinucleosome.bam
			**New in version 2.0:** Alignment with trinucleosome fragments
		f1-14: Merged_nucleosome.bam
			**New in version 2.0:** File containing fragments of nucleosome free and one or more nucleosomes (denoted as NFR +1M, in the HINT-ATAC genome biology paper). Generated by merging files f1-10 to f1-13.

	F2: Out_BigWig
		f2-1: ${PREFIX}.bw 
			bigwig file for track visualization.

	F3: Out_BigWig_NormCov:
		f3-1: ${PREFIX}_NormCov.bw
			bigwig file for track visualization (after normalizing the coverage). Recommended to use this file for visualizing tracks in UCSC genome browser.

	F4: MACS2_Ext_*
		Contains peaks employing MACS2 with the parameters:
			--nomodel --nolambda --shift -100 --extsize 200 --keep-dup all --call-summits
		*Note: this parameter is recommended for ATAC-seq, as mostly followed in existin studies.*

		If the folder name is "*_No_Control", no control BAM file was used to infer the peaks. Otherwise, if the folder name is "*_With_Control", one or more control alignment files were used for inferring the peaks.

			f4-1: *.narrowPeak: narrow peaks with p-value threshold of 0.01
			f4-2: *.narrowPeak_Q0.05filt: narrow peaks with FDR (q-value) threshold = 0.05
			f4-3: *.narrowPeak_Q0.01filt: narrow peaks with FDR threshold = 0.01
			f4-4: *.broadPeak: broad peaks with p-value threshold of 0.01
			f4-5: *.broadPeak_Q0.05filt: broad peaks with FDR threshold = 0.05
			f4-6: *.broadPeak_Q0.01filt: broad peaks with FDR threshold = 0.01
			f4-7: out_FRiP.txt: FRiP (fraction of reads in peaks) statistics for the narrow and broad peaks.
			f4-8: Peak_Statistics.txt: number of peaks in different settings.
			F4-9: Peak_Annotate_Q*:
				HOMER based annotations corresponding to the narrow peaks inferred by the corresponding FDR threshold (0.01 or 0.05). Contains the following files:
				f4-9-1: Out_Summary.log: summary text file containing HOMER annotation.
				f4-9-2: Annotated_Peak_Q*filt.txt: Detailed HOMER annotation of the corresponding peaks.
				f4-9-3: Pie_Chart_Peak_Annotation.pdf: pie chart of peaks containing different annotations.
				f4-9-4: Peak_TSS_Distance.pdf: Histogram of distance between peaks and closest TSS
			f4-10: Files of *.bb extension are big-bed formatted peaks, used to visualize those peaks in UCSC tracks.

	F5: MACS2_Default_*
		Contains peaks employing default MACS2 parameters. (generally not used for ATAC-seq processing, but we've kept it for comparison).
		File and folder structure is similar as F4.

	f8: out_NRF_MAPQ${MAPQ_THR}.txt
		Metric NRF
		
	f9: Read_Count_Stat.txt
		Read count statistics.

	F10: QC_ataqv_ParkerLab_Test
		**New in version 2.0:** Folder containing the summary .json files generated by the package ATAQV, which for diferent samples, can be combined to put a summary statistic and displayed in a Web browser.

	F11: TSS_Enrichment_Peaks
		**New in version 2.0:** Processes the narrow peaks from the folder F4, and computes the TSS enrichment of these peaks. The underlying file structure is:

		MACS2_Ext_*${CONTROLSTR}/macs2_narrowPeak_Q${FDRTHR}filt_Offset_${OFFSETVAL}/${PEAKTYPE}/*.pdf

		where, 
			${CONTROLSTR}: "*_No_Control" or "*_With_Control", depending on the use of control BAM file in inferring the peaks.
			${FDRTHR}: FDR threshold. Can be either 0.01 or 0.05
			${OFFSETVAL}: can be either 1000 (1 Kb) or 5000 (5Kb) (1 Kb or 5 Kb regions surrounding TSS are checked for computing TSS enrichment).
			${PEAKTYPE}: can be either "Complete_Peaks" (means complete set of peaks are experimented), "Promoter_Peaks" (means peaks located within 5 Kb of a TSS site are only considered), or "Enhancer_Peaks" (peaks excluding the promoter peaks).


	F12: Motif_MACS2_Ext_*${CONTROLSTR}_narrowPeak_Q${FDRTHR}filt
		**New in version 2.0:** TF footorinting analysis corresponding to the ChIP-seq peaks stored in F4. Here, ${CONTROLSTR} is either "*_No_Control" or "*_With_Control", depending on the use of control BAM file in inferring the peaks. ${FDRTHR} is either 0.01 or 0.05.

		The principle is to extract the peak summits and surroundings (by some bp, defined as an offset) and compute the TF footprinting regions and underlying motifs within these regions.

		Within this folder, the file structure is as follows:
		Motif_${PEAKS_ANALYZED}_SummitOffset_${OFFSET}/Footprint_HINT_ATAC/${READTYPE}/footprints_HINT_ATAC.bed

		where, 
			${PEAKS_ANALYZED}: can be "Complete_Peaks" (means complete set of peaks) or "Peaks_PvalThr_50" (means peaks with -log10(p-value) > 50 are only considered).
			${OFFSET}: can be either 200 or 500, means the summit +/- offset bp regions are accounted for TF footprinting.
			${READTYPE}: can be one of the following:
			 	"all" (means all de-duplicated reads in the file f1-8 considered), 
			 	"NFR" (means only nucleosome free reads in the file f1-10 are considered), 
			 	"NFRANDNucl" (means NFR regions and +1M reads, indicated by the file f1-14, are considered).

		 	The output file in each occasion, "footprints_HINT_ATAC.bed", contains the TF footprinting regions.


Summarizing a set of ATAC-seq samples
---------------------------------------

Suppose, a directory "/home/sourya/Results" contain within it, the following folders: 
1, 2, 3, 4, ... each corresponding to the output for processing individual ATAC-seq samples.

To get a summarized list of performance metrics for these samples, use the script *Analysis/ResSummary.r*, using the following syntax.

	Rscript ResSummary.r --BaseDir ${BaseDir} --OutDir ${OutDir}

	where,
	1) ${BaseDir}: 
		Directory containing results of all ATAC-seq sample analysis 		
		(like /home/sourya/Results as mentioned above). Mandatory parameter.

	2) ${OutDir}: 
		Output directory to contain the summarized results. Default: current working directory.

	For details of ATAC-seq QC measures, user may check this link:
	https://www.encodeproject.org/atac-seq/

	Upon executing the R script, the following files are created within the specified ${OutDir}:

		1) Results_All_Samples_Summary.txt: summarized statistics for all samples
		2) Field_Description.txt: Summary description of individual fields / parameters.
		3) TotalReadCount_Distribution.html: To be loaded in any web browser. Plot depicting the distribution of total reads for all samples.
		4) Fraction_MappableReadCount_Distribution.html: Fraction of mappability for all samples.
		5) Fraction_MitochondrialReadCount_Distribution.html: Fraction of mitochondrial reads for all samples.
		6) Fraction_UniqueMappReadCount_Distribution.html: Fraction of unique mappability for all samples.
		7) Fraction_LowQualReadCount_Distribution.html: Fraction of low quality reads for all samples.
		8) Fraction_DuplicateReadCount_Distribution.html: Fraction of duplicate reads for all samples.
		9) NRF_Distribution.html: NRF for all samples.
		10) M1_Distribution.html: M1 metric for all samples.
		11) M2_Distribution.html: M2 metric for all samples.
		12) PBC1_Distribution.html: PBC1 metric for all samples.
		13) PBC2_Distribution.html: PBC2 metric for all samples.
		14) FRiP_Def_NoCtrl_Distribution.html: FRiP statistics for MACS2 peaks with default command, and without using any control BAM files.
		15) NumPeak_Def_NoCtrl_Distribution.html: Number of MACS2 peaks with default command, and without using any control BAM files.
		16) FRiP_Ext_NoCtrl_Distribution.html: FRiP statistics for MACS2 peaks with --Extsize option (recommended), and without using any control BAM files.
		17) NumPeak_Ext_NoCtrl_Distribution.html: Number of MACS2 peaks with --Extsize option (recommended), and without using any control BAM files.
		18) FRiP_Def_Ctrl_Distribution.html: FRiP statistics for MACS2 peaks with default command, and when one or more control BAM files are used.
		19) NumPeak_Def_Ctrl_Distribution.html: Number of MACS2 peaks with default command, and when one or more control BAM files are used.
		20) FRiP_Ext_Ctrl_Distribution.html: FRiP statistics for MACS2 peaks with --Extsize option (recommended), and when one or more control BAM files are used.
		21) NumPeak_Ext_Ctrl_Distribution.html: Number of MACS2 peaks with --Extsize option (recommended), and when one or more control BAM files are used.

Command for executing IDR codes
---------------------------------

Current pipeline supports IDR analysis between either a list of ATAC-seq peak files 
or between a list of alignment (BAM) files. In the second case, first the BAM files 
are analyzed and subsampled to contain equal number of reads (minimum number of reads 
contained in the inputs), and subsequently, peaks are estimated from these 
(subsampled) BAM files using MACS2. These peaks are then applied for IDR analysis.

The script "sample_IDRScript.sh" included within this package 
shows calling following two functions (both are included within the folder 
"IDR_Codes"):

	1) IDRMain.sh

	2) IDR_SubSampleBAM_Main.sh

	The first script, IDRMain.sh, performs IDR between two or more 
	input peak files (we have used peaks estimated from MACS2). The parameters 
	corresponding to this script are as follows:

	-I  InpFile        	 
			A list of input peak files (obtained from MACS2 - in .narrowPeak or .narrowPeak.gz format). 
			At least two peak files are required. 
	
	-P 	PathIDRCode		 
			Path of the IDRCode package (Kundaje et. al. after its installation). 
			Please check the "Required packages" section for the details.

	-d  OutDir 		 	 
			Output directory (absolute path preferred) which will store the IDR results.

	-n 	PREFIX 			 
			Prefix of output files. Default 'IDR_ATAC'.

	A sample execution of this script is as follows:

	./IDRMain.sh -I peak1.narrowPeak -I peak2.narrowPeak -I peak3.narrowPeak -P /home/sourya/packages/idrCode/ -d /home/sourya/OutDir_IDR -n 'IDR_test'



	The second script, IDR_SubSampleBAM_Main.sh, takes input of two or more BAM files, 
	estimates peaks from these BAM files, and then performs IDR analysis. The parameters 
	corresponding to this script are as follows:

	-I  InpFile        	 
			A list of input BAM files. At least two BAM files are required. 
	
	-P 	PathIDRCode		 
			Path of the IDRCode package (Kundaje et. al. after its installation). 
			Please check the "Required packages" section for the details.

	-d  OutDir 		 	 
			Output directory (absolute path preferred) which will store the IDR results.

	-n 	PREFIX 			 
			Prefix of output files. Default 'IDR_ATAC'.

	-c  CountPeak		 
			No of peaks in both replicates that will be compared for IDR analysis.
			Default 25000.
		
	-C  CONTROLBAM		 
			Control file (in eiher .BAM or tagalign file in .gz format)	
			used to estimate the peaks from MACS2. User may leave this field 
			blank if no control file is available.

	A sample execution of this script is as follows:

	./IDR_SubSampleBAM_Main.sh -I inpfile1.bam -I inpfile2.bam -P /home/sourya/packages/idrCode/ -d /home/sourya/OutDir_IDR -n 'IDR_test' -c 25000 -C control.bam


Describing output of IDR analysis
----------------------------------

In the specified output directory "OutDir" mentioned in the IDR script, following 
files (f) and folders (F) exist:

	F1: Folders of the name $i$_and_$j$ where 0 <= i < N and 1 <= j <= N, where N is 
	the number of replicates analyzed. Individual folders contain results for 
	pairwise IDR analysis. For example, folder 0_and_1 contain IDR analysis 
	for the sample 0 (first replicate) and the sample 1 (second replicate).

	f1 : "Replicate_Names.txt" : names of the replicate samples used for IDR analysis.

	f2: Input_Peak_Statistics.txt: number of peaks and the peak containing replicates.

	f3: IDR_Batch_Plot-plot.pdf: final IDR plot. Here individual pairs (whose results 
		are stored in the above mentioned folders) are numbered 1, 2, ...
		Consideing N = 3, the number of pairs possible is also 3. Here, 
		the number 1 denotes the folder (pair) 0_and_1, 
		2 denotes the folder (pair) 0_and_2, and 3 denotes the 
		folder (pair) 1_and_2.


Contact
-----------

For any queries, please generate a GitHub issue, or alternatively, e-mail us:

Sourya Bhattacharyya (sourya@lji.org)

Ferhat Ay (ferhatay@lji.org)

Pandurangan Vijayanand (vijay@lji.org)

