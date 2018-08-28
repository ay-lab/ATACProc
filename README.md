# ATACProc - a pipeline for processing ATAC-seq data

Developers
----------

Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Allergy and Immunology

La Jolla, San Diego, CA 92037, USA


#######################

ATACProc is a pipeline to analysis ATAC-seq data, starting from input Fastq/BAM files and generating alignment summary, various quality statistics, peak calling, and BigWig formatted tracks ready for visualization in UCSC genome browser. It also performs IDR analysis between a set of peak 
files or even a set of BAM alignment files (in which case, peaks are estimated first) 
corresponding to a set of biological or technical ATAC-seq replicates. The package also 
integrates a R based library ATACseqQC (Jianhong Ou et. al.) for QC 
analysis of ATAC-seq data.

#######################

Theory
----------

User can check the following papers or links for understanding ATAC-seq QCs:

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

Required packages for executing basic ATAC-seq pipeline
-------------------------------------------------------

When executing basic ATAC-seq pipeline, user should install following 
packages / libraries in the system:

1) Bowtie2 (we have used version 2.3.3.1) http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

2) samtools (we have used version 1.6) http://samtools.sourceforge.net/

3) PICARD tools (we have used 2.7.1 version) https://broadinstitute.github.io/picard/

4) Utilities "bedGraphToBigWig", "bedSort", "bigBedToBed", "hubCheck" and "fetchChromSizes" 
downloaded from UCSC repository. Executables corresponding to the linux system, 
for example, is provided in this link: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

5) deepTools (we have used version 2.0) https://deeptools.readthedocs.io/en/develop/

6) MACS2 (we have used version 2.1.1) https://github.com/taoliu/MACS

7) HOMER (we recommend using the latest version) http://homer.ucsd.edu/homer/

8) R environment (we have used 3.4.3)

9) The R package ATACseqQC (https://bioconductor.org/packages/release/bioc/html/ATACseqQC.html) 
and its associated dependencies.

User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. 
Some of these PATHS are also to be mentioned in a separate configuration file (mentioned below).


Required packages for executing IDR code
------------------------------------------

In addition, when user requires to execute the IDR code, 
following packages / libraries are to be installed in the system:

1) sambamba (we have used version 0.6.7) http://lomereiter.github.io/sambamba/

2) The package IDRCode (available in https://drive.google.com/file/d/0B_ssVVyXv8ZSX3luT0xhV3ZQNWc/view?usp=sharing). Unzip the archieve and store in convenient location. Path of this 
archieve is to be provided for executing IDR code.


Execution of basic ATAC-seq pipeline
------------------------------------

Current package includes a sample script file "pipeline_exec.sh". 
It conains sample commands required to 
invoke the main executable named "pipeline.sh", 
which is provided within the folder "bin".

In general, ATAC-seq pipeline (the executable "pipeline.sh") 
involves following command line options:

Options:

Mandatory parameters:

	-C  ConfigFile		    
	             A configuration file to be separately provided. Mandatory parameter. 
	             Current package includes a sample configuration file named "configfile". 
	             Details of the entries in this file are mentioned later.
	              
	-f  FASTQ1          
	            Read 1 (or forward strand) of paired-end sequencing data  [.fq|.gz|.bz2]. 
				Or, even an aligned genome (.bam file) can be provided.
	        
	-r  FASTQ2          
	            R2 of pair-end sequencing data [.fq|.gz|.bz2]. If not provided, and the -f parameter 
	            is not a BAM file, the input is assumed to be single ended.

	-n  PREFIX           
	            Prefix string of output files. For example, -n "TEST" means that the 
	            output filenames start with the string "TEST".

	-g  BOWTIE2_GENOME   
	            Bowtie2 indexed reference genome. Basically, the folder containing 
	            the bwt2 indices are to be provided. 
	            Mandatory parameter if user provides fastq files as input (-f and -r options).
				If user provides .bam files as input (-f option) then no need to provide this value.

	-d  OutDir 			  
	            Output directory which will contain all the results.

	-c  CONTROLBAM		 
	         	Control file(s) used for peak calling using MACS2. One or more 
				alignment files can be provided to be used 
				as a control. It may not be specified at all, in which 
				case MACS2 operates without any control. 
				Control file can be either in BAM or in  (tagalign.gz) format. 
				If multiple control files are provided, user needs to ensure that all of the 
				control files follow the same format (i.e. either all BAM or all TAGAlign).
				Example: -c control1.bam -c control2.bam puts two control files for using in MACS2.
		
				Conversion from BAM to TagAlign.gz format can be done using the script "TagAlign.sh" 
				provided within the folder "bin".
		
	-w 	BigWigGenome	 
				Reference genome which is used to convert BAM file to a BigWig file. 
				Used for visualization track creation purpose. 
				If -g option is enabled (i.e. the Bowtie2 index genome is provided) 
				then this option is not required. 
				Otherwise, this is a mandatory parameter. Allowed values are 'hg19' 
				(default), 'mm9', 'hg38', and 'mm10'.		
		
	-D  DEBUG_TXT		 
				Binary variable. If 1 (recommended), different statistics corresponding to 
				quality metrics and reads are printed. Useful when a summary of a large set 
				of ChIP-seq samples are to be generated.
		
	-q  MAPQ_THR		 
				Quality value threshold, below which the mapped reads are removed (Default 30).
		
	-p  PEAKCALLGENOMESIZE 
				genome size parameter for MACS2 peak calling ("hs", "mm", "ce", "dm": default "hs")

Optional parameters:

	-O 	Overwrite		 
				Binary variable. If 1, overwrites the existing files (if any). Default = 0.
						 
	-t  NUMTHREADS              
				Number of sorting, Bowtie2 mapping THREADS [Default = 1]. For parallel processing of Bowtie2, 
				user should specify > 1 value such as 4 ot 8.
		
	-m  MAX_MEM          
				Set max memory used for PICARD duplication removal [Default = 8G].
		
	-a  ALIGNVALIDMAX	 
				Set the number of (max) valid alignments which will be searched [Default = 4] 
				for Bowtie2.
		
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

	RPackageExec=
		Installed R package directory.
		Example: /home/sourya/R-3.4.3/bin/Rscript
		If left as blank, default Rscript installed in the system will be invoked.

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
		Can be downloaded from the link:
		http://hgdownload.cse.ucsc.edu/downloads.html
		Example: /home/sourya/genomes/Complete_Genome/hg38/hg38.fa
		
	RefChrAnnotFile=
		file containing reference UCSC annotation (.gtf format) 
		corresponding to the reference Chromosome.
		Can be downloaded from the link:
		http://hgdownload.cse.ucsc.edu/downloads.html
		Example: /home/sourya/genomes/Annotation/hg38/UCSC/hg38_UCSC_Annotation.gtf

Describing output of ATAC-seq analysis
-----------------------------------------

Within the folder "OutDir" (base directory containing all the outputs of 
current ATAC-seq analysis, following files (f) and folders (F) exist):

	F1: Alignment_MAPQ${MAPQ_THR}

		f0: Bowtie2_Init_Align.sam
			Initial alignment by Bowtie2 (if fastq files are provided as the input.)
		f1: UniqMappedRead.bam
			Initial alignment after unique mapping.
		f2: Bowtie2_del_Random.bam
			Alignment after deleting random reads.
		f3: Bowtie2_del_Mitch.bam: 
			After deleting mitochondrial reads.	
		f4: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.bam
			Sorted, and MAPQ thresholded alignment.
		f5: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.bam.bai
			Corresponding index.
		f6: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.rmdup.bam
			De-duplicated alignment (used for subsequent operations)
		f7: ${PREFIX}.align.sort.MAPQ${MAPQ_THR}.picard_metrics.txt
			Corresponding PICARD metrics log file.

		F1: ATAC_Seq_QC_New:

			This folder contains ATAC-seq QC metrics and summary statistics
			when computed by the package ATACseqQC. If the input data is single end, 
			only a file "Summary_QC.txt" is generated, containing the current 
			sequence specific summary metrics. For paired end input data, 
			the package also shows the library complexity, fragment size 
			distribution, nucleosome positioning, footprint analysis, and 
			motif entrichment. For details, please check the reference 
			manual of ATACseqQC in 
			https://bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html
			
	F2: Out_BigWig
		f1: ${PREFIX}.bw 
			bigwig file for track visualization.

	F3: Out_BigWig_NormCov:
		f1: ${PREFIX}_NormCov.bw
			bigwig file for track visualization (after normalizing the coverage).

	F4: MACS2_Default_*
		Contains peaks employing MACS2 with default parameters.
			f1: *.narrowPeak: narowpeak formatted output with P-value threshold of 0.01
			f2: *.narrowPeak_Q0.05filt: peaks with Q-value threshold of 0.05
			f3: *.narrowPeak_Q0.01filt: peaks with Q-value threshold of 0.01
			f4: *.broadPeak: broadpeak formatted output with P-value threshold of 0.01
			f5: *.broadPeak_Q0.05filt: peaks with Q-value threshold of 0.05
			f6: *.broadPeak_Q0.01filt: peaks with Q-value threshold of 0.01
			f7: out_FRiP.txt: FRIP statistics for the narrow and broad peaks.
			f8: Peak_Statistics.txt: number of peaks
			F9: Peak_Annotate_Q*:
				For Q-value thresholds of either 0.01 or 0.05, contains the 
				HOMER based annotations.
			
			In addition, files *.bb denote corresponding big-bd formatted peaks,
			useful for USCSC track visualization.
			
	F5: MACS2_Ext_*
		Contains peaks employing MACS2 with the parameters:
			--nomodel --nolambda --shift 0 --extsize 200
		File structure is similar as above.
		
	f8: out_NRF_MAPQ${MAPQ_THR}.txt
		Metric NRF
		
	f9: Read_Count_Stat.txt
		Read count statistics.

Summarizing a list of ATAC-seq analysis
---------------------------------------

Suppose, a directory "/home/sourya/Results" contain within it, the following folders: 
1, 2, 3, 4, ... Each corresponds to the output results for individual ATAC-seq samples.

To get a summarized list of performance metrics for these samples, use the script "ResSummary.r".

	Rscript ResSummary.r [positional_arguments]:

	1) OutBaseDir: 
		Directory under which results of all the different samples are stored 
		(like /home/sourya/Results as mentioned above)

	2) OutFile: 
		Output file which will store the summarized results. 
		Default: Results_All_Samples_Summary.xls.

Sample execution command:

Rscript ResSummary.r /home/sourya/ChIPResults/ OutSummary.xls


Command for executing IDR codes
---------------------------------

Current pipeline supports IDR analysis between either a list of ChIP-seq peak files 
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
			Prefix of output files. Default 'IDR_ChIP'.

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
			Prefix of output files. Default 'IDR_ChIP'.

	-c  CountPeak		 
			No of peaks in both replicates that will be compared for IDR analysis.
			Default 25000.
	
	-T 	Tagmentation	 
			Binary variable. If 1, the input is a ChiPMentation data 
			where the TAG Align files are created by 
			shifting the strands a bit. Default 0. 
			Tag align files are used for estimating peaks using MACS2.
	
	-C  CONTROLBAM		 
			Control file (in eiher .BAM or tagalign file in .gz format)	
			used to estimate the peaks from MACS2. User may leave this field 
			blank if no control file is available.

	A sample execution of this script is as follows:

	./IDR_SubSampleBAM_Main.sh -I inpfile1.bam -I inpfile2.bam -P /home/sourya/packages/idrCode/ -d /home/sourya/OutDir_IDR -n 'IDR_test' -c 25000 -T 1 -C control.bam


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

For any queries, please e-mail:

Sourya Bhattacharyya (sourya@lji.org)

Ferhat Ay (ferhatay@lji.org)

Pandurangan Vijayanand (vijay@lji.org)

