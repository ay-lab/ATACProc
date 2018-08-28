#!/usr/bin/env python

"""
This program computes the distribution of peak fragment length and the number of aligned reads (read density)
used for benchmarking the ATAC-seq pipeline


Author: Sourya Bhattacharyya
Vijay-AY lab
"""
"""
these two lines force matplotlibv to not choose any X-windows 
This should be declared very first
see http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
"""
import matplotlib
matplotlib.use('Agg') 

import os
from optparse import OptionParser
import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
	parser = OptionParser()
		
	parser.add_option("-I", "--INPFILE", \
				type="string", \
				action="store", \
				dest="INP_BED_FILE", \
				default="", \
				help="Input BED file containing the MACS2 peak detection results")

	parser.add_option("-R", "--REFFILE", \
				type="string", \
				action="store", \
				dest="REF_BAM_FILE", \
				default="", \
				help="Reference BAM alignment file")  
			
	opts, args = parser.parse_args()
	return opts, args
  
#-----------------------------------------------------
"""
main function
"""
def main():  
	opts, args = parse_options()

	INP_BEDFILE = opts.INP_BED_FILE
	REF_BAMFILE = opts.REF_BAM_FILE

	k = INP_BEDFILE.rfind('/')
	if (k == -1):
		INP_BEDFILE_DIR = "./"
		Inp_BED_only_filename = INP_BEDFILE
	else:
		INP_BEDFILE_DIR = INP_BEDFILE[:(k+1)]
		Inp_BED_only_filename = INP_BEDFILE[(k+1):]

	k1 = Inp_BED_only_filename.rfind('.')
	if (k1 == -1):
		OutDN = Inp_BED_only_filename
	else:
		OutDN = Inp_BED_only_filename[(k1+1):]

	"""
	final output directory which will store the plots and data
	"""
	OutDir_Name = INP_BEDFILE_DIR + OutDN
	if (not os.path.exists(OutDir_Name)):
		os.makedirs(OutDir_Name)

	fragment_length_list = []
	aligned_read_count_list = []

	temp_filename = OutDir_Name + "/temp.bed"
	fp_temp = open(temp_filename, "w")

	"""
	scan each line of the input bed file, and note the peak fragment length
	also note the number of input reads (from BAM file) which are mapped in this peak
	"""
	with open(INP_BEDFILE) as fp_inp:
		for line in fp_inp:
			# note the peak fragment length
			curr_line_content = re.split(r'\s', line)
			if 0:
				print '\n\n Current line: ', line, '  Contents: ', curr_line_content
			peak_fragment_len = int(curr_line_content[2]) - int(curr_line_content[1]) + 1
			if 0:
				print 'peak_fragment_len: ', peak_fragment_len
			# write the line to a temporary bed file
			fp_temp.seek(0, os.SEEK_SET)
			fp_temp.write(line)
			# now count the number of mapped reads to this peak
			sys_cmd = "samtools view -cL " + str(temp_filename) + " " + str(REF_BAMFILE)
			read_count = int((subprocess.Popen(sys_cmd, stdout=subprocess.PIPE, shell=True)).stdout.read())
			if 0:
				print 'read_count: ', read_count
			# now append the values in designated lists
			# maintain sorted lists
			n = len(fragment_length_list)
			if (n == 0):
				# very first element
				fragment_length_list.append(peak_fragment_len)
				aligned_read_count_list.append(read_count)
			else:
				flag = False
				for i in xrange((n-1), -1, -1):
					if (peak_fragment_len == fragment_length_list[i]):
						aligned_read_count_list[i] = aligned_read_count_list[i] + read_count
						flag = True
						break
					elif (peak_fragment_len > fragment_length_list[i]):
						if (i == (n-1)):
							fragment_length_list.append(peak_fragment_len)
							aligned_read_count_list.append(read_count)
							flag = True
						else:
							fragment_length_list.insert(i+1, peak_fragment_len)
							aligned_read_count_list.insert(i+1, read_count)
							flag = True
						break

				if (flag == False):
					# condition for insertion at the first location
					fragment_length_list.insert(0, peak_fragment_len)
					aligned_read_count_list.insert(0, read_count)

	# close the temporary file
	fp_temp.close()

	# remove the temporary bed file
	os.system("rm " + temp_filename)

	"""
	open a text file with two columns
	first column will show the peak fragment length
	second column displays the read count
	the plot file is stored in the same directory containing macs2 results
	"""
	plot_data_textfile = OutDir_Name + "/plot.txt"
	fp = open(plot_data_textfile, "w")
	fp.write("Peak_Length" + "\t" + "Read_Count" + "")
	for i in range(len(fragment_length_list)):
		fp.write("\n" + str(fragment_length_list[i]) + "\t" + str(aligned_read_count_list[i]))
	fp.close()

	# """
	# create a read count list which will contain the no of read count in 1K scale
	# """
	# read_count_list_1K_scale = [((aligned_read_count_list[i] * 1.0) / 1000) for i in range(len(aligned_read_count_list))]

	"""
	Now plot the statistics
	"""
	# OutPlotFile = OutDir_Name + "/test_1K_Scale.pdf"
	# f = plt.figure()
	# plt.plot(fragment_length_list, read_count_list_1K_scale, ls='-', lw=2.0)
	# plt.xlabel('Fragment length (bp)')
	# plt.ylabel('Norm Read count')
	# plt.title('ATAC seq - read density vs fragment length')
	# f.savefig(OutPlotFile, bbox_inches='tight')

	OutPlotFile2 = OutDir_Name + "/test_LOG_Scale.pdf"
	f = plt.figure()
	plt.semilogy(fragment_length_list, np.exp(-np.asarray(aligned_read_count_list)/5.0), ls='-', lw=2.0)
	#plt.semilogy(fragment_length_list, np.exp(-np.asarray(aligned_read_count_list)/5.0), ls='-', lw=2.0)
	# plt.plot(fragment_length_list, aligned_read_count_list, ls='-', lw=2.0)
	# plt.yscale('log', basey=10)
	plt.xlabel('Fragment length (bp)')
	plt.ylabel('Norm Read count')
	plt.title('ATAC seq - read density vs fragment length')
	f.savefig(OutPlotFile2, bbox_inches='tight')

#-----------------------------------------------------
if __name__ == "__main__":
	main() 

