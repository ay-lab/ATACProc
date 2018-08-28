#!/usr/bin/env python

"""
This program is for plotting ATAC seq peak distribution

Author: Sourya Bhattacharyya
Vijay-AY lab
"""

import matplotlib
matplotlib.use('Agg') 

import os
from optparse import OptionParser
# import re
import matplotlib.pyplot as plt
# import numpy as np

#-----------------------------------------------------
def parse_options():  
	parser = OptionParser()
		
	parser.add_option("-I", "--INPFILE", \
				type="string", \
				action="store", \
				dest="INP_TEXT_FILE", \
				default="", \
				help="Input TEXT file containing the Picard Insert size results")

	opts, args = parser.parse_args()
	return opts, args

#-----------------------------------------------------
"""
main function
"""
def main():  
	opts, args = parse_options()
	InpFile = opts.INP_TEXT_FILE

	k = InpFile.rfind('/')
	if (k == -1):
		InpDir = "./"
	else:
		InpDir = InpFile[:(k+1)]

	outdir = InpDir + "Plots"
	if (not os.path.exists(outdir)):
		os.makedirs(outdir)

	fragment_length_list = []
	aligned_read_count_list = []

	with open(InpFile) as f:
		for line in f.readlines():
			#curr_line_content = re.split(r'\s', line)
			curr_line_content = line.split()
			if (len(curr_line_content) == 2):
				str1 = str(curr_line_content[0])
				str2 = str(curr_line_content[1])
				if 0:
					print 'str1: ', str1, '  str2: ', str2
				if (len(str1) > 0) and (len(str2) > 0):
					if (str1[0].isdigit() == True) and (str2[0].isdigit() == True):
						fragment_length_list.append(int(str1))
						aligned_read_count_list.append(int(str2))

	if 0:
		print 'fragment_length_list: ',fragment_length_list
		print 'aligned_read_count_list: ',aligned_read_count_list

	total_read_count = sum(aligned_read_count_list)
	if 0:
		print 'total_read_count: ', total_read_count

	
	"""
	normalize the aligned read count
	dividing by the total no of reads and 
	with respect to unit fragment size
	"""
	for i in range(len(aligned_read_count_list)):
		aligned_read_count_list[i] = (aligned_read_count_list[i] * 1.0) / (fragment_length_list[i] * total_read_count)

	"""
	Now plot the statistics
	"""
	OutPlotFile = outdir + "/Fragment_plot_LINEAR.pdf"
	f = plt.figure()
	plt.plot(fragment_length_list, aligned_read_count_list, ls='-', lw=0.3, color='red')
	plt.xlim([0,1400])	# add - sourya - setting view for 1400 Kb
	plt.xlabel('Fragment length (bp)')
	plt.ylabel('Norm Read count')
	plt.title('ATAC seq - read density vs fragment length')
	f.savefig(OutPlotFile, bbox_inches='tight')


	OutPlotFile2 = outdir + "/Fragment_plot_LOG.pdf"
	f = plt.figure()
	plt.semilogy(fragment_length_list, aligned_read_count_list, ls='-', lw=0.3, color='red')
	plt.xlim([0,1400])	# add - sourya - setting view for 1400 Kb
	plt.xlabel('Fragment length (bp)')
	plt.ylabel('Norm Read count')
	plt.title('ATAC seq - read density vs fragment length')
	f.savefig(OutPlotFile2, bbox_inches='tight')


#-----------------------------------------------------
if __name__ == "__main__":
	main() 
