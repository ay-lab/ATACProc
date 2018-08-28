## assign_multimappers.py

# this code is used to process the multi mapped reads
# used after the Bowtie2 output is filtered to remove all the unmapped reads


# code by - Sourya Bhattacharyya
# taken from the standard ATAC seq pipeline

import sys 
import random 
import argparse

# function to parse the arguments
def parse_args(): 
	# Gives options
	parser = argparse.ArgumentParser(description='Saves reads below an alignment threshold and discards all others')
	parser.add_argument('-k', help='Alignment number cutoff') 
	parser.add_argument('--paired-end', dest='paired_ended', action='store_true', help='Data is paired-end')

	# processing the input arguments
	# and return the input parameters to the main function
	args = parser.parse_args()
	alignment_cutoff = int(args.k) 
	paired_ended = args.paired_ended
	return alignment_cutoff, paired_ended

# main function
if __name__ == "__main__": 
	
	# Runs the filtering step of choosing multimapped reads
	[alignment_cutoff, paired_ended] = parse_args()

	# when paired ended input, the cutoff is adjusted
	if paired_ended:
		alignment_cutoff = int(alignment_cutoff) * 2

	# Store each line in sam file as a list of reads,
	# where each read is a list of elements to easily 
	# modify or grab things
	current_reads = []
	current_qname = ''

	# processing individual lines
	for line in sys.stdin:
		read_elems = line.strip().split('\t')
		if read_elems[0].startswith('@'): 
			sys.stdout.write(line)
			continue
	
		# Keep taking lines that have the same qname 
		if read_elems[0] == current_qname:
			# Add line to current reads 
			current_reads.append(line)
			pass

		else:
			# Discard if there are more than the alignment cutoff 
			if len(current_reads) >= alignment_cutoff: 
				current_reads = [line]
				current_qname = read_elems[0]
			elif len(current_reads) > 0:
				# Just output all reads, which are then filtered with samtools
				for read in current_reads:
					sys.stdout.write(str(read))
				# And then discard 
				current_reads = [line] 
				current_qname = read_elems[0] 
			else:
				# First read in file 
				current_reads.append(line) 
				current_qname = read_elems[0]


