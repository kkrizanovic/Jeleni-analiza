#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

# To enable importing from samscripts submodule
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, '/home/kkrizanovic/src/samscripts/src'))

from fastqparser import read_fastq

reads_folder = '/home/kkrizanovic/jeleni/fastq_filtered'
# check_sequences = ['TTCACTGTG', 'TTCGCTGTG']
check_sequences = ['TGTG']


for filename in os.listdir(reads_folder):
	sys.stderr.write('\n\nPROCESSING FILENAME : %s' % filename)
	reads_filename = os.path.join(reads_folder, filename)
	if filename.lower().endswith(".fastq"):
		[headers, seqs, quals] = read_fastq(reads_filename)
		for i in xrange(len(headers)):
			header = headers[i]
			seq = seqs[i]
			contains = False
			for subseq in check_sequences:
				if seq.find(subseq) != -1:
					contains = True

			if contains == False:
				sys.stderr.write('\nNO SUBSTRING FOUND: Filename = %s, Sequence = %s' % (filename, header))

sys.stderr.write("\n")


  
