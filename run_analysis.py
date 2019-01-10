#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

fq_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/fastq'
temp_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/temp'
filtered_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/fastq_filtered'
spoa_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/spoa'

min_len = 293
max_len = 300

for filename in os.listdir(fq_folder):
	sys.stderr.write('\n\nPROCESSING FILENAME : %s' % filename)
	fq_filename = os.path.join(fq_folder, filename)
	if filename.lower().endswith(".fastq"):
		fname, fext = os.path.splitext(filename)
		tmp_filename = os.path.join(temp_folder, fname + '_tmp.fastq')
		flt_filename = os.path.join(filtered_folder, fname + '_filtered.fastq')
		consensus_filename = os.path.join(spoa_folder, fname + '.consensus')
		msa_filename = os.path.join(spoa_folder, fname + '.msa')
		matrix_filename = os.path.join(spoa_folder, fname + '.matrix')

		# Filter short sequences
		cmd = 'python /home/kkrizanovic/src/samscripts/src/fastqfilter.py minlen %d %s %s' % (min_len, fq_filename, tmp_filename)
		sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
		(status, output) = commands.getstatusoutput(cmd)

		#Filter long sequences
		cmd = 'python /home/kkrizanovic/src/samscripts/src/fastqfilter.py maxlen %d %s %s' % (max_len, tmp_filename, flt_filename)
		sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
		(status, output) = commands.getstatusoutput(cmd)

		# SPOA - consensus
		cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 0 -l 1 %s > %s' % (flt_filename, consensus_filename)
		sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
		(status, output) = commands.getstatusoutput(cmd)

		# SPOA - MSA
		cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 1 -l 1 %s > %s' % (flt_filename, msa_filename)
		sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
		(status, output) = commands.getstatusoutput(cmd)

		# SPOA - MSA Matrix
		cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 3 -l 1 %s > %s' % (flt_filename, matrix_filename)
		sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
		(status, output) = commands.getstatusoutput(cmd)		


  
