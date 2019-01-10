#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

import numpy as np

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, '/home/kkrizanovic/src/samscripts/src'))

from fastqparser import read_fastq

if len(sys.argv) != 5:
	sys.stderr.write('Run program: %s position [MSA file] [READS file] [RESULTS folder]!\n' % sys.argv[0])
	exit(1)

position = int(sys.argv[1])
msa_filename = sys.argv[2]
reads_filename = sys.argv[3]
results_folder = sys.argv[4]

if not os.path.exists(results_folder):
	os.mkdir(results_folder)

msafile = open(msa_filename, 'r')
lines = msafile.readlines()

numseq = (len(lines) - 1)/2
reads_dict = {}
seqlen = len(lines[2])	# Length of the first alignment, all should be of the same length

# Looking through MSA
for i in xrange(numseq):
	seqname = lines[i*2 + 1][1:-1]
	seqalign = lines[i*2 + 2][:-1]

	base = seqalign[position]
	if base == '-':
		base = 'N'
	if base in reads_dict:
		reads_dict[base].append(seqname)
	else:
		reads_dict[base] = [seqname]

# Loading reads
[headers, seqs, quals] = read_fastq(reads_filename)
r_fname, r_fext = os.path.splitext(os.path.basename(reads_filename))

# Seaparating reads into files
for base, readname_list in reads_dict.iteritems():
	filename = os.path.join(results_folder, r_fname + '_' + base + r_fext)
	file = open(filename, 'w')
	for i in xrange(len(headers)):
		header = headers[i]
		seq = seqs[i]
		qual = quals[i]

		if header in readname_list:
			file.write('@%s\n%s\n+\n%s\n' % (header, seq, qual))
	file.close()
	


