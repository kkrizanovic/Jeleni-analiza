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

analyze_msa = '/home/kkrizanovic/jeleni/analyze_msa.py'

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

# Separating reads into files, but only if there is a sufficient number of them!
max_coverage = numseq
factor = 0.2
for base, readname_list in reads_dict.iteritems():
	if len(readname_list) < factor*max_coverage:
		continue
	sep_filename = os.path.join(results_folder, r_fname + '_' + base + r_fext)
	file = open(sep_filename, 'w')
	for i in xrange(len(headers)):
		header = headers[i]
		seq = seqs[i]
		qual = quals[i]

		if header in readname_list:
			file.write('@%s\n%s\n+\n%s\n' % (header, seq, qual))
	file.close()
	
	# Calculate consensus and MSA and run MSA analysis on the new reads file
	# SPOA - consensus
	consensus_filename = os.path.join(results_folder, r_fname + '_' + base + ('_pos%d' % position) + '.consensus')
	cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 0 -l 1 %s > %s' % (sep_filename, consensus_filename)
	sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
	(status, output) = commands.getstatusoutput(cmd)

	# SPOA - MSA
	new_msa_filename = os.path.join(results_folder, r_fname + '_' + base + ('_pos%d' % position) + '.msa')
	cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 1 -l 1 %s > %s' % (sep_filename, new_msa_filename)
	sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
	(status, output) = commands.getstatusoutput(cmd)

	# Analyze msa
	new_varpos_filename = os.path.join(results_folder, r_fname + '_' + base + ('_pos%d' % position) + '.varpos')
	cmd = 'python %s %s > %s' % (analyze_msa, new_msa_filename, new_varpos_filename)
	sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
	(status, output) = commands.getstatusoutput(cmd)

sys.stderr.write('\n\n')
