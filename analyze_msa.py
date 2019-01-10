#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

import numpy as np

if (len(sys.argv) < 2) or (len(sys.argv) > 3):
	sys.stderr.write('Run program: %s [MSA file] <factor>!\n' % sys.argv[0])
	exit(1)

max_coverage = 0
factor = 0.20

if len(sys.argv) == 3:
	factor = float(sys.argv[2])

m_filename = sys.argv[1]

msafile = open(m_filename, 'r')
lines = msafile.readlines()
# First line in MSA file is a header, and is followed by line pairs representing
# sequence name and alignment
numseq = (len(lines) - 1)/2
seqs_dict = {}
seqlen = len(lines[2])	# Length of the first alignment, all should be of the same length

# Initializing matrix lines
lineA = np.zeros(seqlen)
lineC = np.zeros(seqlen)
lineG = np.zeros(seqlen)
lineT = np.zeros(seqlen)
lineN = np.zeros(seqlen)

for i in xrange(numseq):
	seqname = lines[i*2 + 1][1:]
	seqalign = lines[i*2 + 2]
	if seqname in seqs_dict:
		sys.stderr.write('\nERROR: seq %s already in the dictonary!' % seqname)
		sys.stderr.write('\n%s' % seqalign)
	seqs_dict[seqname] = seqalign

	for i in xrange(len(seqalign)):
		c = seqalign[i]
		if c == 'A':
			lineA[i] += 1
		if c == 'C':
			lineC[i] += 1
		if c == 'G':
			lineG[i] += 1
		if c == 'T':
			lineT[i] += 1
		if c == '-':
			lineN[i] += 1

# Calculate maximum coverage - equal to number of sequences
max_coverage = numseq

lower_limit = factor*max_coverage
upper_limit = (1-factor)*max_coverage

for i in xrange(len(lineA)-1):
	countA = lineA[i]
	countC = lineC[i]
	countG = lineG[i]
	countT = lineT[i]
	countN = lineN[i]

	# Determining candidates for variants
	if ((countA >= lower_limit and countA <= upper_limit) or \
	    (countC >= lower_limit and countC <= upper_limit) or \
	    (countG >= lower_limit and countG <= upper_limit) or \
	    (countT >= lower_limit and countT <= upper_limit) or \
	    (countN >= lower_limit and countN <= upper_limit)):
		sys.stdout.write('Pos: %d, countA: %d, countC: %d, countG: %d, countT: %d, count-: %d\n' % (i, countA, countC, countG, countT, countN))

sys.stdout.write('\n')



