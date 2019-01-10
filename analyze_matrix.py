#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

if len(sys.argv) != 2:
	sys.stderr.write('Run program: %s [matrix file]!\n' % sys.argv[0])
	exit(1)

max_coverage = 0
factor = 0.2

m_filename = sys.argv[1]

matrixfile = open(m_filename, 'r')
lines = matrixfile.readlines()

lineA = lines[1].split(' ')
lineC = lines[2].split(' ')
lineG = lines[3].split(' ')
lineT = lines[4].split(' ')
lineN = lines[5].split(' ')

# Calculate maximum coverage
# Disregarding the last element because its '\n'
try:
	for s_count in lineA[:-1]:
		count = int(s_count)
		if count > max_coverage:
			max_coverage = count
except Exception:
	import pdb
	pdb.set_trace()

lower_limit = factor*max_coverage
upper_limit = (1-factor)*max_coverage

for i in xrange(len(lineA)-1):
	countA = int(lineA[i])
	countC = int(lineC[i])
	countG = int(lineG[i])
	countT = int(lineT[i])
	countN = int(lineN[i])

	# Determining candidates for variants
	if ((countA >= lower_limit and countA <= upper_limit) or \
	    (countC >= lower_limit and countC <= upper_limit) or \
	    (countG >= lower_limit and countG <= upper_limit) or \
	    (countT >= lower_limit and countT <= upper_limit) or \
	    (countN >= lower_limit and countN <= upper_limit)):
		sys.stdout.write('\nPos: %d, countA: %d, countC: %d, countG: %d, countT: %d, count-: %d' % (i, countA, countC, countG, countT, countN))

sys.stdout.write('\n')



