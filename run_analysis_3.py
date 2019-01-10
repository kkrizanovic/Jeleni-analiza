#! /usr/bin/python

import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

sys.path.append('/home/kkrizanovic/src/samscripts/src')

from fastqparser import read_fastq

results_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/results10'
analyze_msa = '/home/kkrizanovic/jeleni/analyze_msa.py'

factor = 0.10

def parse_vpline(vpline):
	# stripping stuff from end of line 
	vpline = vpline.rstrip('\n ')

	# parsing position
	tpos = vpline.find(',')
	pos = int(vpline[len('pos: '):tpos])
	vpline = vpline[tpos+2:]

	# parsing countA
	tpos = vpline.find(',')
	cntA = int(vpline[len('countA: '):tpos])
	vpline = vpline[tpos+2:]

	# parsing countC
	tpos = vpline.find(',')
	cntC = int(vpline[len('countC: '):tpos])
	vpline = vpline[tpos+2:]

	# parsing countG
	tpos = vpline.find(',')
	cntG = int(vpline[len('countG: '):tpos])
	vpline = vpline[tpos+2:]

	# parsing countT
	tpos = vpline.find(',')
	cntT = int(vpline[len('countT: '):tpos])
	vpline = vpline[tpos+2:]

	# parsing count-
	cntN = int(vpline[len('countA: '):])

	return (pos, cntA, cntC, cntG, cntT, cntN)

for filename1 in os.listdir(results_folder):
	subfolder = os.path.join(results_folder, filename1)
	for filename2 in os.listdir(subfolder):
		if filename2.lower().endswith(".msa"):
			fname, fext = os.path.splitext(filename2)
			msa_filename = os.path.join(subfolder, filename2)
			sys.stderr.write('\n\nPROCESSING FILENAME : %s' % msa_filename)
			reads_filename = os.path.join(subfolder, fname + '.fastq')
			varpos_filename = os.path.join(subfolder, fname + '.varpos')
			consensus_filename = os.path.join(subfolder, fname + '.consensus')
			
			# open varpos file, read all the lines and determine the best split point
			# The one which is close to 50/50 split
			vp_file = open(varpos_filename, 'r')
			vp_lines = vp_file.readlines()
			if len(vp_lines) == 0:
				sys.stderr.write('\nWARNING: no varpos found for: %s' % filename2)
				continue
			firstline = vp_lines[0]
			if firstline.rstrip('\n ') == '':
				sys.stderr.write('\nWARNING: no varpos found for: %s' % filename2)
				continue

			# Looking for position with minimal largest count, this represents the split closest to 50/50
			(bestpos, cntA, cntC, cntG, cntT, cntN) = parse_vpline(firstline)
			best_max_cnt = max(cntA, cntC, cntG, cntT, cntN)

			for vpline in vp_lines[1:]:
				# End on empty line
				if vpline.rstrip('\n ') == '':
					break
				(pos, cntA, cntC, cntG, cntT, cntN) = parse_vpline(firstline)
				max_cnt = max(cntA, cntC, cntG, cntT, cntN)
				if max_cnt < best_max_cnt:
					bestpos = pos
					best_max_cnt = max_cnt

			vp_file.close()


			pos = fname.rfind('_')
			if pos == -1:
				sys.stderr.write('\nWARNING: invalid filename: %s' % fname)
				continue
			base = fname[pos-1]
			splitfolder = os.path.join(subfolder, 'split_%s_pos%d' % (base, bestpos))
			if not os.path.exists(splitfolder):
				os.mkdir(splitfolder)
		
			# Go through MSA file and split reads into groups
			msafile = open(msa_filename, 'r')
			lines = msafile.readlines()

			numseq = (len(lines) - 1)/2
			reads_dict = {}
			seqlen = len(lines[2])	# Length of the first alignment, all should be of the same length

			# Looking through MSA
			for i in xrange(numseq):
				seqname = lines[i*2 + 1][1:-1]
				seqalign = lines[i*2 + 2][:-1]

				base = seqalign[bestpos]
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
			# factor = 0.2
			for base, readname_list in reads_dict.iteritems():
				if len(readname_list) < factor*max_coverage:
					continue
				sep_filename = os.path.join(splitfolder, r_fname + '_' + base + ('_pos%d' % bestpos) + r_fext)
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
				consensus_filename = os.path.join(splitfolder, r_fname + '_' + base + ('_pos%d' % bestpos) + '.consensus')
				cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 0 -l 1 %s > %s' % (sep_filename, consensus_filename)
				sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
				(status, output) = commands.getstatusoutput(cmd)

				# SPOA - MSA
				new_msa_filename = os.path.join(splitfolder, r_fname + '_' + base + ('_pos%d' % bestpos) + '.msa')
				cmd = '/home/kkrizanovic/src/spoa/build/bin/spoa -r 1 -l 1 %s > %s' % (sep_filename, new_msa_filename)
				sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
				(status, output) = commands.getstatusoutput(cmd)

				# Analyze msa
				new_varpos_filename = os.path.join(splitfolder, r_fname + '_' + base + ('_pos%d' % bestpos) + '.varpos')
				cmd = 'python %s %s %f > %s' % (analyze_msa, new_msa_filename, factor, new_varpos_filename)
				sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
				(status, output) = commands.getstatusoutput(cmd)

sys.stderr.write('\n\n') 
