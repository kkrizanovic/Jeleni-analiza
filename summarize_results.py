import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

spoa_folder = '/home/kkrizanovic/jeleni/samo_jeleni/spoa'

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, '/home/kkrizanovic/src/samscripts/src'))

from fastqparser import read_fastq

if len(sys.argv) != 5:
	sys.stderr.write('Run program: %s [INPUT folder] [Consensus fasta file] [ouput file 1] [output file 2]!\n' % sys.argv[0])
	exit(1)

input_folder = sys.argv[1]
if input_folder[-1] == '/':
	input_folder = input_folder[:-1]
consensus_fasta_file = sys.argv[2]
output_file1 = sys.argv[3]
output_file2 = sys.argv[4]

# Openning consensus sequences and placing them in a dictionary
[cheaders, cseqs, cquals] = read_fastq(consensus_fasta_file)
consensus_dict = {}
for i in xrange(len(cheaders)):
	cheader = cheaders[i]
	cseq = cseqs[i]
	consensus_dict[cheader] = cseq

# A dictionary containing summary information
# Key is a reference name, value is a list of samples containing that reference
summary_dict1 = {}

# A dictionary containing summary information
# Key is a sample name, value is a list of sequences contained in that sample
summary_dict2 = {}

for diritem in os.listdir( input_folder):
	full_diritem = os.path.join(input_folder, diritem)
	if os.path.isdir(full_diritem):
		sample_cfilename = os.path.join(full_diritem, diritem + '_consensus_all.fasta')

		[di_headers, di_seqs, di_quals] = read_fastq(sample_cfilename)
		for i in xrange(len(di_headers)):
			di_header = di_headers[i]
			header = di_header[len('Consensus1 - '):]
			di_seq = di_seqs[i]

			if header not in consensus_dict:
				sys.stderr.write('\nERROR: header %s for sample %s not in references!' % (di_header, diritem))

			cseq = consensus_dict[header]
			if cseq != di_seq:
				sys.stderr.write('\nERROR: sequence mismatch for header %s in sample %s!' % (di_header, diritem))

			short_header = header
			pos = header.find(' ')
			if pos != -1:
				short_header = header[:pos]

			if short_header not in summary_dict1:
				summary_dict1[short_header] = [diritem]
			else:
				summary_dict1[short_header].append(diritem)

			if diritem not in summary_dict2:
				summary_dict2[diritem] = [short_header]
			else:
				summary_dict2[diritem].append(short_header)

# Write summary dictionary 1 to output file 1
with open(output_file1, 'w') as ofile1:
	ofile1.write('SEQUENCE NAME; SAMPLE1, SAMPLE2, ...\n')
	for cname in sorted(summary_dict1.keys()):
		sample_list = summary_dict1[cname]
		ofile1.write('%s: %s\n' % (cname, ', '.join(sample_list)))

# Write summary dictionary 2 to output file 2
with open(output_file2, 'w') as ofile2:
	ofile2.write('SAMPLE NAME; SEQUENCE1, SEQUENCE2, ...\n')
	for sname in sorted(summary_dict2.keys()):
		seq_list = summary_dict2[sname]
		ofile2.write('%s: %s\n' % (sname, ', '.join(seq_list)))
