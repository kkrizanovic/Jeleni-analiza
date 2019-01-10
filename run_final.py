import sys, os
import commands
import time
import re
import shutil
from datetime import datetime

spoa_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/spoa'

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, '/home/kkrizanovic/src/samscripts/src'))

from fastqparser import read_fastq

analyze_msa = '/home/kkrizanovic/jeleni/analyze_msa.py'
og_reference_file = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/alela_jelen.fasta'
minimap2 = '/home/kkrizanovic/src/minimap2/minimap2'

if len(sys.argv) < 3:
	sys.stderr.write('Run program: %s [INPUT folder] [OUTPUT folder] <updatable reference file>!\n' % sys.argv[0])
	exit(1)

input_folder = sys.argv[1]
if input_folder[-1] == '/':
	input_folder = input_folder[:-1]
output_folder = sys.argv[2]

start_adapters = ['ATTTCCTG', 'TTTTCCTG']
end_adapters   = ['CAGCGGCG', 'CAGCGGTC']

# If output folder doesn't exist, create it
if not os.path.exists(output_folder):
	os.mkdir(output_folder)

reference_file = og_reference_file
updateReference = False
if len(sys.argv) >= 4:
	reference_file = sys.argv[3]
	updateReference = True
	if not os.path.exists(reference_file):
		sys.stderr.write('\nERROR: given reference file (%s) deas not exist!\n' % reference_file)
		exit(1)

input_basename = os.path.basename(input_folder)
og_varpos_filename = os.path.join(input_folder, input_basename + '.varpos')
og_consensus_filename = os.path.join(spoa_folder, input_basename + '.consensus')
og_msa_filename = os.path.join(spoa_folder, input_basename + '.msa')
id_filename = os.path.join(output_folder, 'identification.txt')

### Load all consensuses into a list
consensus_list = []

# Used to check if a particular consensus is significant enough (build using a significant subset of the original)
maxsize = os.path.getsize(og_msa_filename)
factor = 0.10


# 1. Read consensuses within subfolders within the input folder
for diritem1 in os.listdir(input_folder):
	full_diritem1 = os.path.join(input_folder, diritem1)
	if os.path.isdir(full_diritem1):
		for diritem2 in os.listdir(full_diritem1):
			cfilename = os.path.join(full_diritem1, diritem2)
			if cfilename.lower().endswith('.consensus'):
				msa_filename = cfilename[:-len('consensus')] + 'msa'
				size = os.path.getsize(msa_filename)
				if size < factor*maxsize:
					sys.stderr.write('\nWARINIG: Consensus %s discarded because of low significance!' % diritem2)
					continue

				with open(cfilename, 'r') as cfile:
					cseq = cfile.readline()
					cseq = cfile.readline()[:-1]		# Sequence is in the second line, also removing '\n' at the end
					if cseq not in consensus_list:
						consensus_list.append(cseq)

# 2. Read split consensuses within the input folder
for diritem in os.listdir(input_folder):
	cfilename = os.path.join(input_folder, diritem)
	if cfilename.lower().endswith('.consensus'):
		msa_filename = cfilename[:-len('consensus')] + 'msa'
		size = os.path.getsize(msa_filename)
		if size < factor*maxsize:
			sys.stderr.write('\nWARNINIG: Consensus %s discarded because of low significance!' % diritem)
			continue

		with open(cfilename, 'r') as cfile:
			cseq = cfile.readline()
			cseq = cfile.readline()[:-1]		# Sequence is in the second line, also removing '\n' at the end
			if cseq not in consensus_list:
				consensus_list.append(cseq)

# 3. load original consensus
with open(og_consensus_filename, 'r') as cfile:
	cseq = cfile.readline()
	cseq = cfile.readline()[:-1]		# Sequence is in the second line, also removing '\n' at the end
	if cseq not in consensus_list:
		consensus_list.append(cseq)


min_trim_len = 245
max_trim_len = 255
### Trim adapters from consensuses
trimmed_consensus_list = []
for consensus in consensus_list:
	trimmed_consensus = consensus
	lpos = -1
	for adapter in start_adapters:
		lpos = consensus.find(adapter)
		if lpos > -1:
			trimmed_consensus = trimmed_consensus[lpos+len(adapter):]
			break

	rpos = -1
	for adapter in end_adapters:
		rpos = trimmed_consensus.rfind(adapter)
		if rpos > -1:
			trimmed_consensus = trimmed_consensus[:rpos]
			break

	if lpos == -1 or rpos == -1:
		sys.stderr.write('\nWARNING: could not find adapter in %s' % consensus)

	# Check trimmed consensus for length
	if len(trimmed_consensus) > max_trim_len or len(trimmed_consensus) < min_trim_len:
		sys.stderr.write('\nWARNING: trimmed consensus outside of expected size(%d)!' % len(trimmed_consensus))

	trimmed_consensus_list.append(trimmed_consensus)

### Check if any of the trimmed consensi are equal
new_consensus_list = []
new_trimmed_consensus_list = []
for i in xrange(len(trimmed_consensus_list)):
	consensus = consensus_list[i]
	trimmed_consensus = trimmed_consensus_list[i]

	if trimmed_consensus not in new_trimmed_consensus_list:
		new_trimmed_consensus_list.append(trimmed_consensus)
		new_consensus_list.append(consensus)
	else:
		sys.stderr.write('\nWARNING: found equal consensi after trimming, dropping one!')

consensus_list = new_consensus_list
trimmed_consensus_list = new_trimmed_consensus_list

# Keep only the first 4 consensi
consensus_list = consensus_list[:4]
trimmed_consensus_list = trimmed_consensus_list[:4]

### Try to identify consensuses by exact matching with references
### Writing consensusi and identification to a file
### Running minimap against reference for each untrimmed and trimmed consensus
# Loading references
[headers, seqs, quals] = read_fastq(reference_file)
id_file = open(id_filename, 'w')

# Checking if references already contain new sequences, to determine 
# Accurate serial nubmber of each new sequence

newSeqs = False
newSeqNo = 0


# If new sequences were already detected, determin the number of the last sequence
if headers[-1].find('New_sequence') != -1:
	pos = headers[-1].rfind('_')
	newSeqNo = int(headers[-1][pos+1:])

elog_filename = os.path.join(output_folder, 'error_log.txt')

for index in xrange(len(consensus_list)):
	consensus = consensus_list[index]
	trimmed_consensus = trimmed_consensus_list[index]
	identified = False
	identity = 'Unidentified'

	if trimmed_consensus == '':
		elog_file = open(elog_filename, 'a')
		elog_file.write('ERROR: empty trimmed consensus for (%s)!\n' % consensus)
		elog_file.close()
		# import pdb
		# pdb.set_trace()
		continue


	# Write trimmed and untrimmed consensus to a file
	cfilename = os.path.join(output_folder, input_basename + '_consensus%d_untrimmed.fasta' % (index+1))
	tcfilename = os.path.join(output_folder, input_basename + '_consensus%d_trimmed.fasta' % (index+1))

	for i in xrange(len(headers)):
		header = headers[i]
		seq = seqs[i]
		if trimmed_consensus == seq:
			identified = True
			identity = header
			break

	if not identified:
		newSeqNo += 1
		identity = 'New_sequence_%0d' % (newSeqNo)
		newSeqs = True
		headers.append(identity)
		seqs.append(trimmed_consensus)

	id_file.write('Consensus%d - %s\n' % (index+1, identity))

	with open(cfilename, 'w') as cfile:
		cfile.write('>Consensus%d - %s\n%s\n' % (index+1, identity, consensus))
	with open(tcfilename, 'w') as tcfile:
		tcfile.write('>Consensus%d - %s\n%s\n' % (index+1, identity, trimmed_consensus))


	# Run minimap2
	untrimmed_mm2_output = os.path.join(output_folder, input_basename + '_consensus%d_umapping.sam' % (index+1))
	cmd = '%s --eqx -ax sr %s %s > %s' % (minimap2, cfilename, og_reference_file, untrimmed_mm2_output)
	sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
	(status, output) = commands.getstatusoutput(cmd)

	trimmed_mm2_output = os.path.join(output_folder, input_basename + '_consensus%d_tmapping.sam' % (index+1))
	cmd = '%s --eqx -ax sr %s %s > %s' % (minimap2, tcfilename, og_reference_file, trimmed_mm2_output)
	sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
	(status, output) = commands.getstatusoutput(cmd)


# Collect trimmed consensi into a single file
consensus_summary_file = os.path.join(output_folder, input_basename + '_consensus_all.fasta')
cmd = 'cat %s/*_trimmed.fasta > %s' % (output_folder, consensus_summary_file)
sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
(status, output) = commands.getstatusoutput(cmd)

# Write down all discovered sequences in a separate file
if updateReference and newSeqs:
	with open(reference_file, 'w') as summary_file:
		for i in xrange(len(headers)):
			header = headers[i]
			seq = seqs[i]
			summary_file.write('>%s\n%s\n' % (header, seq))

id_file.close()
sys.stderr.write('\n')



