import sys, os
import commands
import time
import re
import shutil

from datetime import datetime


SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, '/home/kkrizanovic/src/samscripts/src'))

results_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/results10'
final_results_folder = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/final_results10'
reference_file = '/home/kkrizanovic/jeleni/samo_jeleni_250-350/alela_jelen.fasta'

ref_basename = os.path.basename(reference_file)
updatable_reference = os.path.join(final_results_folder, ref_basename)

# Copy reference file to be updated
shutil.copy(reference_file, updatable_reference)

run_final = '/home/kkrizanovic/jeleni/run_final.py'

for diritem in os.listdir(results_folder):
	full_diritem = os.path.join(results_folder, diritem)
	if os.path.isdir(full_diritem):
		if full_diritem[-1] == '/':
			full_diritem = full_diritem[:-1]
		input_basename = os.path.basename(full_diritem)
		output_folder = os.path.join(final_results_folder, input_basename)
		if not os.path.exists(output_folder):
			sys.stderr.write('\nWorking on final results for %s!' % input_basename)
			os.mkdir(output_folder)
			cmd = 'python %s %s %s %s' % (run_final, full_diritem, output_folder, updatable_reference)
			sys.stderr.write('\nRUNNING COMMAND: %s' % cmd)
			(status, output) = commands.getstatusoutput(cmd)
		else:
			sys.stderr.write('\nWARNING: Final results folder %s already exists, skipping!' % input_basename)
			continue

sys.stderr.write('\n')
