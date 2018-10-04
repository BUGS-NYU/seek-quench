from src import global_align0, export0, semiglobal_align0, local_align0
import argparse
import pandas as pd

import csv

def main0():
	"""
	Intake command line parameters, import file, pairwise sequence alignment first column of each pair of rows, and export calculations.
	"""
	args = command_line_parameters0() #Command call inputs
	gap = int(args.gap)
	match = int(args.match)
	mismatch = int(args.mismatch)
	seqdf = pd.read_csv(args.intake,header=None)

	# Eventually be able to label sequence pairs
	# Change export0 to not take in and return the object array every time
	
	# for every row of imported file,
	# 	pairwise sequence align first two columns
	# 	and export to next three columns
	for row in range(0, len(content), 2):
		seq1 = content[row][0]
		seq2 = content[row + 1][0]
		column = 1
		if args.global_:
			score_mat, end1, end2 = global_align0(seq1, seq2, gap, match, mismatch)

			content, column = export0(content, row, column, score_mat, end1, end2,
							method = 'Global Alignment', nowrite = args.nowrite,
							export_matrix = args.export_matrix, toprint = args.print)
		if args.semiglobal:
			score_mat, end1, end2 = semiglobal_align0(seq1, seq2, gap, match, mismatch)

			content, column = export0(content, row, column, score_mat, end1, end2,
							method = 'Semi Global Alignment', nowrite = args.nowrite,
							export_matrix = args.export_matrix, toprint = args.print)
		if args.local:
			score_mat, end1, end2 = local_align0(seq1, seq2, gap, match, mismatch)

			content, column = export0(content, row, column, score_mat, end1, end2,
							method = 'Local Alignment', nowrite = args.nowrite,
							export_matrix = args.export_matrix, toprint = args.print)

	file1 = open(args.intake, 'w')
	writer = csv.writer(file1)
	writer.writerows(content)
	file1.close()

def command_line_parameters0():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
				description="""
	Sequence Alignment Program
	---------------------------
	""")
	parser.add_argument('--global_', action='store_true', help='perform global sequence alignment')
	parser.add_argument('--semiglobal', action='store_true', help='perform semiglobal sequence alignment')
	parser.add_argument('--local', action='store_true', help='perform local sequence alignment')
	parser.add_argument('--gap', default='-1', help='specify gap penalty, defaults to -1')
	parser.add_argument('--match', default='1', help='specify match score, defaults to 1')
	parser.add_argument('--mismatch', default='0', help='specify mismatch score, defaults to 0')

	parser.add_argument('--print', action='store_true', help='print out alignments')
	parser.add_argument('--no_write', action='store_false', help='stop writing of new alignment to import file')
	parser.add_argument('--export_matrix', default='', help='export score matrix to specified file, defaults to no export')
	parser.add_argument('--intake', default='', help='import file with sequences to align')
	args = parser.parse_args()
	return args
