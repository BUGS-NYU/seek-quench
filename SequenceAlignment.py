import argparse
import csv


def main():
	"""
	Intake command line parameters, import file, pairwise sequence alignment first column of each pair of rows, and export calculations.
	"""
	args = command_line_parameters()
	gap = int(args.gap)
	match = int(args.match)
	mismatch = int(args.mismatch)
	
	file1 = open(args.intake, 'r')
	content = list(csv.reader(file1))
	file1.close()
	# for every row of imported file, 
	# 	pairwise sequence align first two columns 
	# 	and export to next three columns
	for row in range(0, len(content), 2):
		seq1 = content[row][0]
		seq2 = content[row + 1][0]
		column = 1
		if args.global_:
			score_mat, end1, end2 = global_align(seq1, seq2, gap, match, mismatch)
			content, column = export(args, content, row, column, score_mat, end1, end2, 'Global Alignment')
		if args.semiglobal:
			score_mat, end1, end2 = semiglobal_align(seq1, seq2, gap, match, mismatch)
			content, column = export(args, content, row, column, score_mat, end1, end2, 'Semi Global Alignment')
		if args.local:
			score_mat, end1, end2 = local_align(seq1, seq2, gap, match, mismatch)
			content, column = export(args, content, row, column, score_mat, end1, end2, 'Local Alignment')
	file1 = open(args.intake, 'w')
	writer = csv.writer(file1)
	writer.writerows(content)
	file1.close()


def command_line_parameters():
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


def global_align(seq1, seq2, gap, match, mismatch):
	"""
	Compute global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = []
	for i in range(len(max(seq1,seq2))):
		if len(start_mat) < 1:
			start_mat.append(-1)
		else:
			start_mat.append(start_mat[-1] - 1)
	score_mat, path_mat = common_align(seq1, seq2, start_mat, gap, match, mismatch)
	end1 = ''
	end2 = ''
	row = len(seq2) - 1
	col = len(seq1) - 1
	while row >= 0 and col >= 0:
		if path_mat[row][col] == 'top':
			end1 = '-' + end1
			end2 = seq2[row] + end2
			row -= 1
		elif path_mat[row][col] == 'left':
			end1 = seq1[col] + end1
			end2 = '-' + end2
			col -= 1
		else:
			end1 = seq1[col] + end1
			end2 = seq2[row] + end2
			row -= 1
			col -= 1
	return score_mat, end1, end2


def semiglobal_align(seq1, seq2, gap, match, mismatch):
	"""
	Compute semi-global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = [0] * len(max(seq1, seq2))
	score_mat, path_mat = common_align(seq1, seq2, start_mat, gap, match, mismatch)
	end1 = ''
	end2 = ''
	row = 0
	col = len(seq1) - 1
	temp_score = score_mat[0][-1]
	for i in range(len(seq2)):
		if score_mat[i][-1] > temp_score:
			row = i
			col = len(seq1) - 1
			temp_score = score_mat[i][-1]
	for i in range(len(seq1)):
		if score_mat[-1][i] > temp_score:
			row = len(seq2) - 1
			col = i
			temp_score = score_mat[-1][i]
	while row >= 0 and col >= 0:
		if path_mat[row][col] == 'top':
			end1 = '-' + end1
			end2 = seq2[row] + end2
			row -= 1
		elif path_mat[row][col] == 'left':
			end1 = seq1[col] + end1
			end2 = '-' + end2
			col -= 1
		else:
			end1 = seq1[col] + end1
			end2 = seq2[row] + end2
			row -= 1
			col -= 1
	return score_mat, end1, end2


def local_align(seq1, seq2, gap, match, mismatch):
	"""
	Compute local sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = [0] * len(max(seq1, seq2))
	score_mat, path_mat = common_align(seq1, seq2, start_mat, gap, match, mismatch)
	end1 = ''
	end2 = ''
	row = 0
	col = 0
	temp_score = score_mat[0][0]
	for i in range(len(seq2)):
		for j in range(len(seq1)):
			if score_mat[i][j] > temp_score:
				temp_score = score_mat[i][j]
				row = i
				col = j
	while (row >= 0 and col >= 0) and temp_score > 0:
		temp_score = score_mat[row][col]
		if path_mat[row][col] == 'top':
			end1 = '-' + end1
			end2 = seq2[row] + end2
			row -= 1
		elif path_mat[row][col] == 'left':
			end1 = seq1[col] + end1
			end2 = '-' + end2
			col -= 1
		else:
			end1 = seq1[col] + end1
			end2 = seq2[row] + end2
			row -= 1
			col -= 1
	return score_mat, end1, end2


def common_align(seq1, seq2, start_mat, gap, match, mismatch):
	"""
	Calculates score matrix based on starting top row and left column scores created by unique alignment method.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:start_mat: starting numbers for computing top row / left column scoring matrix
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: competed score matrix and path matrix
	"""
	score_mat = []
	path_mat = []
	for row in range(len(seq2)):
		new_end = []
		new_path = []
		for col in range(len(seq1)):
			if col == 0 and row == 0:
				corner = 0
			elif row == 0:
				corner = start_mat[col - 1]
			elif col == 0:
				corner = start_mat[row - 1]
			else:
				corner = score_mat[row - 1][col - 1]
			if col == 0:
				left = start_mat[row] + gap
			else:
				left = new_end[col - 1] + gap
			if row == 0:
				top = start_mat[col] + gap
			else:
				top = score_mat[row - 1][col] + gap
			if seq1[col] == seq2[row]:
				corner += match
			else:
				corner += mismatch
			if corner > left and corner > top:
				new_path.append('corner')
			elif left > corner and left > top:
				new_path.append('left')
			else:
				new_path.append('top')
			new_end.append(max(corner, left, top))
		score_mat.append(new_end)
		path_mat.append(new_path)
	return score_mat, path_mat


def export(args, content, row, column, score_mat, end1, end2, method):
	"""
	Perform export functions depending on command line arguments.

	:param args: command line parameters for how to export and print
	:param content: imported file content
	:param row: current row number
	:param column: current column number
	:param score_mat: calculated score matrix
	:param end1: first finalized alignment
	:param end2: second finalized alignment
	:param method: alignment method
	"""
	if args.no_write:
		try:
			content[row][column] = end1
			content[row + 1][column] = end2
		except:
			content[row].append(end1)
			content[row + 1].append(end2)
		column += 1
	if len(args.export_matrix) > 0:
		score_str = '\n'.join('\t'.join(x for x in y) for y in score_mat)
		score_str += '\n'
		file2 = open(args.export_matrix, 'a')
		file2.write(score_str)
		file2.close()
	if args.print:
		print('> ' + method)
		print('\t{}\n\t{}'.format(end1, end2))
	return content, column


main()