import numpy as np
import argparse
import csv


# Global variables and definitions

# Copied from user Basj, linked in sources of LIU_EDITS.md
class Bidict(dict):
    def __init__(self, *args, **kwargs):
        super(Bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.items():
            self.inverse.setdefault(value,[]).append(key) 

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key) 
        super(Bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value,[]).append(key)        

    def __delitem__(self, key):
        self.inverse.setdefault(self[key],[]).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]: 
            del self.inverse[self[key]]
        super(Bidict, self).__delitem__(key)

TOP = 1
CORNER = 0
BOTTOM = -1

encoding = Bidict() # This dictionary holds the encoding from the data in seq to the data in the numpy array

# ---------seqarray.py--------------


def encode(seq): # This converts a sequence into a numpy array
	code = 1 # I'm almost certain that seq is just a string containing AGTC but im not confident
	new_seq = np.zeros(len(seq)) # in that so I'm using a dictionary to encode the sequences first just in case
	global encoding
	encoding[seq[0]] = 0
	for index in range(1,len(seq)):
		element = seq[index]
		if element not in seq:
			encoding[element] = code
			code+=1
		new_seq[index] = encoding[element]
	return new_seq

def decode(seqarray): # Convert numpy array back into seq string
	global encoding
	decoding = encoding.inverse
	seq = []
	for index in range(len(seqarray)):
		seq.append(decoding[seqarray[index]])
	return ''.join(seq)




# ---------global.py--------------

def global_align0(seq1, seq2, gap, match, mismatch):
	"""
	Compute global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = np.arange(-1,-max(len(seq1),len(seq2))-1,-1)
        
	score_mat, path_mat = common_align0(seq1, seq2, start_mat, gap, match, mismatch)
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


# ---------local.py--------------

def local_align0(seq1, seq2, gap, match, mismatch):
	"""
	Compute local sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = np.zeros(max(len(seq1),len(seq2)))
	score_mat, path_mat = common_align0(seq1, seq2, start_mat, gap, match, mismatch)
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

# ---------cmdlineparams.py--------------

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

# ---------export.py--------------

def export0(content, row, column, score_mat, end1, end2, method, nowrite, export_matrix, print):
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
	if nowrite:
		try:
			content[row][column] = end1
			content[row + 1][column] = end2
		except:
			content[row].append(end1)
			content[row + 1].append(end2)
		column += 1
	if len(export_matrix) > 0:
		score_str = '\n'.join('\t'.join(x for x in y) for y in score_mat)
		score_str += '\n'
		file2 = open(export_matrix, 'a')
		file2.write(score_str)
		file2.close()
	if print:
		print('> ' + method)
		print('\t{}\n\t{}'.format(end1, end2))
	return content, column



# ---------semiglobal.py--------------


def semiglobal_align0(seq1, seq2, gap, match, mismatch):
	"""
	Compute semi-global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = np.zeros(max(len(seq1),len(seq2)))
	score_mat, path_mat = common_align0(seq1, seq2, start_mat, gap, match, mismatch)
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

# ---------common.py--------------


def common_align0(seq1, seq2, start_mat, gap, match, mismatch):
	"""
	Calculates score matrix based on starting top row and left column scores created by unique alignment method.
	:param seq1: top row sequence
	:param seq2: left column sequence
	:start_mat: starting numbers for computing top row / left column scoring matrix
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: competed score matrix and path matrix

	Needleman Wunsch algorithm. Formula:
                M(0,0) = 0
                M(i,0) = M(i-1,0) - gap_penalty
                M(0,j) = M(0,j-1) - gap_penalty
                M(i,j) = max{ M(i-1,j-1) + score , M(i-1,j) + gap_penalty, M(i,j-1) + gap_penalty }
            M will be our scoring matrix.

	"""
	score_mat = []
	seq1 = list(np.round(4*np.random.random(500)-.5)) # Random sequence of 4 possible numbers
	seq2 = list(np.round(4*np.random.random(320))) # Random sequence of 4 possible numbers

	cols = len(seq1)
	rows = len(seq2)

	seq1mesh,seq2mesh  = np.meshgrid(seq1,seq2)
	matchBools = 1*np.equal(seq1mesh,seq2mesh) # 1 for match and 0 for not a match

	matchValues = matchBools * match - (matchBools - 1) * mismatch # Number to add when checking match


	score_matrix = np.zeros(shape = (rows+1,cols+1)); # Score matrix

	score_matrix[0,0] = 0 # Corner element M(0,0) = 0

	start_matrix = start_mat #+ np.linspace(0,len(start_mat)*gap,len(start_mat)+1)[0:len(start_mat)]
	score_matrix[1:,0] = start_matrix[0:rows] # First column values
	score_matrix[0,1:] = start_matrix[0:cols] # First row values

	for row in range(1,rows+1):
		for col in range(1,cols+1):
			score_matrix[row,col] = max(score_matrix[row-1,col-1]+matchValues[row-1,col-1],
									   score_matrix[row-1,col]+gap,
									   score_matrix[row,col-1]+gap)
	score_matrix = score_matrix[1:,1:]
	
	for row in range(len(seq2)):
		new_end = []
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
			new_end.append(max(corner, left, top))
		score_mat.append(new_end)

	return score_mat,score_matrix#, path_mat


def getP(left, corner, top):
	if corner > left and corner > top:
		return 0 # corner
	elif left > corner and left > top:
		return -1 # left
	else:
		return 1 # top


# ---------main.py--------------


def main0():
	"""
	Intake command line parameters, import file, pairwise sequence alignment first column of each pair of rows, and export calculations.
	"""
	args = command_line_parameters0() #Command call inputs
	gap = int(args.gap)
	match = int(args.match)
	mismatch = int(args.mismatch)


	file1 = open(args.intake, 'r') # 
	content = list(csv.reader(file1))
	file1.close()
	for row in range(0, len(content), 2):
		seq1 = content[row][0]
		seq2 = content[row + 1][0]
		column = 1
		if args.global_:
			score_mat, end1, end2 = global_align0(seq1, seq2, gap, match, mismatch)
			
			content, column = export(content, row, column, score_mat, end1, end2, 
							method = 'Global Alignment', nowrite = args.nowrite, 
							export_matrix = args.export_matrix, print = args.print)
		if args.semiglobal:
			score_mat, end1, end2 = semiglobal_align0(seq1, seq2, gap, match, mismatch)
			
			content, column = export(content, row, column, score_mat, end1, end2, 
							method = 'Semi Global Alignment', nowrite = args.nowrite, 
							export_matrix = args.export_matrix, print = args.print)
		if args.local:
			score_mat, end1, end2 = local_align0(seq1, seq2, gap, match, mismatch)
			
			content, column = export0(content, row, column, score_mat, end1, end2, 
							method = 'Local Alignment', nowrite = args.nowrite, 
							export_matrix = args.export_matrix, print = args.print)
	
	file1 = open(args.intake, 'w')
	writer = csv.writer(file1)
	writer.writerows(content)
	file1.close()

main0()