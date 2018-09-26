#cython: language_level=3

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

cdef int CORNER = 0
cdef int LEFT = 1
cdef int TOP = 2

DTYPE = np.intc

encoding = Bidict() # This dictionary holds the encoding from the data in seq to the data in the numpy array

#%%main.pyx
# -*- coding: utf-8 -*-


def main0():
	"""
	Intake command line parameters, import file, pairwise sequence alignment first column of each pair of rows, and export calculations.
	"""
	args = command_line_parameters0() #Command call inputs
	cdef int gap = int(args.gap)
	cdef int match = int(args.match)
	cdef int mismatch = int(args.mismatch)

	file1 = open(args.intake, 'r') #
	content = list(csv.reader(file1))
	file1.close()
	# for every row of imported file,
	# 	pairwise sequence align first two columns
	# 	and export to next three columns
	cdef Py_ssize_t row,column
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


#%%global.pyx
# -*- coding: utf-8 -*-

def global_align0(seq1, seq2, int gap, int match,int mismatch):
	"""
	Compute global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""

	global TOP,LEFT,DTYPE
	cdef int length = max(len(seq1),len(seq2))
#	start_mat = np.array([x for x in range(-1,-length-1,-1)],dtype = DTYPE)
	start_mat = np.arange(-1,-length-1,-1,dtype = DTYPE)
	score_mat, path_mat = common_align0(seq1, seq2, start_mat, gap, match, mismatch)
	cdef int [:,:] score_mat_view = score_mat
	cdef int [:,:] path_mat_view = path_mat
	end1 = []
	end2 = []
	cdef char *cseq1 = seq1
	cdef char *cseq2 = seq2
	cdef Py_ssize_t row = len(seq2) - 1
	cdef Py_ssize_t col = len(seq1) - 1
	while row >= 0 and col >= 0:
		if path_mat[row][col] == TOP:
			end1.append('-')
			end2.append(cseq2[row])
			row -= 1
		elif path_mat[row][col] == LEFT:
			end1.append(cseq1[col])
			end2.append('-')
			col -= 1
		else:
			end1.append(cseq1[col])
			end2.append(cseq2[row])
			row -= 1
			col -= 1
	end1.reverse()
	end2.reverse()
	return score_mat, ''.join(end1), ''.join(end2)

#%%cmdlineparams.pyx
# -*- coding: utf-8 -*-

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

#%%export.pyx
# -*- coding: utf-8 -*-

def export0(content,Py_ssize_t row,Py_ssize_t column,
			int [:,:] score_mat, end1, end2,char *method, nowrite, export_matrix, toprint):
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
	if toprint:
		print('> ' + method)
		print('\t{}\n\t{}'.format(end1, end2))
	return content, column


#%%local.pyx
# -*- coding: utf-8 -*-

def local_align0(seq1, seq2, int gap,int match,int mismatch):
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
	score_mat, path_mat = common_align0(seq1, seq2,start_mat,gap,match,mismatch)
	end1 = []
	end2 = []
	cdef int row = 0
	cdef int col = 0
	cdef int temp_score = score_mat[0,0]
	cdef char *cseq1 = seq1,
	cdef char *cseq2 = seq2
	global TOP,LEFT
	cdef int i,j
	for i in range(len(seq2)):
		for j in range(len(seq1)):
			if score_mat[i][j] > temp_score:
				temp_score = score_mat[i][j]
				row = i
				col = j
	while (row >= 0 and col >= 0) and temp_score > 0:
		temp_score = score_mat[row][col]
		if path_mat[row][col] == TOP:
			end1.append('-')
			end2.append(cseq2[row])
			row -= 1
		elif path_mat[row][col] == LEFT:
			end1.append(cseq1[col])
			end2.append('-')
			col -= 1
		else:
			end1.append(cseq1[col])
			end2.append(cseq2[row])
			row -= 1
			col -= 1
	end1.reverse()
	end2.reverse()
	return score_mat, ''.join(end1), ''.join(end2)

#%%seqarray.pyx
# -*- coding: utf-8 -*-


def encode(seq): # This converts a sequence into a numpy array
	code = 1 # I'm almost certain that seq is just a string containing AGTC but im not confident
	new_seq = np.zeros(len(seq)) # in that so I'm using a dictionary to encode the sequences first just in case
	global encoding
	cdef char *cseq = seq
	encoding[cseq[0]] = 0
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




#%%semiglobal.pyx
# -*- coding: utf-8 -*-


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
	end1 = []
	end2 = []
	cdef char *cseq1 = seq1
	cdef char *cseq2 = seq2
	cdef int row = 0
	cdef int col = len(seq1) - 1
	cdef int temp_score = score_mat[0][-1]
	global TOP,LEFT
	cdef int i
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
		if path_mat[row][col] == TOP:
			end1.append('-')
			end2.append(cseq2[row])
			row -= 1
		elif path_mat[row][col] == LEFT:
			end1.append(cseq1[col])
			end2.append('-')
			col -= 1
		else:
			end1.append(cseq1[col])
			end2.append(cseq2[row])
			row -= 1
			col -= 1
	end1.reverse()
	end2.reverse()
	return score_mat, ''.join(end1), ''.join(end2)


#%%common.pyx
# -*- coding: utf-8 -*-


def common_align0(seq1, seq2,int [:] start_matrix, int gap, int match, int mismatch):
	"""
	Calculates score matrix based on starting top row and left column scores created by unique alignment method.
	:param seq1: top row sequence
	:param seq2: left column sequence
	:start_mat: starting numbers for computing top row / left column scoring matrix
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: competed score matrix and path matrix

	# Copied from benhid
	Needleman Wunsch algorithm. Formula:
                M(0,0) = 0
                M(i,0) = M(i-1,0) - gap_penalty
                M(0,j) = M(0,j-1) - gap_penalty
                M(i,j) = max{ M(i-1,j-1) + score , M(i-1,j) + gap_penalty, M(i,j-1) + gap_penalty }
            M will be our scoring matrix.

	"""

	# Number of rows and columns
	global DTYPE
	cdef int cols, rows
	cols,rows = len(seq1),len(seq2)

	# Combine these statements into 1 eventually
	matchBools = 1*np.equal(*np.meshgrid(seq1,seq2)) # 1 for match and 0 for not a match

	matchValues = matchBools * match - (matchBools - 1) * mismatch # Number to add when checking match
	cdef int [:,:] matchValues_view = matchValues
	# Matrices to hold data
	score_matrix = np.zeros(shape = (rows+1,cols+1),dtype = DTYPE); # Score matrix
	path_matrix = np.zeros(shape = (rows,cols),dtype = DTYPE); # path matrix

	cdef int[:,:] score_matrix_view = score_matrix
	cdef int[:,:] path_matrix_view = path_matrix
	score_matrix[0,0] = 0 # Corner element M(0,0) = 0
	
	score_matrix_view[1:,0] = start_matrix[0:rows] # First column values
	score_matrix_view[0,1:] = start_matrix[0:cols] # First row values
	cdef Py_ssize_t row,col
	cdef int corner,left,top,path,max_val
	for row in range(1,rows+1):
		for col in range(1,cols+1):
			corner  = score_matrix_view[row-1,col-1]+matchValues_view[row-1,col-1]
			left = score_matrix_view[row,col-1]+gap
			top = score_matrix_view[row-1,col]+gap
			path,max_val = getP(corner,left,top)
			path_matrix_view[row-1,col-1] = path
			score_matrix_view[row,col] = max_val
	score_matrix = score_matrix[1:,1:]
	return score_matrix, path_matrix

def getP(int corner,int left, int top):
	global CORNER,LEFT,TOP
	if corner >= left and corner >= top:
		return CORNER,corner # corner
	elif left >= corner and left >= top:
		return LEFT,left # left
	else:
		return TOP,top # top