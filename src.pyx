#cython: language_level=3

import numpy as np


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
cdef int NONE = -1

DTYPE = np.intc

encoding = Bidict() # This dictionary holds the encoding from the data in seq to the data in the numpy array
encoding[-1] = '-'

#%%global_align.pyx

def global_align0(int [:] seq1, int len1, int [:] seq2, int len2, int gap, int match, int mismatch):
	"""
	Compute global sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""

	global TOP,LEFT,DTYPE, NONE
	cdef int length = max(len1,len2)
	start_mat = np.arange(-1,-length-1,-1,dtype = DTYPE)
	score_mat, path_mat = common_align0(seq1, len1, seq2, len2, start_mat, gap, match, mismatch)
	cdef int [:,:] score_mat_view = score_mat
	cdef int [:,:] path_mat_view = path_mat
	
	
	cdef Py_ssize_t col = len1 - 1
	cdef Py_ssize_t row = len2 - 1
	
	length = 0
	while row >= 0 and col >= 0:
		if path_mat_view[row][col] == TOP:
			length+=1
			row -= 1
		elif path_mat_view[row][col] == LEFT:
			length+=1
			col -= 1
		else:
			length+=1
			row -= 1
			col -= 1
	
	align1 = np.zeros(length, dtype = DTYPE);
	align2 = np.zeros(length, dtype = DTYPE);
	
	int [:] align1_view = align1
	int [:] align2_view = align2
	
	cdef Py_ssize_t index = length - 1
	col = len1 - 1
	row = len2 - 1
	
	while row >= 0 and col >= 0:
		if path_mat_view[row][col] == TOP:
			align1_vew[index] = NONE
			align2_view[index] = seq2[row]
			row -= 1
			index-=1
		elif path_mat_view[row][col] == LEFT:
			align1_view[index] = seq[col]
			align2_view[index] = NONE
			col -= 1
			index-=1
		else:
			align1_view[index] = seq1[col]
			align2_view[index] = seq2[row]
			row -= 1
			col -= 1
			index-=1
	
	return score_mat, align1, align2


#%%export.pyx
def export0(content,Py_ssize_t row,Py_ssize_t column,
			int [:,:] score_mat, end1, end2,char *method,int nowrite,char * export_matrix,int toprint):
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
	if nowrite != 0:
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
	if toprint != 0:
		print('> ' + method)
		print('\t{}\n\t{}'.format(end1, end2))
	return content, column


#%%seqarray.pyx

def encode(seq): # This converts a sequence into a numpy array
	cdef int code = 1, length = len(seq), index
	 # I'm almost certain that seq is just a string containing AGTC but im not confident
	cdef int [:] new_seq = np.zeros(length) # in that so I'm using a dictionary to encode the sequences first just in case
	global encoding
	cdef char *cseq = seq
	cdef char element
	encoding[cseq[0]] = 0
	for index in range(1,length):
		element = cseq[index]
		if element not in encoding:
			encoding[element] = code
			code+=1
		new_seq[index] = encoding[element]
	return new_seq

def decode(int [:] seqarray): # Convert numpy array back into seq string
	global encoding
	decoding = encoding.inverse
	cdef int length = len(seqarray),index
	cdef char [:] seq = np.zeros(length)
	for index in range(length):
		seq[index] = decoding[seqarray[index]]
	return ''.join(seq)


#%%common.pyx


def common_align0(int [:] seq1,int len1,int [:] seq2, int len2, int [:] start_matrix, int gap, int match, int mismatch):
	"""
	Calculates score matrix based on starting top row and left column scores created by unique alignment method.
	:param seq1: top row sequence
	:param seq2: left column sequence
	:start_mat: starting numbers for computing top row / left column scoring matrix
	
	:param gap: gap penalty THIS IS NEGATIVE
	
	:param match: match score THIS IS POSITIVE
	
	:param mismatch: mismatch score THIS IS NEGATIVE
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
	cdef int cols = len1, rows = len2

	matchBools = 1*np.equal(*np.meshgrid(seq1,seq2)) # 1 for match and 0 for not a match

	matchValues = matchBools * match - (matchBools - 1) * mismatch # Number to add when checking match
	
	cdef int [:,:] matchValues_view = matchValues
	
	# Matrices to hold data
	score_matrix = np.zeros(shape = (rows,cols),dtype = DTYPE); # Score matrix
	path_matrix = np.zeros(shape = (rows,cols),dtype = DTYPE); # path matrix

	cdef int[:,:] score_matrix_view = score_matrix
	cdef int[:,:] path_matrix_view = path_matrix

	cdef Py_ssize_t row,col
	cdef int corner,left,top,path,max_val
	
	path,max_val = getP(0+matchValues_view[0,0],start_matrix[0]+gap,start_matrix[0]+gap)#corner, left, top
	path_matrix_view[0,0] = path
	score_matrix_view[0,0] = max_val
	
	for row in range(1,rows):
		left = start_matrix[row] + gap
		top = score_matrix_view[row,0]
		corner = start_matrix[row-1] + matchValues_view[row,0]
		path,max_val = getP(corner,left,top)
		path_matrix_view[row,0]= path
		score_matrix_view[row,0] = max_val

	for col in range(1,cols):
		top = start_matrix[col] + gap
		left = score_matrix_view[0,col]
		corner = start_matrix[col-1] + matchValues_view[0,col]
		path,max_val = getP(corner,left,top)
		path_matrix_view[0,col]= path
		score_matrix_view[0,col] = max_val

	for row in range(1,rows):
		for col in range(1,cols):
			corner  = score_matrix_view[row,col]+matchValues_view[row,col]
			left = score_matrix_view[row,col-1]+gap
			top = score_matrix_view[row-1,col]+gap
			path,max_val = getP(corner,left,top)
			path_matrix_view[row,col] = path
			score_matrix_view[row,col] = max_val
	return score_matrix, path_matrix

def getP(int corner,int left, int top):
	global CORNER,LEFT,TOP
	if corner >= left and corner >= top:
		return CORNER,corner # corner
	elif left >= corner and left >= top:
		return LEFT,left # left
	else:
		return TOP,top # top
