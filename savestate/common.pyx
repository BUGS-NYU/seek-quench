import numpy as np


def common_align0( seq1, seq2, int [:] start_matrix, int gap, int match, int mismatch):
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
	cdef int cols = len(seq1), rows = len(seq2)
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
