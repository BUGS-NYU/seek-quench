import numpy as np

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
