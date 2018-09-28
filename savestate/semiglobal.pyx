import numpy as np

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
