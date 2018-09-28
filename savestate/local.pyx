import numpy as np

def local_align0(char *seq1, int len1,char *seq2, int len2, int gap, int match, int mismatch):
	"""
	Compute local sequence alignment on pair.

	:param seq1: top row sequence
	:param seq2: left column sequence
	:param gap: gap penalty
	:param match: match score
	:param mismatch: mismatch score
	:return: completed score matrix and final sequence alignments
	"""
	start_mat = np.zeros(max(len1,len2))
	score_mat, path_mat = common_align0(seq1,len1,seq2,len2,start_mat,gap,match,mismatch)
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
