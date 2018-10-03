# -*- coding: utf-8 -*-
import numpy as np
from random import getrandbits

def check_equality(matrix1, matrix2, encoding = None):
	'''
	Checks to see if 2 matrices are equivalent. If an encoding is specified,
	first converts elements in matrix1 using encoding, then checks for equality
	with matrix2

	Parameters:
	matrix1 (array-like) - The first matrix to compare
	matrix2 (array-like) - The second matrix to compare
	encoding (dictionary, default = None) - A dictionary used to convert values in matrix1 into
	their equivalents in matrix2. If a value in matrix1 is not in the encoding, it will
	not be converted
	'''

	if encoding is None:
		return __check_naive(matrix1,matrix2)
	for row in range(len(matrix1)):
		for col in range(len(matrix1[0])):
			if (encoding[matrix1[row][col]] != matrix2[row][col]):
				print('indices: ({},{})'.format(row,col))
				print('Values: m1 - {}, m2 - {}'.format(encoding[matrix1[row][col]],matrix2[row][col]))
				return False
	return True

def __check_naive(matrix1,matrix2): # Helper function to checkEquality()
	for row in range(len(matrix1)):
		for col in range(len(matrix1[0])):
			if (matrix1[row][col] != matrix2[row][col]):
				print('indices: ({},{})'.format(row,col))
				print('Values: m1 - {}, m2 - {}'.format(matrix1[row][col],matrix2[row][col]))
				return False
	return True

VALUES = {0:'A',1:'G',2:'C',3:'T'}
BASES = ['A','G','T','C']
def generate_seq(length):
	"""
	Generates a sequence of given length using black magic
	Sequence is a string of AGCT
	"""
	global VALUES
	seq = []
	intseq = np.floor(np.random.rand(length) * 4)
	for i in range(length):

		seq.append(VALUES[intseq[i]])
	return ''.join(seq)

def mutate_seq(seq,mutation_rate = .01): # Given a sequence string, mutate it randomly according to a ratio
	# Use boolean mask to determine random mutation locations
	intseq = np.random.rand(len(seq)) < mutation_rate
	newseq = []
	global BASES
	for index in range(len(seq)):
		if intseq[index]:
			newseq.append(np.random.choice(BASES))
		else:
			newseq.append(seq[index])
	return ''.join(newseq)
