# -*- coding: utf-8 -*-

def checkEquality(matrix1, matrix2, encoding = None):
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
	
	#TODO Do this function to test the changes made
	if encoding is None:
		return checkEqualityNaive(matrix1,matrix2)
	for row in range(len(matrix1)):
		for col in range(len(matrix1[0])):
			if (encoding[matrix1[row][col]] != matrix2[row][col]):
				print(row,col)
				return False
	return True

def checkEqualityNaive(matrix1,matrix2):
	for row in range(len(matrix1)):
		for col in range(len(matrix1[0])):
			if (matrix1[row][col] != matrix2[row][col]):
				print(row,col)
				return False
	return True