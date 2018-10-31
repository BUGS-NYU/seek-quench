# seek-quench
Classic bioinformatic pairwise sequence alignment algorithms implemented into a simple open source web app that anyone can easily use. The project is named Seek Quench because it sounds like sequence and we are trying to *quench* answers that people are *seek*ing about their sequences.

Direct questions, comments, and concerns to [David](https://github.com/KEYS248).

## Long Term TODO
1. Build a web application using Django that takes in two nucleotide sequences, aligns them based on classic bioinformatic algorithms such as the Needleman–Wunsch algorithm for global sequence alignment, and outputs the results.
2. Allow more options!
	* different kinds of sequence alignments such as semi-global and local
	* multiple sequence alignment using Burrows–Wheeler transform
	* sequence conservation calculations such as sequence logos
	* phylogenetic tree construction
	* possibly other types of algorithms and options such as the profiles of PSI-Blast and Position Specific Scoring Matrices (PSSM)
3. Implement other algorithms such as k-means clustering, cluster coefficients, and UPGMA tree construction. Possibly into a new repo and just leave this repo for pairwise and multiple sequence alignments.

## Short Term TODO
1. Pairwise Sequence Alignment
	- [x] create algorithm functions
		- [x] global
		- [x] semi global
		- [x] local
		- [ ] have to figure out how to allow gaps at the front of the sequence
		- [ ] compute an overall alignment score (probably simply the highest number from score matrix)
	- [ ] more complex decisions to be made and implemented
		- [ ] decide if to run jobs as objects or simple function calls (currently latter)
		- [ ] be able to take two types of input: 
			- [ ] a CSV file of many sequences formatted correctly
			- [ ] two sequences passed to the program as variables
		- [ ] be able to output:
			- [ ] overall alignment score
			- [ ] final sequence alignment with gaps as dashes
			- [ ] score matrix with sequences on top and left hand column
			- [ ] same score matrix but exclude all scores not part of the final alignment
		- [ ] be able to output the results:
			- [ ] straight to the screen as a print
			- [ ] as an output file(s)
		- [ ] decide how best to format the input file:
			- each row is a separate sequence and the program will align every two sequential rows (e.g. align row 1 with 2, 3 with 4, etc)
			- each row is one pair, first sequence in column 0 and second sequence in column 1
			- some other option that's easy for people to understand
		- [ ] decide how best to output the results:
			- in the same input file, on the same rows as the sequences, just in the columns to the right
			- some other option that's easy for people to understand
	- [ ] recommendations on changes or other improvements are welcome
2. Django Implementation
	- [ ] I need to take the Django lynda course
	- [x] Create wireframe of front end
	- [ ] Implement wireframe with empty sockets for plugging functions into later
	- [ ] Incorporate algorithms into sockets for back end connection
	
