# seek-quench
Classic bioinformatic sequence alignment algorithms implemented into a simple open source web app. The project is named Seek Quench because it sounds like sequence and we are trying to *quench* answers that people are *seek*ing about their sequences.

## Goals
Projected to finish by the end of 2018.
1. Implement a web application using Django that takes in two nucleotide sequences, aligns them based on classic bioinformatic algorithms such as the Needleman–Wunsch algorithm for global sequence alignment, and outputs the results.
2. Allow more options!
	* different kinds of sequence alignments such as semi-global and local
	* multiple sequence alignment using Burrows–Wheeler transform
	* sequence conservation calculations such as sequence logos
	* phylogenetic tree construction
	* possibly other types of algorithms and options such as the profiles of PSI-Blast and Position Specific Scoring Matrices (PSSM)
3. Implement other algorithms such as k-means clustering, cluster coefficients, and UPGMA tree construction. Possibly into a new repo and just leave this repo for pairwise and multiple sequence alignments.

## TODO
1. Pairwise Sequence Alignment
	- [x] create algorithm functions
		- [x] global
		- [x] semi global
		- [x] local
		- [ ] have to figure out how to allow gaps at the front of the sequence
		- [ ] compute an overall alignment score (probably simply the highest number from score matrix)
	- [ ] decide how to best implement for later Django
		- [ ] e.g. objects or command line file calling with functions (currently latter)
		- [ ] take in file of many sequences formatted correctly and uploaded to front end or take in two sequences via command line from long response front end form box (currently former)
			- we may also want both to give users options for large and slow or small and quick alignments
		- [ ] these can be easily changed later but will take some time to ensure error-free algorithms are preserved
2. Django Implementation
	- [ ] I need to take the Django lynda course
	- [ ] Create wireframe of front end
	- [ ] Implement wireframe with empty sockets for plugging functions into later
	- [ ] Incorporate algorithms into sockets for back end connection
