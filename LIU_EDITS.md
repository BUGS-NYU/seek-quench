# Liu Edits
Hi! I'm Albert Liu. Here's a summary of the edits I'm planning to make/have made.

## Checklist
* [ ] numpy
* [ ] cython
* [ ] pandas
* [ ] README TODO's

## Overview
In general the edits that I made were for method specificity, runtime, and bug fixing.

## NumPy
I used numpy to make operations faster.
* [ ] Alignment Functions
	* [ ] global
	* [ ] semi-global
	* [ ] local
	* [ ] common

# Pandas
I used pandas to make reading in data simpler and faster.
* [ ] importing
* [ ] exporting

## Cython
I used cython to reduce the python/numpy overhead with accessing arrays.
* [ ] Variables
	* [ ] global_align
	* [ ] semi-global_align
	* [ ] local_align
	* [ ] common_align
* [ ] Arrays
	* [ ] global_align
	* [ ] semi-global_align
	* [ ] local_align
	* [ ] common_align

## Sources and Acknowledgements
Here's what I used to help me learn the bioinformatics and python packages I used
* [NeedlemanWunschPy by benhid](https://github.com/benhid/NeedlemanWunschPy/blob/master/NeedlemanWunschPy/algorithms.py) - This was what I used to get a sense of what the Needleman-Wunsch algorithm is.
* [SciPy Docs](https://docs.scipy.org/) - Just helpful information on how to use numPy.
* [Pairwise Sequence Alignment](https://towardsdatascience.com/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f) - Background info on pairwise sequence alignment
* [Distances between Aligned Sequences](https://www.inf.ethz.ch/personal/gonnet/papers/Distance/Distance.html) - Got a general introduction to the problem from this paper by Gaston Gonnet and Chantal Korostensky. (I didn't read much of it but it still helped some.)
* [PEP 305](https://www.python.org/dev/peps/pep-0305/#reading-csv-files) - This PEP describes the general concept of the csv package
* [pandas](https://github.com/pandas-dev/pandas) - I don't know yet, but I'm probably going to replace the csv stuff with pandas.
* [Bidirectional Hash Table](https://stackoverflow.com/questions/3318625/efficient-bidirectional-hash-table-in-python) - Used this for encoder
* Denis Kaydanov - My CL helped me understand the underlying biology better
* [Cython for Numpy Users](https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html) - Introduction for Cython
  * [Kurt Smith's Cython Tutorial](https://www.youtube.com/watch?v=gMvkiQ-gOW8&t=4730s&ab_channel=Enthought) - Lecture at SciPy 2015
  * [Cython Tutorial Slides](https://github.com/kwmsmith/scipy-2015-cython-tutorial)
* [BioPython](https://biopython.org/wiki/Documentation) - Biopython is an option that I haven't yet explored
* [Recursive Vectorization with SciPy](https://stackoverflow.com/questions/21336794/python-recursive-vectorization-with-timeseries/21338665#21338665) - SciPy is an option I haven't yet explored
