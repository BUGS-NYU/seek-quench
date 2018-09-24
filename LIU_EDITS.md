# Liu Edits
Hi! I'm Albert Liu. Here's a summary of the edits I'm planning to make/have made.

## Checklist
* [x] bugs
* [x] .gitignore
* [x] Reorganizing functions
* [x] Code Simplicity and Efficiency

## Overview
In general the edits that I made were for method specificity, runtime, and bug fixing.

## Bugs
I tried to fix any bugs that I found. I don't really know what exactly I should be looking for, but I think I fixed a bug with the generation of the starting matrix, where the matrix would have an incorrect length.

## Gitignore
I added a gitignore file so that the repository wouldn't have files related to my workflow but unrelated to the actual script.

## Re-Organization of Functions
I tried to move the argparse stuff away from everything else to make testing individual functions easier.

## Code simplification and Efficiency
Changed a few things to make code simpler, i.e. used list comprehension instead of for loop. Also improved efficiency by using lists instead of strings in global, semi-global, and local alignment functions.

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
