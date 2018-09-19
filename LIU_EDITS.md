# Liu Edits
Hi! I'm Albert Liu. Here's a summary of the edits I'm planning to make.

## Checklist
* [ ] .gitignore
* [ ] numPy
* [ ] Reorganizing functions

## Gitignore
I added a gitignore file so that the repository wouldn't have files related to my workflow but unrelated to the actual script.

## Integration of numPy
From what I saw, you were using python lists to store integers. I changed your sequence alignment functions to use numpy arrays instead of lists. Hopefully that'll increase performance by using vectorized operations and stuff like that.

## Re-Organization of Functions
I tried to move the argparse stuff away from everything else to make testing individual functions easier.

## Sources
Here's what I used to help me learn the bioinformatics and python packages I used
* [NeedlemanWunschPy by benhid](https://github.com/benhid/NeedlemanWunschPy/blob/master/NeedlemanWunschPy/algorithms.py) - This was what I used to get a sense of what the Needleman-Wunsch algorithm is.
* [SciPy Docs](https://docs.scipy.org/) - Just helpful information on how to use numPy.
* [Pairwise Sequence Alignment](https://towardsdatascience.com/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f) - Background info on pairwise sequence alignment
