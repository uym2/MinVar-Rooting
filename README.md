# MP_reroot
This is a fast implementation of the midpoint rerooting. The code was designed to be reusable for implementation of a new idea of rerooting, called the minVar reroot, which reroots the tree at the point that minimizes variances of the root to tip distances.

Complexity: both rerooting methods are linear (with the number of tips) in time and memory.


Installation:
Work out of the box. All you need is python 2.7

Usage:
python MP_reroot.py <tree_file_in_newick>
python minVAR_reroot.py <tree_file_in_newick>

Output:
In the same directory as the input tree, you will see the output tree with the same name as the input added suffix MP_rooted or VAR_rooted
