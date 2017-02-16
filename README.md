This is a dendropy-based implementation for a new idea of rooting a phylogenetic tree, called the minVar (MV) root. Our method roots the tree at the point that minimizes the variance of the root to tip distances. The code was designed to be easily extensible for a class of 'optimization-based' rooting methods, which includes a fast version of midpoint (MP) rooting as well.

Complexity: the rooting methods are linear (with the number of species) in time and space.

Dependencies:
- python (version 2.7 recommended)
- dendropy (version 4 recommended)

Usage:

*
FastRoot.py [-h] -i INPUT -m METHOD [-o OUTFILE] [-s SCHEMA]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file
  -m METHOD, --method METHOD
                        method: MP for midpoint and MV for minVAR
  -o OUTFILE, --outfile OUTFILE
                        specify output file
  -s SCHEMA, --schema SCHEMA
                        schema of your input treefile. Default is newick

NOTE: FastRoot.py works with tree list. It will reroot all trees in the treefile by the method of your choice. 

Besides, the following scripts might seem easier to use, but only work with single-tree input, and don't accept optional options.

*
python MP_reroot.py \<tree_file_in_newick\>

*
python minVAR_reroot.py \<tree_file_in_newick\>

Output:
If you called FastRoot.py with -o, you know what you named the output. Otherwise, in the same directory as the input tree you will find the (re)rooted tree with the same name as the input with suffix MP_rooted or MV_rooted.
