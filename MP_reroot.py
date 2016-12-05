# usage: python MP_reroot.py <tree_file>

from Tree_extend import MPR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MPR_Tree(tree_file=tree_file)
a_tree.Reroot()
a_tree.tree_as_newick(outfile=base_name+"_MP_rooted"+ext,restore_label=True)


