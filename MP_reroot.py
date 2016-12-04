# usage: python MP_reroot.py <tree_file>

from Tree_extend import MPR_Tree
from sys import argv

tree_file = argv[1]

a_tree = MPR_Tree(tree_file=tree_file)
a_tree.Reroot()
a_tree.tree_as_newick(outfile=tree_file.split(".tre")[0]+"_MP_rooted.tre",restore_label=True)


