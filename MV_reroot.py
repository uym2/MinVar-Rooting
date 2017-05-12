#! /usr/bin/env python

# usage: python MV_reroot.py <tree_file>

from Tree_extend import MV00_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MV00_Tree(tree_file=tree_file)
d2currRoot,br2currRoot = a_tree.Reroot()
print("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
a_tree.tree_as_newick(outfile=base_name+"_MV_rooted"+ext)
