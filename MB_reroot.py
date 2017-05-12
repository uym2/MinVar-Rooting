#! /usr/bin/env python

# usage: python MB_reroot.py <tree_file>

from Tree_extend import MBR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MBR_Tree(tree_file=tree_file)

d2currRoot,br2currRoot = a_tree.Reroot()
print("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
a_tree.tree_as_newick(outfile=base_name+"_MB_rooted"+ext)


