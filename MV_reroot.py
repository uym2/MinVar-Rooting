#! /usr/bin/env python

# usage: python MV_reroot.py <tree_file>

from Tree_extend import MVR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MVR_Tree(tree_file=tree_file)
#head_id, tail_id, edge_length, x = a_tree.Reroot()
d2currRoot,br2currRoot = a_tree.Reroot()
#print("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x))
print("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
a_tree.tree_as_newick(outfile=base_name+"_MV_rooted"+ext)


