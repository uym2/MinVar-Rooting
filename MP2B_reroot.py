#! /usr/bin/env python

# usage: python MP2B_reroot.py <tree_file>

from Tree_extend import MPR2B_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MPR2B_Tree(tree_file=tree_file)
head_id, tail_id, edge_length, x = a_tree.Reroot()
print("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x))
a_tree.tree_as_newick(outfile=base_name+"_MP2B_rooted"+ext)

