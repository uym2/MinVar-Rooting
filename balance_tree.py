#! /usr/bin/env python

from Tree_extend import MBR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MBR_Tree(tree_file=tree_file)
bt = a_tree.build_balance_tree()
bt.write_to_path(base_name+"_balance"+ext,"newick")
#a_tree.list_balance_points()
#head_id, tail_id, edge_length, x = a_tree.Reroot()
#print("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x))
#a_tree.tree_as_newick(outfile=base_name+"_MD_rooted"+ext)

