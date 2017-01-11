#! /usr/bin/env python

# usage: python reroot_at_edge.py <tree_file> <head_node> <d2head> <out_file>

from Tree_extend import Tree_extend
from sys import argv
from os.path import splitext

tree_file = argv[1]
head = argv[2]
x = float(argv[3])
out_file = argv[4]
base_name,ext = splitext(tree_file)

a_tree = Tree_extend(tree_file=tree_file)
for edge in a_tree.ddpTree.preorder_edge_iter():
	#print edge.head_node.label
	if edge.head_node.label == head or (edge.head_node.is_leaf() and edge.head_node.taxon.label == head):
		a_tree.reroot_at_edge(edge,edge.length-x,x)
		#print("found!")
		break
#head_id, tail_id, edge_length, x = a_tree.Reroot()
#print("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x))

a_tree.tree_as_newick(outfile=out_file)


