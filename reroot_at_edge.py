#! /usr/bin/env python

# usage: python reroot_at_edge.py <tree_file> <head_node> <d2head> <out_file>

from Tree_extend import Tree_extend
#from sys import argv
from os.path import splitext

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-n','--hnode',required=True,help="the label of the head (child) node)")
parser.add_argument('-d','--distance',required=True,help="distance of the new root to head (child) node)")
parser.add_argument('-o','--outfile',required=False,help="specify output file")

args = vars(parser.parse_args())

tree_file = args["input"]
head = args["hnode"]
x = float(args["distance"])

if args["outfile"]:
	out_file = args["outfile"]
else:
	out_file = None

base_name,ext = splitext(tree_file)

a_tree = Tree_extend(tree_file=tree_file)
for edge in a_tree.ddpTree.preorder_edge_iter():
	if (edge.head_node.label == head) or (edge.head_node.is_leaf() and edge.head_node.taxon.label == head):
		if (edge is not None) and edge.length:
			a_tree.reroot_at_edge(edge,edge.length-x,x)
		break

a_tree.tree_as_newick(outfile=out_file)


