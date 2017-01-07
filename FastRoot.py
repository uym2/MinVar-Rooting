#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree,minVAR_Tree,MDR_Tree
try:
	from dendropy4 import Tree,TreeList
except:
	from dendropy import Tree,TreeList

from sys import argv
from os.path import splitext
import argparse
import optparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-m','--method',required=True,help="method: MP for midpoint, MV for minVAR, MD for mean-diff")
parser.add_argument('-o','--outfile',required=False,help="specify output file")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")
parser.add_argument('-f','--infofile',required=False,help="write info of the new root to file. This is important only for research and debugging purposes. Default is to write NOTHING.")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)
schema=args["schema"] if args["schema"] else "newick"
if args["outfile"]:
	outfile = args["outfile"]
else:
	outfile = base_name + "_" + args["method"] + "_rooted" + ext
try:
	os.remove(outfile)
except:
	pass

try:
	os.remove(args["infofile"])
except:
	pass

i = 1

#trees = TreeList.get(path=tree_file,schema=schema)
#fin = open(tree_file,'r')
#for tree in trees:
tree = Tree.get(path=tree_file,schema=schema)
if args["method"] == "MP":
	a_tree = MPR_Tree(ddpTree=tree)
elif args["method"] == "MD":
	a_tree = MDR_Tree(ddpTree=tree)
elif args["method"] == "MV":
	a_tree = minVAR_Tree(ddpTree=tree)
else:
	a_tree = None

if a_tree is None:
	print "Invalid method!"
else:
	head_id, tail_id, edge_length, x = a_tree.Reroot()

	if args["infofile"]:
		with open(args["infofile"],'w') as f:
			f.write("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x) + "\n")
	a_tree.tree_as_newick(outfile=outfile)
