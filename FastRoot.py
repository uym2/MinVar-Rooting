#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree,MVR_Tree,MDR_Tree,MPR2_Tree
try:
	from dendropy4 import Tree,TreeList
except:
	from dendropy import Tree,TreeList

from os.path import splitext
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-m','--method',required=True,help="method: MP for midpoint, MV for minVAR, MD for mean-diff, MP2 for midpoint2")
parser.add_argument('-o','--outfile',required=False,help="specify output file. Default is stdout")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")
parser.add_argument('-f','--infofile',required=False,help="write info of the new root to file. This is important only for research and debugging purposes. Default is to write NOTHING.")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)
schema=args["schema"] if args["schema"] else "newick"

if args["outfile"]:
	outfile = args["outfile"]
else:
	outfile = None

if args["infofile"]:
	f_info = open(args["infofile"],'w')
else: 
	f_info = None

tree = Tree.get_from_path(tree_file,schema)	
if args["method"] == "MP":
	a_tree = MPR_Tree(ddpTree=tree)
elif args["method"] == "MD":
	a_tree = MDR_Tree(ddpTree=tree)
elif args["method"] == "MV":
	a_tree = MVR_Tree(ddpTree=tree)
elif args["method"] == "MP2":
	a_tree = MPR2_Tree(ddpTree=tree)
else:
	a_tree = None
if a_tree is None:
	print("Invalid method!")
else:
	#head_id, tail_id, edge_length, x = a_tree.Reroot()
	d2currRoot, br2currRoot = a_tree.Reroot()

if f_info:
	#f_info.write("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x)+ "\n")
	f_info.write("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
#else:	
#	print("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot))
a_tree.tree_as_newick(outfile=outfile)

if f_info:
	f_info.close()
