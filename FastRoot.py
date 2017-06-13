#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree,MV00_Tree,MBR_Tree,MV0F_Tree,MVD0_Tree,MVDF_Tree
from dendropy import Tree,TreeList

from os.path import splitext
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-m','--method',required=True,help="method: MP for midpoint, MV for minVAR")
parser.add_argument('-o','--outfile',required=True,help="specify output file")
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

try:
	os.remove(outfile)
except:
	pass

try:
	os.remove(args["infofile"])
except:
	pass

if args["infofile"]:
	f_info = open(args["infofile"],'w')
else: 
	f_info = None

method = "MV"
if args["method"]:
    method = args["method"]

with open(tree_file,'r') as f:
	for line in f:
		tree = Tree.get(data=line,schema=schema,preserve_underscores=True)	
		if method == "MP":
			a_tree = MPR_Tree(ddpTree=tree)
		elif method == "MV":
			a_tree = MV00_Tree(ddpTree=tree)
		elif method == "MV2":
			a_tree = MVDF_Tree(ddpTree=tree)
		else:
			a_tree = None
		if a_tree is None:
			print("Invalid method!")
			break
		else:
			d2currRoot,br2currRoot = a_tree.Reroot()

		if f_info:
			#f_info.write("Head: " + str(head_id) + "\nTail: " + str(tail_id) + "\nEdge_length: " + str(edge_length) + "\nx: " + str(x)+ "\n")
			f_info.write("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
			
		a_tree.tree_as_newick(outfile=outfile,append=True)

if f_info:
	f_info.close()
