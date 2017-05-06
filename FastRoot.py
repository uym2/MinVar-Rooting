#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import MPR_Tree,MVR_Tree
from dendropy import Tree,TreeList

from os.path import splitext
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i','--input',required=True,help="input file")
parser.add_argument('-m','--method',required=False,help="method: MP for midpoint, MV for minVAR. Default is MV")
parser.add_argument('-o','--outfile',required=False,help="specify output file")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")
parser.add_argument('-f','--infofile',required=False,help="write info of the new root to file; mostly for research and debugging purposes")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)

schema=args["schema"] if args["schema"] else "newick"

outfile = args["outfile"] if args["outfile"] else None

f_info = open(args["infofile"],'w') if args["infofile"] else None

method = args["method"] if args["method"] else "MV"

try:
	os.remove(outfile)
except:
	pass

try:
	os.remove(args["infofile"])
except:
	pass

	

with open(tree_file,'r') as f:
	for line in f:
		tree = Tree.get(data=line,schema=schema)	
		if method == "MP":
			a_tree = MPR_Tree(ddpTree=tree)
		elif method == "MV":
			a_tree = MVR_Tree(ddpTree=tree)
		else:
			print("Invalid method!")
			break
		d2currRoot,br2currRoot = a_tree.Reroot()

		if f_info:
			f_info.write("d2currRoot: " + str(d2currRoot) + "\nbr2currRoot: " + str(br2currRoot) + "\n")
			
		a_tree.tree_as_newick(outfile=outfile,append=True)

if f_info:
	f_info.close()
