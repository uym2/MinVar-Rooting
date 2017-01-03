#! /usr/bin/env python

# usage: python MP_reroot.py <tree_file>

import os
from Tree_extend import Tree_extend
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
parser.add_argument('-o','--outfile',required=False,help="specify output file")
parser.add_argument('-s','--schema',required=False,help="schema of your input treefile. Default is newick")

args = vars(parser.parse_args())

tree_file = args["input"]
base_name,ext = splitext(tree_file)
schema=args["schema"] if args["schema"] else "newick"

if args["outfile"]:
	outfile = args["outfile"]
else:
	outfile = base_name + "_labeled" + ext

try:
	os.remove(outfile)
except:
	pass

try:
	os.remove(args["infofile"])
except:
	pass


trees = TreeList.get_from_path(tree_file,schema)

for tree in trees:
	a_tree = Tree_extend(ddpTree=tree)
	a_tree.Topdown_label()
	a_tree.tree_as_newick(outfile=outfile,append=True)
