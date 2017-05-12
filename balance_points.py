#! /usr/bin/env python

from Tree_extend import MBR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MBR_Tree(tree_file=tree_file)
a_tree.list_balance_points()
