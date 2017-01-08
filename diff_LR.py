#! /usr/bin/env python


from Tree_extend import MDR_Tree
from sys import argv
from os.path import splitext

tree_file = argv[1]
base_name,ext = splitext(tree_file)

a_tree = MDR_Tree(tree_file=tree_file)

print a_tree.diff_of_means()




