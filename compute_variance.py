from treeswift import *
from sys import argv
from RTT import *
from Tree_extend import *

def compute_variance(tree):  # tree is rooted
    tree.root.droot = 0
    D = []
    for v in tree.traverse_preorder():
        if not v.is_root():
            # v.set_edge_length(4)  -- must have defined edge lengths
            v.droot = v.get_parent().droot + v.edge_length
            if v.is_leaf():
                D.append(v.droot)
    SSD, ST = 0, 0
    n = len(D)
    for d in D:
        SSD += (d ** 2)
        ST += d
    var = (SSD / n) - ((ST / n) ** 2)
    return var
    # output: non-negative float # which is the variance of the root-to-tip distance of the input

'''
myTreeFile = argv[1]
myTree = read_tree_newick(myTreeFile)
#print(compute_variance(myTree))
x = Tree_extend(myTree)
y = minVAR_Base_Tree(x)
#print(y.report_score())

z = RTT_Tree(x)
z.find_root()
'''