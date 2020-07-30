from RTT import *
from sys import argv

def RTT_score(tree,time):
    smplTimes = {}
    for line in time:
        sp, t = line.strip().split()
        smplTimes[sp] = float(t)
    tree.root.droot, tree.root.SSD, tree.root.SDT, tree.SST, n = 0,0,0,0,0
    for v in tree.traverse_preorder():
        if not v.is_root():
            v.droot = v.parent.droot + v.edge_length
            if v.is_leaf():
                n += 1
                tree.root.SSD += (v.droot ** 2)
                tree.SST += (smplTimes[v.label] ** 2)
                tree.root.SDT += (v.droot * smplTimes[v.label])

    b, e, f = tree.SST, (-2 * tree.root.SDT), tree.root.SSD
    mu_star = -e / (2 * b)
    RTT = b * mu_star * mu_star + e * mu_star + f
    return RTT/n


def check_score(s,use_as=False,use_qp=False):
    i = 1
    while (i <= 10):
        r_ = "Tests/R_RTT/s" + str(s) + "/Trees/s" + str(s) + "_r" + str(i) + ".tre"  # R
        as_ = "Tests/R_RTT/s" + str(s) + "/MyCodeTest/s" + str(s) + "_r" + str(i) + "test_as.tre"  # AS
        qp_ = "Tests/R_RTT/s" + str(s) + "/MyCodeTest/s" + str(s) + "_r" + str(i) + "test_qp.tre"  # QP
        time1 = open("Tests/R_RTT/s" + str(s) + "/Time/s" + str(s) + "_lst" + str(i) + ".txt", "r")
        time2 = open("Tests/R_RTT/s" + str(s) + "/Time/s" + str(s) + "_lst" + str(i) + ".txt", "r")
        if use_as and use_qp:
            tree1 = qp_
            s1 = RTT_score(read_tree_newick(tree1), time1)
            tree2 = as_
            s2 = RTT_score(read_tree_newick(tree2), time2)
        else:
            tree1 = r_
            s1 = RTT_score(read_tree_newick(tree1), time1)
            if use_as:
                tree2 = as_
            else:
                tree2 = qp_
            s2 = RTT_score(read_tree_newick(tree2), time2)
        #print("RTT Score",s2)
        if abs(s2-s1) > 0.0001:
            print(abs(s2-s1))
        else:
            print("same")
        i += 1

check_score(1000,use_as=True)

'''
myTreeFile = argv[1]
timeFile = argv[2]
time = open(timeFile,"r")
myTree = read_tree_newick(myTreeFile)
print(RTT_score(myTree,time))
'''