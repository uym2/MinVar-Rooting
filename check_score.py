from RTT import *

def RTT_score(tree,time,use_as=False,use_qp=False):
    smplTimes = {}
    for line in time:
        sp, t = line.strip().split()
        smplTimes[sp] = float(t)
    if use_qp:
        a_tree = RTT_Tree(smplTimes, ddpTree=tree, solver="QP")
    else:
        a_tree = RTT_Tree(smplTimes, ddpTree=tree, solver="AS")
    a_tree.Bottomup_update()
    a_tree.prepare_root()
    a_tree.Topdown_update()
    return a_tree.opt_score()/a_tree.total_leaves

def check_score(s,use_as=False,use_qp=False):
    i = 1
    while (i <= 20):
        r_ = "Tests/R_RTT/s" + str(s) + "/Trees/s" + str(s) + "_r" + str(i) + ".tre"  # R
        as_ = "Tests/R_RTT/s" + str(s) + "/MyCodeTest/s" + str(s) + "_r" + str(i) + "test_as.tre"  # AS
        qp_ = "Tests/R_RTT/s" + str(s) + "/MyCodeTest/s" + str(s) + "_r" + str(i) + "test_qp.tre"  # QP
        time1 = open("Tests/R_RTT/s" + str(s) + "/Time/s" + str(s) + "_lst" + str(i) + ".txt", "r")
        time2 = open("Tests/R_RTT/s" + str(s) + "/Time/s" + str(s) + "_lst" + str(i) + ".txt", "r")
        if use_as and use_qp:
            tree1 = qp_
            s1 = RTT_score(read_tree_newick(tree1), time1,use_qp=True)
            tree2 = as_
            s2 = RTT_score(read_tree_newick(tree2), time2,use_as=True)
        else:
            tree1 = r_
            s1 = RTT_score(read_tree_newick(tree1), time1)
            if use_as:
                tree2 = as_
            else:
                tree2 = qp_
            s2 = RTT_score(read_tree_newick(tree2), time2)
        print("RTT Score",s2)
        '''
        if abs(s2-s1) > 0.00001:
            print(abs(s2-s1))
        else:
            print("same")
        '''
        i += 1

#check_score(4000,use_as=True)