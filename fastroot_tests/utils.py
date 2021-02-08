from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from treeswift import *
from fastroot.MinVar import MV00_Tree, minVAR_Base_Tree
from fastroot.RTT import RTT_Tree

# check if two newick strings are equal (same clades):
def compute_clades(tree):
    clades = []
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.clade = [node.label]
        else:
            node.clade = []
            for child in node.children:
                node.clade += child.clade
        clades.append(node.clade)
    return clades

def check_clades_include(cl1, cl2):
    # check if cl1 is a subset of cl2
    for c1 in cl1:
        found = False
        for c2 in cl2:
            if sorted(c1) == sorted(c2):
                found = True
                break
        if not found:
            return False
    return True

def check_two_nwk_str(nwk1, nwk2):
    tree1 = read_tree_newick(nwk1)
    tree2 = read_tree_newick(nwk2)

    L1 = sorted([x.label for x in tree1.traverse_leaves()])
    L2 = sorted([x.label for x in tree2.traverse_leaves()])
    if L1 != L2:
        return False

    cld1 = compute_clades(tree1)
    cld2 = compute_clades(tree2)

    return check_clades_include(cld1, cld2)


# Root input trees
def root_trees(inputFile,method='MV',timeFile=None,OGFile=None,t0=0.0):
    # define OGs in file, not list
    # reading outgroups
    if OGFile:
        method = 'OG'
        OGs = []
        f = open(OGFile, 'r')
        for line in f:
                OGs.append(line.strip())
        f.close()

    # reading sampling times
    if timeFile:
        smplTimes = {}
        f = open(timeFile, 'r')
        for line in f:
            sp, t = line.strip().split()
            smplTimes[sp] = float(t) + t0
        method = 'RTT'
        f.close()

    METHOD2FUNC = {'MP': MPR_Tree, 'MV': MV00_Tree, 'OG': OGR_Tree, 'RTT': RTT_Tree}

    # read and root each tree
    score = {}  # dictionary of tree number and its score
    tree_list = []
    branches = {}
    inputTrees = open(inputFile)
    for i, line in enumerate(inputTrees.readlines()):
        tree = read_tree(line, 'newick')
        if method == 'OG':
            a_tree = OGR_Tree(OGs, ddpTree=tree, logger_id=i + 1)
        elif method == 'RTT':
            a_tree = RTT_Tree(smplTimes, ddpTree=tree, logger_id=i + 1)
        else:
            a_tree = METHOD2FUNC[method](ddpTree=tree, logger_id=i + 1)

        outTrees, annTree = a_tree.Reroot()

        if method == 'RTT':
            score[i+1] = a_tree.return_values()  # tree number: RTT score, mu, t0
        else:
            score[i+1]= a_tree.opt_score()

        tree_list.append(outTrees[0])

        branches[i + 1] = sorted([child.edge_length for child in read_tree_newick(outTrees[0]).root.child_nodes()])
    inputTrees.close()
    return score, tree_list, branches


def branch_lengths(file):
    branch_lengths = {}
    f = open(file,'r')
    for i, line in enumerate(f):
        l1,l2 = line.strip().split()
        branch_lengths[i+1] = sorted([float(l1), float(l2)])
    f.close()
    return branch_lengths

def score_from_file(file, RTT=False):
    score = {}
    f = open(file, 'r')
    if RTT:
        for i,line in enumerate(f):
            RTT_score, mu, t0 = line.strip().split()
            score[i+1] = [float(RTT_score), float(mu), float(t0)]
    else:
        for i, line in enumerate(f):
            score[i+1] = float(line.strip())
    f.close()
    return score
