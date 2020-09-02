from fastroot.RTT import *
from sys import argv


EPSILON_mu = 1e-5

logger = logging.getLogger("compute_RTT")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False

def RTT_score(tree,time):
    smplTimes = {}
    for line in time:
        sp, t = line.strip().split()
        smplTimes[sp] = float(t)
    tree.root.droot, tree.root.SSD, tree.root.SD, tree.root.SDT, tree.SST, tree.ST, n = 0,0,0,0,0,0,0
    tmin = min(smplTimes.values())
    for v in tree.traverse_preorder():
        if not v.is_root():
            v.droot = v.parent.droot + v.edge_length
            if v.is_leaf():
                n += 1
                tree.root.SD += v.droot
                tree.root.SSD += (v.droot ** 2)
                tree.SST += (smplTimes[v.label] ** 2)
                tree.root.SDT += (v.droot * smplTimes[v.label])
                tree.ST += smplTimes[v.label]
    b, e, h, m, r, f = tree.SST, (-2 * tree.root.SDT), n, (-2*tree.ST), 2*tree.root.SD, tree.root.SSD
    y_star = ((m * e)/(2*b) - r) / (2*h - (m*m)/(2*b))
    mu_star = - (e + m * y_star)/(2*b)
    if mu_star < 0:
        mu_star = EPSILON_mu
        y_star = -r/(2*h)
    logger.info("t0:" + str(y_star/mu_star) + " mu:" + str(mu_star))
    RTT = b*mu_star*mu_star + e*mu_star + f + h*y_star*y_star + m*mu_star*y_star + r*y_star
    return RTT/n


myTreeFile = argv[1]
timeFile = argv[2]
time = open(timeFile,"r")
myTree = read_tree_newick(myTreeFile)
logger.info(RTT_score(myTree,time))
