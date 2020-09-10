# from dendropy import Tree
import logging

from treeswift import *
import sys
import math
from fastroot import new_logger

'''
logger = logging.getLogger("Tree_extend.py")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False
'''

class Tree_extend(object):
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):#,logger_id=1,logger_stream=sys.stderr):
        #self.logger = new_logger(__name__+ "_" + str(logger_id),myStream=logger_stream)
        if tree_file:
            self.ddpTree = read_tree(tree_file, schema)
        else:
            self.ddpTree = ddpTree

    def Bottomup_label(self):
        # assign each node a label so that we can later relate to it
        i = 0
        for node in self.ddpTree.traverse_postorder():
            if node.is_leaf():
                node.name = 'L' + str(i)
            else:
                node.name = 'I' + str(i)
            i += 1

    def Topdown_label(self, label_type="all"):
        # assign each node a label so that we can later relate to it
        i = 0

        for node in self.ddpTree.traverse_preorder():
            if node.is_leaf():
                if label_type == "all" or label_type == "leaves":
                    node.name = 'L' + str(i)
                else:
                    node.name = node.label
            else:
                if label_type == "all" or label_type == "internal":
                    node.name = 'I' + str(i)
                else:
                    node.name = node.label
            i += 1

    def Bottomup_update(self):
        for node in self.ddpTree.traverse_postorder():
            self.Node_init(node)
            self.bUp_update(node)

    def Topdown_update(self):
        for node in self.ddpTree.traverse_preorder():
            self.tDown_update(node, self.Opt_function)

    def compute_distances(self):
        D = {}

        def __compute_dRoot__(node, cumm_l):
            if node.is_leaf():
                D[node.name] = cumm_l
            else:
                for child in node.child_nodes():
                    __compute_dRoot__(child, cumm_l + child.edge_length)

        __compute_dRoot__(self.ddpTree.root, 0)
        return D

    def compute_ingroup_distances(self):
        D = []

        def __compute_dLeaf__(node, cumm_l):
            if node.is_leaf():
                D.append(cumm_l)
            else:
                for child in node.child_nodes():
                    __compute_dLeaf__(child, cumm_l + child.edge_length)

        children = self.ddpTree.root.child_nodes()
        crowded_child = None
        maxleaf = -1

        for node in children:
            if node.nleaf > maxleaf:
                maxleaf = node.nleaf
                crowded_child = node

        __compute_dLeaf__(children[1], 0)

        return D

    def filter_branch(self, threshold=None):
        # filter out abnormally long branches
        i = 1
        self.logger.info("Iteration: " + str(i))
        self.Reroot()
        while 1:
            check = self.filter_by_threshold(threshold=threshold)
            if (not check):
                self.logger.info("I could not remove anything more! I stop here!")
                break
            i += 1
            self.logger.info("Iteration: " + str(i))
            self.reset()
            self.Reroot()

    def filter_by_threshold(self, threshold=None, k=3.5):
        if threshold is None:
            threshold = self.compute_threshold(k=k)

        def __filter__(node, cumm_l):
            removed = False
            node.child_removed = False
            for child in node.child_nodes():
                check = __filter__(child, cumm_l + child.edge_length)
                removed = removed or check

            p = node.parent_node
            # if ( cumm_l > threshold ) or ( node.child_removed and len(node.child_nodes()) == 0 ):
            if (cumm_l > threshold) or (node.child_removed and node.num_children() == 0):
                # remove node
                p.remove_child(node)
                # update parent node
                p.child_removed = True
                removed = True
                try:
                    self.logger.info(node.label + " removed")
                except:
                    self.logger.info(node.name + " removed")
            # elif len(node.child_nodes()) == 1:
            elif node.num_child_nodes() == 1:
                # remove node and attach its only child to its parent
                e1 = node.edge_length
                child = node.child_nodes()[0]
                e2 = child.edge_length
                p.remove_child(node)
                node.remove_child(child)
                p.add_child(child)
                child.edge_length = e1 + e2
            return removed

        return __filter__(self.get_root(), 0)

    def compute_threhold(self, k=3.5):
        self.logger.warning("Abstract class! Should never be called")
        return 0

    def reset(self):
        self.logger.warning("Abstract class! Should never be called")

    def find_root(self):
        self.Topdown_label()  # temporarily included for debugging
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

    def opt_score(self):
        self.logger.warning("Abstract class! Should never be called")

    def report_score(self):
        self.logger.warning("Abstract class! Should never be called")

    def Reroot(self):
        self.find_root()
        #self.report_score()
        # d2currRoot = 0
        # br2currRoot = 0
        if self.opt_root != self.ddpTree.root:
            # d2currRoot,br2currRoot = self.reroot_at_edge(self.opt_root.edge, self.opt_root.edge_length-self.opt_x, self.opt_x)
            self.reroot_at_edge(self.opt_root, self.opt_x)
            #self.ddpTree.reroot(self.opt_root,self.opt_x)
        
        # return head_id, tail_id, edge_length, self.opt_x
        # return d2currRoot,br2currRoot

    def Opt_function(self, node):
        self.logger.warning("Abstract method! Should never be called")

    def tree_as_newick(self, outstream=sys.stdout, label_by_name=False):
        # dendropy's method to write newick seems to have problem ...
        self.__write_newick(self.ddpTree.root, outstream, label_by_name=label_by_name)
        outstream.write(";\n")

    #            outstream.write(bytes(";\n", "ascii"))

    def __write_newick(self, node, outstream, label_by_name=False):
        if node.is_leaf():
            if label_by_name:
                outstream.write(str(node.name))
            #                    outstream.write(bytes(str(node.name), "ascii"))
            else:
                try:
                    outstream.write(node.label)
                #                        outstream.write(bytes(node.label, "ascii"))
                except:
                    outstream.write(node.label)
        #                        outstream.write(bytes(str(node.label), "ascii"))
        else:
            outstream.write('(')
            # outstream.write(bytes('(', "ascii"))
            is_first_child = True
            for child in node.child_nodes():
                if is_first_child:
                    is_first_child = False
                else:
                    outstream.write(',')
                #                        outstream.write(bytes(',', "ascii"))
                self.__write_newick(child, outstream, label_by_name=label_by_name)
            outstream.write(')')
        #                outstream.write(bytes(')', "ascii"))
        if not node.is_leaf():
            if label_by_name:
                outstream.write(str(node.name))
            #                    outstream.write(bytes(str(node.name), "ascii"))
            elif node.label is not None:
                outstream.write(str(node.label))
        #                    outstream.write(bytes(str(node.label), "ascii"))

        if not node.edge_length is None:
            outstream.write(":" + str(node.edge_length))

    #                outstream.write(bytes(":" + str(node.edge_length), "ascii"))

    def reroot_at_edge(self, node, length):
    # the method provided by dendropy DOESN'T seem to work ...
    # change edge to opt_root
        length1 = node.edge_length-length
        length2 = length
        if not node:
            return
        head = node #opt_root = v = node
        tail = node.parent #u parent of opt_root
        if not tail:
            return

        if (length2 == 0) and head.is_leaf():
            return 0, 0

        #new_root = self.ddpTree.node_factory()
        new_root = Node()

        tail.remove_child(head)

        new_root.add_child(head)
        head.edge_length=length2

        p = tail.parent
        l = tail.edge_length

        new_root.add_child(tail)
        tail.edge_length = length1

        br2currRoot = 0
        d2currRoot = length1

#            if tail.label == self.ddpTree.root.label:
        if (tail is self.ddpTree.root):
            head = new_root


        while tail is not self.ddpTree.root:
# MAD@ add
            #q = tail.parent #tail should have 2 parents right now: new_root and its old parent
            q = head.parent
# End MAD@ add
            head = tail
            tail = p
            p = tail.parent

            br2currRoot += 1
            d2currRoot += l

            l1 = tail.edge_length
            tail.remove_child(head)
# MAD@ add
            head.parent = q
# End MAD@ add

            head.add_child(tail)
            tail.edge_length=l
            l = l1

        # out of while loop: tail IS now tree.root
        if tail.num_children() == 1:
            # merge the 2 branches of the old root and adjust the branch length
            #sis = [child for child in tail.child_nodes()][0]
            sis = tail.child_nodes()[0]
            l = sis.edge_length
            tail.remove_child(sis)    
            head.add_child(sis)
            sis.edge_length = l + tail.edge_length
            head.remove_child(tail)
            #tail.remove_child(head)

        new_root.name = self.ddpTree.root.name
        self.ddpTree.root.name = "OLD"
        self.ddpTree.root = new_root

### MAD@ add
#            for node in self.ddpTree.traverse_postorder():
#                for child in node.child_nodes():
#                    if child.parent_node is not node:
#                        logger.info("Error found!")
#                        child.parent_node = node
### MAD@ add

        return d2currRoot,br2currRoot

    def get_root(self):
        return self.ddpTree.root


class OGR_Tree(Tree_extend):
    # supportive class to implement outgroup-reroot (OGR = outgroup reroot, hence the name)
    # this rooting method solve the difficulty in finding the root when there are mulitple outgroups
    # and they are not monophyletic. It seeks for the rooting place that maximizes the triplet score
    # of the specified outgroups.
    def __init__(self, outgroups, ddpTree=None, tree_file=None, schema="newick",logger_id=1,logger_stream=sys.stderr):
        super(OGR_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("OGR_Tree_" + str(logger_id),myStream=logger_stream)
        #L = self.ddpTree.leaf_nodes()
        L = []
        for leaf in self.ddpTree.traverse_leaves():
            L.append(leaf)
        self.OGs = set([x.label for x in L if x.label in set(outgroups)])
        self.nOGs = len(self.OGs)
        self.nIGs = len(L) - self.nOGs
        self.max_nTrpls = self.nIGs * self.nOGs * (self.nOGs - 1) / 2 + self.nOGs * self.nIGs * (self.nIGs - 1) / 2
        self.reset()

    def reset(self):
        self.opt_root = self.ddpTree.root
        self.opt_nTrpls = 0

    def Node_init(self, node, nTrpl_in=0, nTrpl_out=0, nOGs=0, nIGs=0):
        node.nTrpl_in = nTrpl_in
        node.nTrpl_out = nTrpl_out
        node.nOGs = nOGs
        node.nIGs = nIGs

    def Opt_function(self, node):
        curr_nTrpls = node.nTrpl_in + node.nTrpl_out
        if curr_nTrpls > self.opt_nTrpls:
            self.opt_nTrpls = curr_nTrpls
            self.opt_root = node
            self.opt_x = node.edge_length / 2  # NOTE: this method does not consider branch length, the *middle point* of the edge is just arbitrarily chosen

    def bUp_update(self, node):
        if node.is_leaf():
            node.nOGs = 1 if node.label in self.OGs else 0
            node.nIGs = 1 if node.nOGs == 0 else 0
        else:
            C = node.child_nodes()

            node.nOGs = sum([c.nOGs for c in C])
            node.nIGs = sum([c.nIGs for c in C])

            node.nTrpl_in = sum([c.nTrpl_in for c in C])

            for i, c1 in enumerate(C):
                for c2 in C[i + 1:]:
                    IG_trpls = c1.nIGs * c2.nIGs * (self.nOGs - node.nOGs)
                    OG_trpls = c1.nOGs * c2.nOGs * (self.nIGs - node.nIGs)
                    node.nTrpl_in += IG_trpls + OG_trpls

    def tDown_update(self, node, opt_function):
        C = node.child_nodes()

        for child in C:
            C1 = [c for c in C if c is not child]
            child.nTrpl_out = node.nTrpl_out

            for i, c1 in enumerate(C1):
                child.nTrpl_out += c1.nTrpl_in
                child.nTrpl_out += (self.nIGs - node.nIGs) * c1.nIGs * child.nOGs
                child.nTrpl_out += (self.nOGs - node.nOGs) * c1.nOGs * child.nIGs

                for c2 in C1[i + 1:]:
                    IG_trpls = c1.nIGs * c2.nIGs * child.nOGs
                    OG_trpls = c1.nOGs * c2.nOGs * child.nIGs

                    child.nTrpl_out += IG_trpls + OG_trpls

            opt_function(child)

    def prepare_root(self):
        pass

    def opt_score(self):
        return self.opt_nTrpls / float(self.max_nTrpls) if self.max_nTrpls != 0 else None

    def report_score(self):
        myScore = self.opt_score()
        if myScore is None:
            self.logger.warning("OG rooting failed because the tree has no outgroup")
        return "Triplet score: " + str(self.opt_score())


class MPR_Tree(Tree_extend):
    # supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)G
    def __init__(self, ddpTree=None, tree_file=None, schema="newick",logger_id=1,logger_stream=sys.stderr):
        super(MPR_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("MPR_Tree_" + str(logger_id),myStream=logger_stream)
        self.reset()

    def reset(self):
        self.max_distance = -1
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, max_in=None, max_out=-1):
        node.max_in = max_in if max_in else [0, 0]
        node.max_out = max_out

    def Opt_function(self, node):
        m = max(node.max_in)
        curr_max_distance = m + node.max_out
        x = (node.max_out - m) / 2
        if curr_max_distance > self.max_distance and x >= 0 and x <= node.edge_length:
            self.max_distance = curr_max_distance
            self.opt_x = x
            self.opt_root = node

    def bUp_update(self, node):
        if not node.is_leaf():
            node.max_in = []
            for child in node.child_nodes():
                node.max_in.append(max(child.max_in) + child.edge_length)

    def tDown_update(self, node, opt_function):
        child_idx = 0
        for child in node.child_nodes():
            child.max_out = max([node.max_out] + [node.max_in[k] for k in range(len(node.max_in))
                                                  if k != child_idx]) + child.edge_length
            opt_function(child)
            child_idx += 1

    def prepare_root(self):
        pass

    def compute_threhold(self, k=3.5):
        self.logger.warning("Trying to compute threshold for MPR_Tree, which is not supported.")
        return 0

    def opt_score(self):
        return self.max_distance / 2

    def report_score(self):
        return "Tree height: " + str(self.opt_score())




