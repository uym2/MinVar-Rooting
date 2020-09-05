from fastroot.Tree_extend import *

'''
logger = logging.getLogger("MinVar")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False
'''

class minVAR_Base_Tree(Tree_extend):
    # supportive base class to implement VAR-reroot, hence the name
    def __init__(self, ddpTree=None, tree_file=None, schema="newick",logger_id=1,logger_stream=sys.stderr):
        super(minVAR_Base_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("MinVar_Tree_" + str(logger_id),myStream=logger_stream)
        self.reset()

    def reset(self):
        self.minVAR = None
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, sum_in=0, sum_total=0, var=-1):
        node.sum_in = sum_in
        node.sum_total = sum_total
        node.nleaf = nleaf
        node.var = var

    def Opt_function(self, node, a, b, c):
        self.logger.info("Abstract method! Should never be called")

    def compute_dRoot_VAR(self):
        cumm = {'ssq': 0, 'sum': 0}

        def compute_dRoot(node, cumm_l):
            if node.is_leaf():
                cumm['ssq'] += cumm_l ** 2
                cumm['sum'] += cumm_l
            else:
                for child in node.child_nodes():
                    compute_dRoot(child, cumm_l + child.edge_length)

        compute_dRoot(self.get_root(), 0)
        N = self.get_root().nleaf
        root_var = cumm['ssq'] / N - (cumm['sum'] / N) ** 2
        self.get_root().var = root_var

    def bUp_update(self, node):
        if node.is_leaf():
            node.nleaf = 1
            node.sum_in = 0
        else:
            node.nleaf = 0
            node.sum_in = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.sum_in += child.sum_in + child.nleaf * child.edge_length

    def Update_var(self, child, node, edge_length):
        alpha = 2 * (node.sum_total - 2 * (child.sum_in + child.nleaf * edge_length)) / self.total_leaves
        beta = 1 - 2 * float(child.nleaf) / self.total_leaves
        a = 1 - beta * beta
        b = alpha - 2 * node.sum_total * beta / self.total_leaves
        c = node.var
        child.var = a * edge_length * edge_length + b * edge_length + c
        return a, b, c

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.sum_total = node.sum_total + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            a, b, c = self.Update_var(child, node, child.edge_length)
            opt_function(child, a, b, c)

    def prepare_root(self):
        root = self.get_root()
        root.sum_total = root.sum_in
        self.compute_dRoot_VAR()
        self.total_leaves = root.nleaf

    def opt_score(self):
        return self.minVAR

    def report_score(self):
        return "MinVar score: " + str(self.opt_score())

class MVDF_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + deepest node + factorization
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MVDF_Tree, self).__init__(ddpTree, tree_file, schema)
        self.deep_node = None

    def reset(self):
        super(MVDF_Tree, self).reset()
        self.deep_node = None

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            #                curr_minVAR = a*x*x + b*x + c
            factor = float(node.nleaf) / self.total_leaves
            factor = factor * (1 - factor)
            curr_minVAR = (a * x * x + b * x + c) / factor

            if node.var < node.parent_node.var:
                deep_node = node
            else:
                deep_node = node.parent_node

            updateNeed = False
            if (self.deep_node is None) or (deep_node.var < self.deep_node.var):
                self.deep_node = deep_node
                self.minVAR = curr_minVAR
                updateNeed = True
            elif (self.deep_node is deep_node) and (curr_minVAR < self.minVAR):
                self.minVAR = curr_minVAR
                updateNeed = True

            if updateNeed:
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        factor = float(self.opt_root.nleaf) / self.total_leaves
        factor = factor * (1 - factor)
        rootVar = self.minVAR * factor
        self.logger.info(mean)
        self.logger.info(rootVar)
        std = math.sqrt(rootVar)
        return mean + k * std


class MVD0_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + deepest node + no factorization
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MVD0_Tree, self).__init__(ddpTree, tree_file, schema)
        self.deep_node = None

    def reset(self):
        super(MVD0_Tree, self).reset()
        self.deep_node = None

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            curr_minVAR = a * x * x + b * x + c

            if node.var < node.parent_node.var:
                deep_node = node
            else:
                deep_node = node.parent_node

            updateNeed = False
            if (self.deep_node is None) or (deep_node.var < self.deep_node.var):
                self.deep_node = deep_node
                self.minVAR = curr_minVAR
                updateNeed = True
            elif (self.deep_node is deep_node) and (curr_minVAR < self.minVAR):
                self.minVAR = curr_minVAR
                updateNeed = True

            if updateNeed:
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        self.logger.info(mean)
        self.logger.info(self.minVAR)
        std = math.sqrt(self.minVAR)
        return mean + k * std


class MV0F_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + no deepest node + factorization
    #        def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
    #            super().__init__(ddpTree, tree_file, schema)

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            #                curr_minVAR = a*x*x + b*x + c
            factor = float(node.nleaf) / self.total_leaves
            factor = factor * (1 - factor)
            curr_minVAR = (a * x * x + b * x + c) / factor
            if self.minVAR is None or curr_minVAR < self.minVAR:
                self.minVAR = curr_minVAR
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        factor = float(self.opt_root.nleaf) / self.total_leaves
        factor = factor * (1 - factor)
        rootVar = self.minVAR * factor
        self.logger.info(mean)
        self.logger.info(rootVar)
        std = math.sqrt(rootVar)
        return mean + k * std


class MV00_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + no deepest node + no factorization
    #        def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
    #            super().__init__(ddpTree, tree_file, schema)

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            curr_minVAR = a * x * x + b * x + c
            if self.minVAR is None or curr_minVAR < self.minVAR:
                self.minVAR = curr_minVAR
                self.opt_root = node
                self.opt_x = node.edge_length - x

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        self.logger.info(mean)
        self.logger.info(self.minVAR)
        std = math.sqrt(self.minVAR)
        return mean + k * std


class MBR_Tree(Tree_extend):
    # supportive class to implement midpoint balance root
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MBR_Tree, self).__init__(ddpTree, tree_file, schema)

        self.BPs = []  # BPs : balance points
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, sum_in=0, sum_out=-1):
        self.nleaf = nleaf
        self.sum_in = sum_in
        self.sum_out = sum_out

    def Opt_function(self, node):
        nleaf = node.nleaf
        mean_in = node.sum_in / nleaf
        mean_out = node.sum_out / (self.total_leaves - nleaf)
        x = (mean_out - mean_in) / 2
        if x >= 0 and x <= node.edge_length:
            self.BPs.append((node, x, mean_in + x))
            node.x = x
            node.mean = mean_in + x
        else:
            node.x = None
            node.mean = None

    def bUp_update(self, node):
        node.sum_in = 0
        if node.is_leaf():
            node.nleaf = 1
        else:
            node.nleaf = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.sum_in += child.sum_in + child.nleaf * child.edge_length

    def tDown_update(self, node, opt_function):
        child_idx = 0
        for child in node.child_nodes():
            child.sum_out = (node.sum_out + node.sum_in + child.edge_length *
                             (self.total_leaves - 2 * child.nleaf) - child.sum_in)
            opt_function(child)
            child_idx += 1

    def prepare_root(self):
        root = self.get_root()
        root.sum_out = 0
        self.total_leaves = root.nleaf
        root.x = None
        root.mean = None

    def list_balance_points(self):
        self.Topdown_label()
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

        for (node, x, mean) in self.BPs:
            if node.is_leaf():
                #                  self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))
                self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))
            else:
                self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))

    def build_balance_tree(self):
        self.Topdown_label()  # keep this step for now for debugging purpose
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

        # self.list_balance_points()

        self.balance_tree = self.ddpTree.extract_tree()

        # bottom up pruning
        for node in self.balance_tree.traverse_postorder():
            node.type = "real"
            node.BPbelow = False

            '''if node.is_leaf():
                self.logger.info("parent: " + node.label)# + "\t" + str(node.extraction_source.x))
            else:
                self.logger.info("parent: " + node.label)#+ "\t" + str(node.extraction_source.x))'''

            for ch in node.child_nodes():
                '''try:
                    self.logger.info("child: " + ch.label)# + "\t" + str(ch.extraction_source.x))
                except:
                    self.logger.info("child: " + ch.label) #+ "\t" + str(ch.extraction_source.x))'''

                if ch.BPbelow or (ch.extraction_source.x is not None):
                    node.BPbelow = True
                # node.BPbelow = node.BPbelow or ch.BPbelow or (ch.extraction_source.x is not None)

                if not ch.BPbelow:
                    # remove the whole clade under ch
                    # for ch1 in ch.child_nodes():
                    #    ch.remove_child(ch1)
                    edgelen = ch.edge_length
                    node.remove_child(ch)

                    if ch.extraction_source.x is not None:
                        # add a new node p at the balance point
                        # set p to be a child of node (edge length ch.edge_length - x)
                        # add a new node ch1 to be another child of p (edge length ch.mean)
                        edgelen = ch.edge_length

                        # p = self.ddpTree.node_factory()
                        # ch1 = self.ddpTree.node_factory()
                        p = Node()
                        ch1 = Node()

                        p.type = "bp"  # bp: balance-point
                        p.ref_child = ch.extraction_source  # link p to the original tree (for later use after finding midpoint)
                        ch1.type = "dm"  # dm: dummy

                        # node.remove_child(ch)
                        node.add_child(p)
                        p.add_child(ch1)

                        p.edge_length = edgelen - ch.extraction_source.x
                        ch1.edge_length = ch.extraction_source.mean

                elif ch.extraction_source.x is not None:
                    # add a new node p at the balance point
                    # set p to be a child of node (edge length ch.edge_length - x)
                    # set ch to be a child of p (edge length x)
                    # add a new node ch1 to be another child of p (edge length ch.mean)

                    edgelen = ch.edge_length

                    # p = self.ddpTree.node_factory()
                    p = Node()
                    # ch1 = self.ddpTree.node_factory()
                    ch1 = Node()

                    p.type = "bp"
                    p.ref_child = ch.extraction_source  # link p to the original tree (for later use after finding midpoint)
                    ch1.type = "dm"

                    node.remove_child(ch)
                    node.add_child(p)
                    p.add_child(ch)
                    p.add_child(ch1)

                    ch.edge_length = ch.extraction_source.x
                    p.edge_length = edgelen - ch.extraction_source.x
                    ch1.edge_length = ch.extraction_source.mean

                    # topdown pruning
        node = self.balance_tree.root
        nchild = len(node.child_nodes())
        while nchild > 0 and nchild < 2:
            # node has less than 2 children
            temp = node
            node = node.child_nodes()[0]
            temp.remove_child(node)
            if node.type == "dm":
                node = temp
                break
            nchild = len(node.child_nodes())

        self.balance_tree.root = node
        self.balance_tree.root.edge_length = None
        # balance_tree.root = None

        # mptre = MPR_Tree(ddpTree=balance_tree)
        # mptre.tree_as_newick()

        # return balance_tree

    def find_root(self):
        self.build_balance_tree()
        mptre = MPR_Tree(ddpTree=self.balance_tree)
        mptre.tree_as_newick()
        mptre.find_root()

        self.logger.info(mptre.opt_root.type)

        if mptre.opt_root.type == "bp":
            self.opt_root = mptre.opt_root.ref_child
            self.opt_x = mptre.opt_root.ref_child.x + mptre.opt_x
        elif mptre.opt_root.type == "dm":
            self.logger.info("Hmm... Is it possible that a dummy was found as the opt_root?")
        else:
            self.opt_root = mptre.opt_root.extraction_source
            self.opt_x = mptre.opt_x

        self.logger.info(self.opt_root.label)
        self.logger.info(self.opt_x)
