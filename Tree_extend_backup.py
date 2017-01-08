from dendropy import Tree,Node
import copy

class Tree_extend(object):
		def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
			if tree_file:
				self.ddpTree = Tree.get_from_path(tree_file,schema)
			else:
				self.ddpTree = ddpTree
			self.Tree_records = Tree_records

		def New_record(self):
			print("Abstract method! Should never be called!")

		def Bottomup_label(self):
			# assign each node a label so that we can later relate to it
			i = 0	
			for node in self.ddpTree.postorder_node_iter():
				if not node.is_leaf():
					node.label = 'I'+str(i)
				else:
					node.label = 'L'+str(i)
				i = i+1

		
		def Topdown_label(self):
			# assign each node a label so that we can later relate to it
			i = 0	
			for node in self.ddpTree.preorder_node_iter():
				if not node.is_leaf():
					node.label = 'I'+str(i)
				else:
					node.label = 'L'+str(i)
				i = i+1

		def Bottomup_update(self):
			i = 0
			for node in self.ddpTree.postorder_node_iter():
				node_record = self.New_record()
				node.idx = i
				node_record.Bottomup_update(node,self.Tree_records)
				self.Tree_records.append(node_record)
				i = i+1
			
		def Topdown_update(self):
			for node in self.ddpTree.preorder_node_iter():
				self.Tree_records[node.idx].Topdown_update(node,self.Tree_records,self.Opt_function)

		def Reroot(self):
			self.Bottomup_update()
			self.prepare_root()
			self.Topdown_update()

			if self.opt_root.is_leaf():
				head_id = self.opt_root.taxon.label
			else:
				head_id = self.opt_root.label
			tail_id = self.opt_root.parent_node.label if self.opt_root.parent_node else None
			edge_length = self.opt_root.edge_length		

			if self.opt_root != self.ddpTree.seed_node:
				self.__reroot_at_edge(self.opt_root.edge,self.opt_root.edge_length-self.opt_x,self.opt_x)
			
			return head_id, tail_id, edge_length, self.opt_x
			
		def Opt_function(self,node):
			print("Abstract method! Should never be called")


		def tree_as_newick(self,outfile=None,append=False):
		# dendropy's method to write newick seems to have problem ...
			if outfile:
				outstream = open(outfile,'a') if append else open(outfile,'w')
			else:
				outstream = VARout
		
			self.__write_newick(self.ddpTree.seed_node,outstream)
			outstream.write(";\n")
			if outfile:
				outstream.close()	

		def __write_newick(self,node,outstream):
			if node.is_leaf():
					outstream.write(node.taxon.label)
					#outstream.write(str(node.label))
			else:
				outstream.write('(')
				is_first_child = True
				for child in node.child_node_iter():
					if is_first_child:
						is_first_child = False
					else:
						outstream.write(',')
					self.__write_newick(child,outstream)
				outstream.write(')')
			if not node.is_leaf() and node.label is not None:
					outstream.write(str(node.label))
			
			if not node.edge_length is None:
				outstream.write(":"+str(node.edge_length))

		def __reroot_at_edge(self,edge,length1,length2,new_root=None):
		# the method provided by dendropy DOESN'T seem to work ...
			head = edge.head_node
			tail = edge.tail_node
		
			if not new_root:
				#new_root = Node()		
				new_root = self.ddpTree.node_factory()

			tail.remove_child(head)
			
			new_root.add_child(head)
			head.edge_length=length2

			p = tail.parent_node
			l = tail.edge_length

			new_root.add_child(tail)
			tail.edge_length=length1

			if tail == self.ddpTree.seed_node:
				head = new_root

			while tail != self.ddpTree.seed_node:
				head = tail
				tail = p
				p = tail.parent_node

				l1 = tail.edge_length
				tail.remove_child(head)

				head.add_child(tail)
				tail.edge_length=l
				l = l1
				
			# out of while loop: tail IS now tree.seed_node
			if tail.num_child_nodes() < 2:
				# merge the 2 branches of the old root and adjust the branch length
				sis = [child for child in tail.child_node_iter()][0]
				l = sis.edge_length
				tail.remove_child(sis)	
				head.add_child(sis)
				sis.edge_length=l+tail.edge_length
				head.remove_child(tail)
			
			new_root.label = self.ddpTree.seed_node.label
			self.ddpTree.seed_node = new_root

		def get_root_idx(self):
			return self.ddpTree.seed_node.idx

		def get_root(self):
			return self.ddpTree.seed_node

class MPR_Tree(Tree_extend):
	# supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)
		def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
			if tree_file:
				self.ddpTree = Tree.get_from_path(tree_file,schema)
			else:
				self.ddpTree = copy.deepcopy(ddpTree)
			self.Tree_records = Tree_records
			self.max_distance = -1
			self.opt_root = self.ddpTree.seed_node
			self.opt_x = 0

		def New_record(self):
			return MPR_Node_record()

		def Opt_function(self,node):
			m = max(self.Tree_records[node.idx].max_in) 
			curr_max_distance = m + self.Tree_records[node.idx].max_out
			x = (self.Tree_records[node.idx].max_out - m)/2
			if curr_max_distance > self.max_distance and x >= 0 and x <= node.edge_length:
				self.max_distance = curr_max_distance
				self.opt_x = x
				self.opt_root = node

		def prepare_root(self):
			pass

class minVAR_Tree(Tree_extend):
	# supportive class to implement VAR-reroot, hence the name
		def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
			if tree_file:
				self.ddpTree = Tree.get_from_path(tree_file,schema)
			else:
				self.ddpTree = copy.deepcopy(ddpTree)
			self.Tree_records = Tree_records
			self.minVAR = None
			self.opt_root = self.ddpTree.seed_node
			self.opt_x = 0

		def New_record(self):
			return minVAR_Node_record()

		def Opt_function(self,node,a,b,c):
			x = -b/(2*a)
			if x >= 0 and x <= node.edge_length:
				curr_minVAR = a*x*x + b*x + c
				if self.minVAR is None or curr_minVAR < self.minVAR:
					self.minVAR = curr_minVAR
					self.opt_root = node
					self.opt_x = node.edge_length-x

		def compute_dRoot_VAR(self):
			cumm = {'ssq':0,'sum':0}
			def compute_dRoot(node,cumm_l):
				if node.is_leaf():
					cumm['ssq'] += cumm_l**2
					cumm['sum'] += cumm_l
				else:
					for child in node.child_node_iter():
						compute_dRoot(child,cumm_l+child.edge_length)

			compute_dRoot(self.get_root(),0)
			N = self.Tree_records[self.get_root_idx()].nleaf
			root_var = cumm['ssq']/N-(cumm['sum']/N)**2	
			self.Tree_records[self.get_root_idx()].var = root_var
			self.minVAR = root_var

		def prepare_root(self):
			self.Tree_records[self.get_root_idx()].sum_total = self.Tree_records[self.get_root_idx()].sum_in
			self.compute_dRoot_VAR()

class MDR_Tree(Tree_extend):
# supportive class to implement mean difference root (mdr = mean difference reroot, hence the name)
        def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
                if tree_file:
                        self.ddpTree = Tree.get_from_path(tree_file,schema)
                else:
                        self.ddpTree = copy.deepcopy(ddpTree)
                self.Tree_records = Tree_records
                self.min_MD = None
                self.opt_root = self.ddpTree.seed_node
                self.opt_x = 0

        def New_record(self):
                return MDR_Node_record()

        def Opt_function(self,node):
		nleaf = self.Tree_records[node.idx].nleaf
		mean_in = sum(self.Tree_records[node.idx].sum_in)/nleaf
        	mean_out = self.Tree_records[node.idx].sum_out/(MDR_Node_record.total_leaves-nleaf)         
		x = (mean_out - mean_in)/2 
                if x < 0:
			x = 0
		elif x > node.edge_length:
			x = node.edge_length
                curr_MD = abs(mean_out-mean_in-2*x)
		
		if self.min_MD is None or curr_MD < self.min_MD: 
			self.min_MD = curr_MD
                        self.opt_x = x
                        self.opt_root = node

        def prepare_root(self):
		ridx = self.get_root_idx()
                self.Tree_records[ridx].sum_out = 0
		#child_idx = 0
		#means = []
		#for child in self.get_root().child_node_iter():
		#	means.append(self.Tree_records[ridx].sum_in[child_idx]/self.Tree_records[child.idx].nleaf)
		#self.min_MD = abs(means[0]-means[1]) # temporary solution: assume 


class MPR2_Tree(Tree_extend):
	# supportive class to implement MP2 rooting (extension of midpoint)
		def __init__(self,ddpTree=None,tree_file=None,schema="newick",Tree_records=[]):
			if tree_file:
				self.ddpTree = Tree.get_from_path(tree_file,schema)
			else:
				self.ddpTree = copy.deepcopy(ddpTree)
			self.Tree_records = Tree_records
			self.max_distance = -1
			self.opt_root = self.ddpTree.seed_node
			self.opt_x = 0

		def New_record(self):
			return MPR2_Node_record()

		#def score_reroot(self,node):
			# compute the new score if the tree was rerooted at the node specified
		#	new_score = self.Tree_records[node.idx].cumm_score

		def Opt_function(self,node):
			#m = max(self.Tree_records[node.idx].max_in) 
			#curr_max_distance = m + self.Tree_records[node.idx].max_out
			#x = (self.Tree_records[node.idx].max_out - m)/2
			#if curr_max_distance > self.max_distance and x >= 0 and x <= node.edge_length:
			#	self.max_distance = curr_max_distance
			#	self.opt_x = x
			#	self.opt_root = node
			pass
		def prepare_root(self):
			pass


class Node_record(object):
	def __init__(self):
		pass
	
	def Bottomup_update(self,node,Tree_records):
		print ("Just an abstract method! You should never see this message. Otherwise please check your code!")

class MPR_Node_record(Node_record):
# supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)
	def __init__(self,max_in=[0,0],max_out=-1):
#		self.old_label=old_label
		self.max_in = max_in
		self.max_out = max_out

	def Bottomup_update(self,node,Tree_records):
		if not node.is_leaf():
			self.max_in=[]
			for child in node.child_node_iter():
				self.max_in.append(max(Tree_records[child.idx].max_in) + child.edge_length)	
	
	def Topdown_update(self,node,Tree_records,opt_function):
		child_idx = 0
		for child in node.child_node_iter():	
			Tree_records[child.idx].max_out = max([self.max_out]+[self.max_in[k] for k in range(len(self.max_in)) if k != child_idx])+child.edge_length
			opt_function(child)
			child_idx = child_idx+1


class minVAR_Node_record(Node_record):
# supportive class to implement VAR-reroot, hence the name
	total_leaves = 0
	def __init__(self,nleaf=1,sum_in=0,sum_total=0,var=-1):
		self.sum_in = sum_in
		self.sum_total = sum_total
		self.nleaf = nleaf
		self.var = var

	def Bottomup_update(self,node,Tree_records):
		if node.is_leaf():
			self.nleaf = 1
			self.sum_in = 0
		else:
			self.nleaf = 0
			self.sum_in = 0
			for child in node.child_node_iter():
				self.nleaf += Tree_records[child.idx].nleaf
				self.sum_in += Tree_records[child.idx].sum_in + Tree_records[child.idx].nleaf*child.edge_length
			minVAR_Node_record.total_leaves = max(minVAR_Node_record.total_leaves,self.nleaf)	
	
	def Update_var(self,p_record,edge_length):
		alpha = 2*( p_record.sum_total-2*(self.sum_in+self.nleaf*edge_length) )/minVAR_Node_record.total_leaves
		beta = 1-2*float(self.nleaf)/minVAR_Node_record.total_leaves
		a = 1-beta*beta
		b = alpha-2*p_record.sum_total*beta/minVAR_Node_record.total_leaves
		c = p_record.var
		self.var = a*edge_length*edge_length + b*edge_length + c
		return a,b,c

	def Topdown_update(self,node,Tree_records,opt_function):
		for child in node.child_node_iter():	
			Tree_records[child.idx].sum_total = Tree_records[node.idx].sum_total + (minVAR_Node_record.total_leaves-2*Tree_records[child.idx].nleaf)*child.edge_length
			a,b,c = Tree_records[child.idx].Update_var(self,child.edge_length)
			opt_function(child,a,b,c)


class MDR_Node_record(Node_record):
# supportive class to implement mean-difference reroot (mdr = mean difference reroot, hence the name)
	total_leaves = 0
	def __init__(self,nleaf=1,sum_in=[0,0],sum_out=-1):
		self.nleaf = nleaf
		self.sum_in = sum_in
		self.sum_out = sum_out

	def Bottomup_update(self,node,Tree_records):
		if node.is_leaf():
			self.nleaf = 1
			self.sum_in = [0,0]
		else:
			self.nleaf = 0
			self.sum_in=[]
			for child in node.child_node_iter():
				self.nleaf += Tree_records[child.idx].nleaf
				s = sum(Tree_records[child.idx].sum_in) + Tree_records[child.idx].nleaf*child.edge_length
				self.sum_in.append(s)	
			MDR_Node_record.total_leaves = max(MDR_Node_record.total_leaves,self.nleaf)	
	def Topdown_update(self,node,Tree_records,opt_function):
		child_idx = 0
		for child in node.child_node_iter():	
			Tree_records[child.idx].sum_out = self.sum_out + sum([self.sum_in[k] for k in range(len(self.sum_in)) if k != child_idx]) + (MDR_Node_record.total_leaves - Tree_records[child.idx].nleaf)*child.edge_length
			opt_function(child)
			child_idx = child_idx+1

class MPR2_Node_record(Node_record):
# supportive class to implement MPR2
	def __init__(self,max_in=[[0]],max_out=[0,0]):
		self.max_in = max_in
		self.max_out = max_out
		self.score = 0
		self.cumm_score = 0

	def MDM_score(lists):
		# MDM = Min of Difference of Max
		n = len(lists)
		score = None		
		for i in range(n-1):
			for j in range(i+1,n):
				delta = min([abs(x-y) for x in lists[i] for y in lists[j]])
				if score is None or delta < score:
					score = delta
		return score

	def score_clade():
		self.score = MDM_score(self.max_in)

	def Bottomup_update(self,node,Tree_records):
		if not node.is_leaf():
			self.max_in=[]
			self.score_clade()
			self.cumm_score = self.score

			for child in node.child_node_iter():
				child_max_in = [ max(L)+child.edge_length for L in Tree_records[child.idx].max_in ]
				self.max_in.append(child_max_in)	
				self.cumm_score += Tree_records[child.idx].cumm_score
	

	def Topdown_update(self,node,Tree_records,opt_function):
		child_idx = 0
		for child in node.child_node_iter():	
			if self.max_out:
				Tree_records[child.idx].max_out = [ max(self.max_in[k])+child.edge_length for k in range(len(self.max_in)) if k != child_idx ] + [ max(self.max_out)+child.edge_length ]
			else:
				if len(self.max_in) > 2:
					Tree_records[child.idx].max_out = [ max(self.max_in[k])+child.edge_length for k in range(len(self.max_in)) if k != child_idx ]  
				else:
					k = 0 if child_idx else 1
					Tree_records[child.idx].max_out = [ x + child.edge_length for x in self.max_in[k] ]

			opt_function(child)
			child_idx = child_idx+1
