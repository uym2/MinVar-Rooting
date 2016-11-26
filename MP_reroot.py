from dendropy import Tree,Node
from sys import argv,stdout

tree_file = argv[1]

class record(object):
	def __init__(self,max_left=0,max_right=0,max_outside=0):
		self.max_left = max_left
		self.max_right = max_right
		self.max_outside = max_outside

def tree_as_newick(tree,outfile=None):
	if outfile:
		outstream = open(outfile,'w')
	else:
		outstream = stdout
	
	_write_newick(tree.seed_node,outstream)
	outstream.write(";")

	if outfile:
		outstream.close()	

def _write_newick(node,outstream):
	if node.is_leaf():
			outstream.write(str(node.taxon))
	else:
		outstream.write('(')
		is_first_child = True
		for child in node.child_node_iter():
			if is_first_child:
				is_first_child = False
			else:
				outstream.write(',')
			_write_newick(child,outstream)
		outstream.write(')')
	if not node.is_leaf() and node.label:
		outstream.write(str(node.label))
	if node.edge_length:
		outstream.write(":"+str(node.edge_length))

def deroot_tree(tree,parent_side):
	root = tree.seed_node
	left_child,right_child = [child for child in node.child_node_iter()]
	combined_edge = left_child.edge_length + right_child.edge_length
	root.remove_child(left_child)
	root.remove_child(right_child)
	if parent_side == 'r':
		right_child.add_child(left_child,edge_length=combined_edge)
	else:
		left_child.add_child(right_child,edge_length=combined_edge)
	
def reroot_at_edge(tree,edge,length1,length2,new_root=None):
	head = edge.head_node
	tail = edge.tail_node
	
	#left_child,right_child = [child for child in opt_root.child_node_iter()]
	#print left_child,right_child

	#print head.child_node_iter()
	#print head.child_node_iter()
	#print head.label
	#print tail.label

	if not new_root:
		#new_root = tree.node_factory()
		new_root = Node()		

	tail.remove_child(head)
		
	new_root.add_child(head)
	head.edge_length=length2
	head.edge.length = length2

	p = tail.parent_node
	l = tail.edge_length

	new_root.add_child(tail)
	tail.edge_length=length1
	tail.edge.length=length1
	#print p.label
	is_root_edge = False
	if tail.label == tree.seed_node.label:
		head = new_root
        	is_root_edge = True
		l = tail.edge_length

	while tail.label != tree.seed_node.label:
		#print l	
		head = tail
		tail = p
		p = tail.parent_node
		
		#print tail.label
		if tail.label == tree.seed_node.label:
			break
		l1 = tail.edge_length

		tail.remove_child(head)

		p = tail.parent_node
		l1 = tail.edge_length

		head.add_child(tail)
		tail.edge.length=l
		tail.edge_length=l
		
		l = l1
		
	# out of while loop: tail IS now tree.seed_node
	#print tail.label
	#print head.label
	sis = [child for child in tail.child_node_iter() if child.label != head.label][0]
	#print sis.label
	l1 = sis.edge_length
	print l,l1
	if not is_root_edge:
		tail.remove_child(head)
	tail.remove_child(sis)	
	head.add_child(sis)
	sis.edge.length=l+l1
	sis.edge_length=l+l1	

	tree.seed_node = new_root
	#tree.update_bipartitions()
	#tree.write(path="rooted.tre",schema="newick")	
	
	#tree.print_plot()
	#print tree.as_string("newick")

a_tree = Tree.get_from_path(tree_file,"newick")

i = 0
record_arr = []


for node in a_tree.postorder_node_iter():
	#print i
	node.label = i
	if node.is_leaf():
		node_record = record()
		record_arr.append(node_record)
	else:
		left_child,right_child = [child for child in node.child_node_iter()]
		max_left = max(record_arr[left_child.label].max_left,record_arr[left_child.label].max_right)+left_child.edge_length
		max_right = max(record_arr[right_child.label].max_left,record_arr[right_child.label].max_right)+right_child.edge_length
		node_record = record(max_left=max_left,max_right=max_right)
		record_arr.append(node_record)
	i = i+1
#print a_tree.seed_node.label
record_arr[a_tree.seed_node.label].max_outside = 0
max_distance = 0
opt_root = a_tree.seed_node
#left_child,right_child = [child for child in opt_root.child_node_iter()]
#print left_child,right_child
opt_x = 0

for node in a_tree.preorder_node_iter():
	child_idx = 1
	node_record = record_arr[node.label]
	node_maxl = node_record.max_left
	node_maxr = node_record.max_right
	node_maxo = node_record.max_outside

	for child in node.child_node_iter():
		if child_idx == 1:
			record_arr[child.label].max_outside = max(node_maxr,node_maxo)+child.edge_length
		else:
			record_arr[child.label].max_outside = max(node_maxl,node_maxo)+child.edge_length
		m = max(record_arr[child.label].max_left,record_arr[child.label].max_right) 
		curr_max_distance = m + record_arr[child.label].max_outside
		x = (record_arr[child.label].max_outside - m)/2
		if curr_max_distance > max_distance and x >= 0 and x <= child.edge_length:
			max_distance = curr_max_distance
			opt_x = x
			opt_root = child
		child_idx = child_idx+1

#for node in a_tree.postorder_node_iter():
	#print str(node.label) + " " + str(node.edge_length) + " " + str(record_arr[node.label].max_left) + " " + str(record_arr[node.label].max_right) + " " + str(record_arr[node.label].max_outside)
	#print max_left[edge.label]
	#print max_right[edge.label]
#print max_distance
#print opt_x
opt_idx = opt_root.label
#print str(opt_root.label) + " " + str(opt_root.edge_length) + " " + str(record_arr[opt_idx].max_left) + " "+ str(record_arr[opt_idx].max_right) + " " + str(record_arr[opt_idx].max_outside)
#print a_tree.seed_node.label
reroot_at_edge(a_tree,opt_root.edge,opt_root.edge_length-opt_x,opt_x)

'''
for node in a_tree.levelorder_node_iter():
	#if not node.is_leaf():
	#	print node._as_newick_string()
	#print node.label, node.edge.length, node.taxon
	#print node.taxon
	print node.label
	if not node.is_leaf():
		left,right = [child for child in node.child_node_iter()]
		print left.label, left.edge_length
		print  right.label, right.edge_length
	else:
		print node.taxon
	#	print [child.label for child in node.child_node_iter()]
'''		
#p = opt_root.parent_node
#a_tree.reroot_at_edge(opt_root.edge)

#a_tree.update_bipartitions()

#a_tree.write(path=tree_file.split(".tre")[0]+"_rooted.tre",schema="newick")

tree_as_newick(a_tree,outfile=tree_file.split(".tre")[0]+"_rooted.tre")
