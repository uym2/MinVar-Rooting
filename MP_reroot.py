from dendropy import Tree,Node
from sys import argv

tree_file = argv[1]

class record(object):
	def __init__(self,max_left=0,max_right=0,max_outside=0):
		self.max_left = max_left
		self.max_right = max_right
		self.max_outside = max_outside

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
	
def reroot_at_edge(tree,edge,length1,length2,parent_side):
	head = edge.head_node
	tail = edge.tail_node
	new_root = Node()
	new_root.add_child(head,length2)
	new_root.add_child(tail,length1)
	tail.remove_child(head)
	if tail != tree.seed_node:
		deroot_tree(tree,parent_side)
	tree.seed_node = new_root		

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
print a_tree.seed_node.label
record_arr[a_tree.seed_node.label].max_outside = 0
max_distance = 0
opt_root = a_tree.seed_node
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

for node in a_tree.postorder_node_iter():
	print str(node.label) + " " + str(node.edge_length) + " " + str(record_arr[node.label].max_left) + " " + str(record_arr[node.label].max_right) + " " + str(record_arr[node.label].max_outside)
	#print max_left[edge.label]
	#print max_right[edge.label]
print max_distance
print opt_x
opt_idx = opt_root.label
print str(opt_root.label) + " " + str(opt_root.edge_length) + " " + str(record_arr[opt_idx].max_left) + " "+ str(record_arr[opt_idx].max_right) + " " + str(record_arr[opt_idx].max_outside)

a_tree.reroot_at_edge(opt_root.edge,length1=opt_root.edge_length-2,length2=2,update_bipartitions=False)

a_tree.write(path=tree_file.split(".tre")[0]+"_rooted.tre",schema="newick")
