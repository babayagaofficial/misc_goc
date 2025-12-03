from ete3 import Tree, TreeStyle, NodeStyle, random_color, CircleFace, TextFace, faces

filepath = "/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_dated_95.nwk"
filepath2 = "/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_95.nwk"
t = Tree(filepath, format=1)
t2 = Tree(filepath2, format=1)

root=t.get_tree_root()
child_1, child_2 = root.get_children()
sub_1 = child_1.copy()
sub_2 = child_2.copy()
if len(sub_1.get_descendants())>len(sub_2.get_descendants()):
    outgroup = sub_2.get_leaf_names()
else:
    outgroup = sub_1.get_leaf_names()

root=t2.get_tree_root()
child_1, child_2 = root.get_children()
sub_1 = child_1.copy()
sub_2 = child_2.copy()
if len(sub_1.get_descendants())>len(sub_2.get_descendants()):
    outgroup2 = sub_2.get_leaf_names()
else:
    outgroup2 = sub_1.get_leaf_names()

print(outgroup)
print(outgroup2)
print(t.compare(t2))