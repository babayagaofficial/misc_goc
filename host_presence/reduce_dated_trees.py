from ete3 import Tree, TreeStyle, NodeStyle, random_color, CircleFace, TextFace, faces

filepath = "/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_dated_73.nwk"
t = Tree(filepath, format=1)

root=t.get_tree_root()
child_1, child_2 = root.get_children()
sub_1 = child_1.copy()
sub_2 = child_2.copy()
if len(sub_1.get_descendants())>len(child_2.get_descendants()):
    t = sub_1
else:
    t = sub_2

root=t.get_tree_root()
child_1, child_2 = root.get_children()
sub_1 = child_1.copy()
sub_2 = child_2.copy()
if len(sub_1.get_descendants())>len(child_2.get_descendants()):
    t = sub_1
else:
    t = sub_2

print(len(t))

ts = TreeStyle()
ts.show_leaf_name = False

t.render(f'reduced.png', tree_style=ts)

#t.write(format=1, outfile=f"reduced_73.nwk")