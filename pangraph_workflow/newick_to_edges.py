from ete3 import Tree, NodeStyle, TreeStyle, TextFace, faces


filepath = "/home/daria/Documents/projects/ABC/clades/trees/st131_cl416_community_0_subcommunity_503.nw"

t = Tree(filepath, format=1)

i=1
for n in t.traverse():
    if not n.is_leaf():
        n.name = f"I_{i}"
        i+=1

with open("st131_cl416_community_0_subcommunity_503_tree.txt", "w") as f:
    for n in t.traverse():
        children = n.children
        for child in children:
            f.write(f"{child.name}\t{n.name}\n")