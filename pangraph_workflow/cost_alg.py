from ete3 import Tree
#import pypangraph as pp
import numpy as np

def read_block_intervals(pangraph_path, tree, min_len):
    graph = pp.Pangraph.from_json(pangraph_path)
    stats_df = graph.to_blockstats_df()
    bl_count = graph.to_blockcount_df()
    bl_count = bl_count.loc[bl_count.index.isin(stats_df[(stats_df["len"]>min_len) & (stats_df["core"]==False)].index)]
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in bl_count.columns if chr.name in plasmid}
    all_trees = {block:tree.copy() for block in bl_count.index}
    for block in all_trees.keys(): 
        for leaf in all_trees[block].iter_leaves():
            #states = {1 if bl_count.loc[block, chr_to_plasmid[leaf.name]]>0 else 0}
            count = bl_count.loc[block, chr_to_plasmid[leaf.name]]
            states = (count, count)
            leaf.add_feature("states", states)
    return list(stats_df[stats_df["core"]==True].index), all_trees


def recon_bottom_up(node, cost):
    if node.is_leaf():
        node.add_feature("possible_states", node.states)  # Assign leaf state
    else:
        left_states = recon_bottom_up(node.children[0], cost)
        right_states = recon_bottom_up(node.children[1], cost)

        # Ensure the possible_states attribute exists
        node.add_feature("possible_states", dict())

        if left_states[0]<=right_states[0]:
            node.possible_states = (min(left_states[1], right_states[0]), max(left_states[1], right_states[0]))
            cost[0] += max(0,right_states[0]-left_states[1])
        else:
            node.possible_states = (min(left_states[0], right_states[1]), max(left_states[0], right_states[1]))
            cost[0] += max(0,left_states[0]-right_states[1])

    return node.possible_states

def recon_top_down(node, parent=None):

    if node.is_leaf():
        node.add_feature("final_states", node.possible_states)
    elif parent is None:
        node.add_feature("final_states", node.possible_states)
    else:
        node.add_feature("final_states", dict())
        # Step I: Check if preliminary nodal set ⊆ parent final set, i.e. check if parent was formed through intersection
        if parent.final_states[0]>=node.possible_states[0] and parent.final_states[1]<=node.possible_states[1]:
            # Step II: Eliminate markers not in parent’s final set
            node.final_states = parent.final_states
        else:
            # Step III: Determine whether it was a union or intersection
            left_states = node.children[0].possible_states
            right_states = node.children[1].possible_states
            was_union = (left_states[0]<=right_states[0] and left_states[1]<right_states[0]) or (right_states[0]<left_states[0] and right_states[1]<left_states[0])
            
            if was_union:
                # Step IV: Union case: add ALL parental states
                if parent.final_states[1]<=node.possible_states[0]: #parent interval is to the left of node interval
                    final_set = (parent.final_states[0], node.possible_states[1])
                else: 
                    final_set = (node.possible_states[0], parent.final_states[1])
            else:
                # Step V: Intersection case: add parental states shared with at least one child
                if parent.final_states[1]<=node.possible_states[0]: #to the left
                    if left_states[0]<=right_states[0]: 
                        final_set = (max(left_states[0], parent.final_states[0]), node.possible_states[1])
                    else:
                        final_set = (max(right_states[0], parent.final_states[0]), node.possible_states[1])
                else:
                    if left_states[0]<=right_states[0]:
                        final_set = (node.possible_states[0], min(right_states[1], parent.final_states[1]))
                    else:
                        final_set = (node.possible_states[0], min(left_states[1], parent.final_states[1]))

            node.final_states = final_set

    # Recursively process child nodes
    for child in node.children:
        recon_top_down(child, node)

def reconstruction(tree):
    """Executes both phases on an ete3 Tree."""
    root = tree.get_tree_root()
    score = [0]
    recon_bottom_up(root, score)
    #recon_top_down(root)
    return score

def intervals_to_tsv(t, recon_trees, core, tsv):
    with open(tsv, "w") as f:
        f.write("genome\tfamily\tlower\thigher\n")
        for block in recon_trees.keys():
            tree = recon_trees[block]
            for node in tree.traverse():
                if not node.is_leaf():
                    counts = (tree&node.name).final_states
                    lower = counts[0]
                    upper = counts[1]
                    if not (lower==0 and upper==0):
                        f.write(f"{node.name}\t{block}\t{lower}\t{upper}\n")
        for node in t.traverse():
            if not node.is_leaf():
                for block in core:
                    f.write(f"{node.name}\t{block}\t1\t1\n")

def print_tree(node, level=0):
    """Utility function to print the tree with assigned states."""
    print("  " * level + f"{node.name}: {node.possible_states}")
    for child in node.children:
        print_tree(child, level + 1)



newick = "(((A:1,B:1):1,(D:1,C:1):1):1,(E:1,F:1):1)Root;" 
tree = Tree(newick, format=1)

node = tree&"A"
node.add_feature("states", (0,0))
node = tree&"B"
node.add_feature("states", (2,2))
node = tree&"C"
node.add_feature("states", (4,4))
node = tree&"D"
node.add_feature("states", (7,7))
node = tree&"E"
node.add_feature("states", (6,6))
node = tree&"F"
node.add_feature("states", (8,8))


score = reconstruction(tree)

print(score)
print_tree(tree)
