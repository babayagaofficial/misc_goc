from ete3 import Tree
import pypangraph as pp
import glob
import os
import pandas as pd
import numpy as np
from itertools import chain


def fitch_bottom_up(node, score):
    """Bottom-up pass to determine possible character states for each node."""
    if node.is_leaf():
        node.add_feature("possible_states", node.states)  # Assign leaf state
    else:
        left_states = fitch_bottom_up(node.children[0], score)
        right_states = fitch_bottom_up(node.children[1], score)

        # Ensure the possible_states attribute exists
        node.add_feature("possible_states", dict())

        # Intersection if possible, otherwise union
        node.possible_states = left_states & right_states if left_states & right_states else left_states | right_states
        score[0] = score[0] if left_states & right_states else score[0] + 1

    return node.possible_states

def fitch_top_down(node, parent=None):
    
    if node.is_leaf():
        node.add_feature("final_states", node.possible_states)
    elif parent is None:
        node.add_feature("final_states", node.possible_states)
    else:
        node.add_feature("final_states", dict())
        # Step I: Check if preliminary nodal set ⊆ parent final set
        if parent.final_states.issubset(node.possible_states):
            # Step II: Eliminate markers not in parent’s final set
            node.final_states = node.possible_states.intersection(parent.final_states)
        else:
            # Step III: Determine whether it was a union or intersection
            left_states = node.children[0].possible_states
            right_states = node.children[1].possible_states
            was_union = not left_states.intersection(right_states)
            
            if was_union:
                # Step IV: Union case
                final_set = node.possible_states.union(parent.final_states)
            else:
                # Step V: Intersection case
                final_set = node.possible_states.copy()
                for state in parent.final_states:
                    if state in left_states or state in right_states:
                        final_set.add(state)
            node.final_states = final_set
                
    
    # Recursively process child nodes
    for child in node.children:
        fitch_top_down(child, node)


def fitch_reconstruction(tree):
    """Executes both phases of Fitch's algorithm on an ete3 Tree."""
    root = tree.get_tree_root()
    score = [0]
    fitch_bottom_up(root, score)
    fitch_top_down(root)
    return score

def print_tree(node, level=0):
    """Utility function to print the tree with assigned states."""
    print("  " * level + f"{node.name}: {node.final_states}, {node.possible_states}")
    for child in node.children:
        print_tree(child, level + 1)

def read_tree(filepath): #read in Newick tree with unlabelled internal nodes, and label the internal nodes
    t = Tree(filepath, format=1)
    i=1
    for n in t.traverse():
        if not n.is_leaf():
            n.name = f"I_{i}"
            i+=1
    return t

def read_block_counts(pangraph_path, tree, min_len):
    graph = pp.Pangraph.from_json(pangraph_path)
    stats_df = graph.to_blockstats_df()
    bl_count = graph.to_blockcount_df()
    bl_count = bl_count.loc[bl_count.index.isin(stats_df[(stats_df["len"]>min_len) & (stats_df["core"]==False)].index)]
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in bl_count.columns if chr.name in plasmid}
    all_trees = {block:tree.copy() for block in bl_count.index}
    for block in all_trees.keys(): 
        for leaf in all_trees[block].iter_leaves():
            #states = {1 if bl_count.loc[block, chr_to_plasmid[leaf.name]]>0 else 0}
            states = {bl_count.loc[block, chr_to_plasmid[leaf.name]]}
            leaf.add_feature("states", states)
    return list(stats_df[stats_df["core"]==True].index), all_trees

def read_core_polymorphisms(pangraph_path, tree):
    graph = pp.Pangraph.from_json(pangraph_path)
    core_aln = graph.core_genome_alignment()
    plasmids = [record.id for record in core_aln]
    # turn the alignment in a numpy matrix
    A = np.array(core_aln)
    # exclude sites with gaps
    non_gap = np.all(A != "-", axis=0)
    A = A[:, non_gap]
    # whether a site is polymorphic
    is_polymorphic_w_recomb = np.any(A != A[0, :], axis=0)
    indices = [i for i in range(0,len(is_polymorphic_w_recomb),500) if is_polymorphic_w_recomb[i:i+500].sum()<=3]
    iterator = []
    for i in range(0,len(is_polymorphic_w_recomb),500):
        if i in indices:
            iterator.append(is_polymorphic_w_recomb[i:i+500])
        else:
            iterator.append([False for i in range(500)])
    is_polymorphic = list(chain.from_iterable(iterator))
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in plasmids if chr.name in plasmid}
    all_trees = {site:tree.copy() for site in range(len(is_polymorphic)) if is_polymorphic[site]}
    for site in all_trees.keys():
        for leaf in all_trees[site].iter_leaves():
            row = plasmids.index(chr_to_plasmid[leaf.name])
            states = {A[row][site]}
            leaf.add_feature("states", states)
    return all_trees
            


def tree_to_tsv(t, recon_trees, core, tsv):
    with open(tsv, "w") as f:
        f.write("genome\tfamily\tlower\thigher\n")
        for block in recon_trees.keys():
            tree = recon_trees[block]
            for node in tree.traverse():
                if not node.is_leaf():
                    counts = list((tree&node.name).final_states)
                    lower = min(counts)
                    upper = max(counts)
                    if not (lower==0 and upper==0):
                        f.write(f"{node.name}\t{block}\t{lower}\t{upper}\n")
        for node in t.traverse():
            if not node.is_leaf():
                for block in core:
                    f.write(f"{node.name}\t{block}\t1\t1\n")

def output_scores(scores, tsv):
    with open(tsv, "w") as f:
        f.write("block\tscore\n")
        for block in scores.keys():
            f.write(f"{block}\t{scores[block][0]}\n")

def root_to_tip_scores(recon_trees, filepath):
    root_to_tip = {block:{} for block in recon_trees.keys()}
    for block in recon_trees.keys():
        tree = recon_trees[block]
        for leaf in tree.iter_leaves():
            score = 0
            path = leaf.get_ancestors()  # from leaf to root
            path = list(reversed(path)) + [leaf]
            states = list(path[0].final_states)
            state = states[0]
            for i in range(1,len(path)):
                if not state in path[i].final_states:
                    score+=1
                    state = path[i].final_states.pop()
                    path[i].final_states.add(state)
            root_to_tip[block][leaf.name] = score
    df = pd.DataFrame(data=root_to_tip)
    df.to_csv(filepath, sep="\t")
    return root_to_tip



# Example usage: Constructing a tree in Newick format
newick = "((A:1,B:1):1,(D:1,C:1):1)Root;"  # A simple binary tree
tree = Tree(newick, format=1)

# Assign states to leaf nodes manually (assuming names represent states)
#for leaf in tree.iter_leaves():
#    leaf.add_feature("state", )  # Treat leaf names as states

'''
node = tree&"A"
node.add_feature("state", {1,2,3})
node = tree&"B"
node.add_feature("state", {1,2})
node = tree&"C"
node.add_feature("state", {1,2,3,4})
node = tree&"D"
node.add_feature("state", {1,2,4})

# Run Fitch's method
fitch_reconstruction(tree)

# Print reconstructed tree states
print_tree(tree)
'''

scores = {}
cluster = snakemake.params.cluster
filepath = snakemake.input.tree
tree = read_tree(filepath) 

pangraph = snakemake.input.pangraph

if snakemake.params.mode == "pangraph":
    core, all_trees = read_block_counts(pangraph, tree, 200)
else:
    all_trees = read_core_polymorphisms(pangraph, tree)

# Run Fitch's method
for block in all_trees.keys():
    k = fitch_reconstruction(all_trees[block])
    scores[block] = k
    # Print reconstructed tree states
    #print_tree(all_trees[block])


if snakemake.params.mode == "pangraph":
    root_to_tip = root_to_tip_scores(all_trees, f"fitch/scores/{cluster}_rtt_scores.tsv")
    tree_to_tsv(tree, all_trees,core,f"fitch/ranges/{cluster}_range.tsv")

    new_root = Tree()
    new_root.name = "I_1"
    new_root.add_child(tree)
    tree.write(format=8, outfile=f"fitch/relabelled_trees/{cluster}.nw")
else:
    root_to_tip = root_to_tip_scores(all_trees, f"fitch/scores/filtered/{cluster}_rtt_snp_scores.tsv") #filtered as in accounting for recombination

