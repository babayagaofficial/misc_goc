from ete3 import Tree
import pypangraph as pp
import glob
import os
import pandas as pd


def fitch_bottom_up(node):
    """Bottom-up pass to determine possible character states for each node."""
    if node.is_leaf():
        node.add_feature("possible_states", node.states)  # Assign leaf state
    else:
        left_states = fitch_bottom_up(node.children[0])
        right_states = fitch_bottom_up(node.children[1])

        # Ensure the possible_states attribute exists
        node.add_feature("possible_states", dict())

        # Intersection if possible, otherwise union
        node.possible_states = left_states & right_states if left_states & right_states else left_states | right_states

    return node.possible_states

def single_fitch_top_down(node, parent=None):
    """Top-down pass using the six-step Fitch algorithm to assign final states."""        
    
    if node.is_leaf():
        node.add_feature("final_states", node.possible_states)
    elif parent is None:
        pass
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
        single_fitch_top_down(child, node)

def fitch_top_down(tree):
    root = tree.get_tree_root()
    recon_trees = []
    while root.possible_states: #Want to reconstruct all possible trees, so go through all elements in candidate set of root
        recon_tree = tree.copy()
        # Step 1: For the root, choose arbitrary count for ambiguous blocks
        recon_root = recon_tree.get_tree_root()
        recon_root.add_feature("final_states", dict())
        recon_root.final_states = {root.possible_states.pop()}
        single_fitch_top_down(recon_root)
        recon_trees.append(recon_tree)
    return recon_trees



def fitch_reconstruction(tree):
    """Executes both phases of Fitch's algorithm on an ete3 Tree."""
    root = tree.get_tree_root()
    fitch_bottom_up(root)
    recon_trees = fitch_top_down(root)
    return recon_trees

def print_tree(node, level=0):
    """Utility function to print the tree with assigned states."""
    print("  " * level + f"{node.name}: {node.final_states}")
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
            states = {bl_count.loc[block, chr_to_plasmid[leaf.name]]}
            leaf.add_feature("states", states)
    return list(stats_df[stats_df["core"]==True].index), all_trees

def tree_to_df(tree): #after reconstruction, get counts as a DataFrame for both ancestral genomes and leaves
    root = tree.get_tree_root()
    counts = {block:[] for block in root.final_states.keys()}
    index = []
    for node in tree.traverse():
        index.append(node.name)
        for block in node.final_states.keys():
            counts[block].append(node.final_states[block].pop())
    df = pd.DataFrame(counts, index=index)
    return df


def tree_to_tsv(recon_trees, core, tsv): #unfinished
    with open(tsv, "w") as f:
        f.write("genome\tfamily\tlower\thigher\n")
        for node in tree.traverse():
            for block in recon_trees.keys():
                f.write(f"{node.name}\t{block}\t")
                counts = [(tree&node.name).final_states.pop() for tree in recon_trees[block]]
                lower = min(counts)
                upper = max(counts)
                f.write(f"{lower}\t{upper}\n")
            for block in core:
                f.write(f"{node.name}\t{block}\t1\t1\n")





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

cluster_path = "/home/daria/Documents/projects/ABC/clades/lists"
clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]
#clusters = ["st131_cl416_community_0_subcommunity_503"]
clusters.remove("st69_cl461_community_0_subcommunity_499")
clusters.remove("st95_cl461_community_0_subcommunity_502")

for cluster in clusters:
    print(cluster)
    filepath = f"/home/daria/Documents/projects/ABC/clades/trees/{cluster}.nw"
    tree = read_tree(filepath)

    pangraph = f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json"
    core, all_trees = read_block_counts(pangraph, tree, 200)

    # Run Fitch's method
    recon_trees = {}
    for block in all_trees.keys():
        recon_trees[block] = fitch_reconstruction(all_trees[block])
        # Print reconstructed tree states
        #for tree in recon_trees[block]:
        #    print(block)
        #    print_tree(tree)

    tree_to_tsv(recon_trees,core,f"fitch/ranges/{cluster}_range.tsv")
    #tree.write(format=1, outfile=f"fitch/relabelled_trees/{cluster}.nw")
