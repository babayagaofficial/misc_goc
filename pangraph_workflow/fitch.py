from ete3 import Tree
import pypangraph as pp
import glob
import os


def fitch_bottom_up(node):
    """Bottom-up pass to determine possible character states for each node."""
    if node.is_leaf():
        node.add_feature("possible_states", node.state)  # Assign leaf state
    else:
        left_states = fitch_bottom_up(node.children[0])
        right_states = fitch_bottom_up(node.children[1])

        # Ensure the possible_states attribute exists
        node.add_feature("possible_states", set())

        # Intersection if possible, otherwise union
        node.possible_states = left_states & right_states if left_states & right_states else left_states | right_states

    return node.possible_states

def fitch_top_down(node, parent=None):
    """Top-down pass using the six-step Fitch algorithm to assign final states."""        
    
    # Step 1: For the root, the final set is simply the preliminary set
    if parent is None or node.is_leaf():
        node.add_feature("final_states", node.possible_states.copy())
    else:
        # Step I: Check if preliminary nodal set ⊆ parent final set
        if node.possible_states.issubset(parent.final_states):
            # Step II: Eliminate markers not in parent’s final set
            node.add_feature("final_states", node.possible_states.intersection(parent.final_states))
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

            node.add_feature("final_states", final_set)
    
    # Recursively process child nodes
    for child in node.children:
        fitch_top_down(child, node)

def fitch_reconstruction(tree):
    """Executes both phases of Fitch's algorithm on an ete3 Tree."""
    root = tree.get_tree_root()
    fitch_bottom_up(root)
    fitch_top_down(root)

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

def read_block_counts(pangraph_path, tree):
    graph = pp.Pangraph.from_json(pangraph_path)
    bl_count = graph.to_blockcount_df()
    stats_df = graph.to_blockstats_df()
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in bl_count.columns if chr.name in plasmid}
    for leaf in tree.iter_leaves():
        present = list(bl_count[chr_to_plasmid[leaf.name]][bl_count[chr_to_plasmid[leaf.name]]>0].index)
        present = set([block for block in present if stats_df.loc[block,"len"]>200])
        leaf.add_feature("state", present)
    return bl_count

def tree_to_tsv(tree, block_counts, tsv):
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in bl_count.columns if chr.name in plasmid}
    with open(tsv, "w") as f:
        f.write("genome\tfamily\tlower\thigher\n")
        for node in tree.traverse():
            state_string = ','.join([str(state) for state in list(node.final_states)])
            f.write(f"{node.name}\t{state_string}")
            subtree = tree&node.name
            leaves = [leaf.name for leaf in subtree.get_leaves()]
            lower = len(node.final_states)
            #upper = max([block_counts[block_counts.index.isin(list(node.final_states))][chr_to_plasmid[leaf]].sum() for leaf in leaves])
            upper = sum([max(block_counts[[chr_to_plasmid[leaf] for leaf in leaves]].loc[state].values) for state in list(node.final_states)])
            print(upper)
            f.write(f"\t{lower}\t{upper}\n")


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
'''
cluster_path = "/home/daria/Documents/projects/ABC/clades/lists"
#clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]
clusters = ["st131_cl416_community_0_subcommunity_503"]
#clusters.remove("st69_cl461_community_0_subcommunity_499")
#clusters.remove("st95_cl461_community_0_subcommunity_502")

for cluster in clusters:
    print(cluster)
    filepath = f"/home/daria/Documents/projects/ABC/clades/trees/{cluster}.nw"
    tree = read_tree(filepath)

    pangraph = f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json"
    bl_count = read_block_counts(pangraph, tree)

    # Run Fitch's method
    fitch_reconstruction(tree)

    # Print reconstructed tree states
    #print_tree(tree)

    for n in tree.traverse():
        print(n.name, len(n.final_states))

    tree_to_tsv(tree,bl_count,f"fitch/ranges/{cluster}_range.tsv")
    tree.write(format=1, outfile=f"fitch/relabelled_trees/{cluster}.nw")
