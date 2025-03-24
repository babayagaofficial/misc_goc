import pypangraph as pp
import os
import glob

cluster_path = "/home/daria/Documents/projects/ABC/clades/lists"
clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]
clusters.remove("st69_cl461_community_0_subcommunity_499")
clusters.remove("st95_cl461_community_0_subcommunity_502")

for cluster in clusters:
    graph = pp.Pangraph.from_json(f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json")
    block_stats = graph.to_blockstats_df()
    path_dict = graph.to_path_dictionary()

    with open(f"unimogs/{cluster}.unimog", "w") as f:
        for plasmid in path_dict.keys():
            f.write(f">{plasmid}\n")
            for block in path_dict[plasmid]:
                strand = '' if block[1] else '-'
                length = block_stats.loc[block[0],"len"]
                if length>200:
                    f.write(f"{strand}{block[0]} ")
            f.write(")\n")