import pypangraph as pp
import os
import glob

graph = pp.Pangraph.from_json(snakemake.input.pangraph)
block_stats = graph.to_blockstats_df()
path_dict = graph.to_path_dictionary()

with open(snakemake.output.unimog, "w") as f:
    for plasmid in path_dict.keys():
        f.write(f">{plasmid}\n")
        for block in path_dict[plasmid]:
            strand = '' if block[1] else '-'
            length = block_stats.loc[block[0],"len"]
            if length>200:
                f.write(f"{strand}{block[0]} ")
        f.write(")\n")
