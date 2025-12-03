import pypangraph as pp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from collections import defaultdict

graph = pp.Pangraph.from_json("/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/st131_cl343_community_0_subcommunity_496/pangraph.json")

print(graph.to_blockstats_df())

path_dict = graph.to_path_dictionary()

block_stats = graph.to_blockstats_df()


# dictionary to assign a new random color to each block
block_color = defaultdict(lambda: plt.cm.rainbow(np.random.rand()))

fig, ax = plt.subplots(figsize=(8, 6))

y = 0
for path_name, path in path_dict.items():
    x = 0
    for block_id, block_strand in path:

        L = block_stats.loc[block_id, "len"] # block consensus length
        is_core = block_stats.loc[block_id, "core"]

        # block color
        color = block_color[block_id] if is_core else "lightgray"
        block_color[block_id] = mpl.colors.to_hex(color)

        height = 0.8 if is_core else 0.6 # block thickness

        ax.barh(y, L, left=x, height=height, color=color)

        x += L
    y += 1

ax.set_yticks(range(len(path_dict)))
ax.set_yticklabels(path_dict.keys())
ax.set_xlabel("length (bp)")
plt.show()



# find MSUs
threshold_len = 200  # minimal length of core blocks to consider
MSU_mergers, MSU_paths, MSU_len = pp.minimal_synteny_units(graph, threshold_len)

# dictionary to assign colors to MSUs
cmap = mpl.colormaps["rainbow"]
color_generator = (cmap(i / len(MSU_len)) for i in range(len(MSU_len)))
colors = defaultdict(lambda: next(color_generator))

fig, ax = plt.subplots(figsize=(8, 5))

for i, (iso, path) in enumerate(MSU_paths.items()):
    for j, node in enumerate(path.nodes):
        ax.barh(i, 1, left=j, color=colors[node.id])
        if not node.strand:
            ax.arrow(j + 1, i, -0.8, 0, head_width=0.2, head_length=0.2)
ax.set_yticks(range(len(MSU_paths)))
ax.set_yticklabels(list(MSU_paths.keys()))
plt.show()