import pypangraph as pp
import numpy as np
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree

def root_to_tip_scores(tree, dists, filepath):
    root_to_tip = {}
    for leaf in tree.iter_leaves():
        path = leaf.get_ancestors()  # from leaf to root
        path = list(reversed(path)) + [leaf]
        score = sum([dists[node.name] for node in path if not node.is_root()])
        root_to_tip[leaf.name] = score
    with open(filepath, "w") as f:
            for key in root_to_tip.keys():
                f.write(f"{key}\t{root_to_tip[key]}\n")
    return root_to_tip

def edgewise_distances(filepath):
    dists = {}
    with open(filepath, "r") as f:
        for line in f:
            child, parent, dist = line.strip().split(" ")
            dists[child] = float(dist)
    return dists

def plot_fcn_compare(rtt_rel, name, colours, type = "clade"):
    fig, ax = plt.subplots(figsize=(10,8))
    #ax.axline(xy1=(0, 0), slope=1, alpha=0.8)

    for _, row in rtt_rel.iterrows():
        plt.plot([row['Fitch'], row['DCJ-Indel']], [row['SNPs'], row['SNPs']], color='gray', linewidth=1, zorder=0)

    df_long = rtt_rel.melt(id_vars=['SNPs','subcommunity','core_size','clade'], value_vars=['DCJ-Indel','Fitch'], var_name='type', value_name='rates')
    sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
    if type == "total":
        sns.scatterplot(data=df_long, y="SNPs", x="rates", style="type", ax=ax, hue="subcommunity", alpha=0.8, palette=colours)
    elif type == "subcomm":
        sns.scatterplot(data=df_long, y="SNPs", x="rates", style="type", ax=ax, hue="clade", alpha=0.8)
    else:
        sns.scatterplot(data=df_long, y="SNPs", x="rates", style="type", ax=ax)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    plt.axis('square')
    plt.ylabel("snps", fontsize=12)
    plt.xlabel("structural", fontsize=12)
    ax.set_title("Parsimonous root-to-tip rates", fontsize=14)
    plt.savefig(f"compare_rates/{name}.png")