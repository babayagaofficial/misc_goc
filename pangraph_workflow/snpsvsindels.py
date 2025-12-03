import pypangraph as pp
import numpy as np
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree

cluster_path = "/home/daria/Documents/projects/ABC/clades/lists_phylofactor"
clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]
clusters = sorted(clusters)
#clusters = ["st131_cl416_community_0_subcommunity_503"]
#clusters.remove("st69_cl461_community_0_subcommunity_499")
#clusters.remove("st95_cl461_community_0_subcommunity_502")
snps = []
indels = []

colours = {
    "c_0_sc_503": "#CC6677",
    "c_0_sc_496": "#332288",
    "c_0_sc_502": "#AA4499",
    "c_0_sc_500": "#117733",
    "c_0_sc_493": "#88CCEE",
    "c_0_sc_494": "#882255",
    "c_0_sc_501": "#44AA99",
    "c_0_sc_499": "#999933",
    "c_0_sc_498": "#DDCC77",
    "c_0_sc_485": "#f3b2ff",
    "c_0_sc_490": "#006575",
    "c_0_sc_487": "#6d75d7",
    "c_0_sc_489": "#613d00"
}

total = pd.DataFrame(columns=["indels", "snps"])
avg_indel_rate = {"clade":[], "rate":[]}
scores = []

for cluster in clusters:
    print(cluster)

    pangraph = f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json"
    graph = pp.Pangraph.from_json(pangraph)
    core_aln = graph.core_genome_alignment()
    A = np.array(core_aln)
    core_size = len(A[0])/1000

    split = cluster.split("_")
    name = '_'.join(['c',split[3],'sc',split[5]])
    clade = '_'.join([split[0],split[1]])

    filepath = f"/home/daria/Documents/projects/ABC/clades/trees/{cluster}.nw"
    tree = Tree(filepath, format=1)

    indel_per_block = pd.read_csv(f"fitch/scores/{cluster}_rtt_scores.tsv", sep="\t", index_col=0)
    snp_per_site = pd.read_csv(f"fitch/scores/filtered/{cluster}_rtt_snp_scores.tsv", sep="\t", index_col=0)
    parsimony_score = pd.read_csv(f"fitch/scores/{cluster}_scores.tsv", sep="\t", index_col=0)
    print(parsimony_score)
    try:
        scores.append(parsimony_score["score"].sum(numeric_only=True))
    except:
        scores.append(0)

    time = {}
    root = tree.get_tree_root()
    for leaf in tree.iter_leaves():
        time[leaf.name] = tree.get_distance(root,leaf)

    time = pd.Series(time, name="time")

    indel_rtt = indel_per_block.sum(axis=1, numeric_only=True)
    snp_rtt = snp_per_site.sum(axis=1, numeric_only=True)

    rtt = pd.concat([indel_rtt,snp_rtt], axis=1)
    rtt = pd.concat([rtt,time], axis=1)
    rtt.columns = ["indels", "snps", "time"]
    rtt.fillna(0, inplace=True)

    rtt_rel_indels = rtt["indels"]/rtt["time"]
    rtt_rel_indels.name = "indels"
    rtt_rel_snps = rtt["snps"]/rtt["time"]
    rtt_rel_snps.name = "snps"
    rtt_rel = pd.DataFrame([rtt_rel_indels,rtt_rel_snps]).transpose()
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([name for el in rtt_rel.index], columns=["subcommunity"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([core_size for el in rtt_rel.index], columns=["core_size"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([clade for el in rtt_rel.index], columns=["clade"],index=rtt_rel.index)],axis=1)



    fig, ax = plt.subplots(figsize=(8,8))
    sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
    sns.scatterplot(data=rtt_rel, y="snps", x="indels", ax=ax)
    plt.ylabel("snps", fontsize=12)
    plt.xlabel("indels", fontsize=12)
    ax.set_title("Parsimonous root-to-tip number of snps vs indels", fontsize=14)
    plt.savefig(f"snpsvsindels/{cluster}.png")
    print(rtt_rel['snps'].corr(rtt_rel['indels']))

    avg_indel_rate["rate"].append(rtt_rel_indels.mean())
    avg_indel_rate["clade"].append(cluster)

    total = pd.concat([total,rtt_rel])
'''
print(total)

fig, ax = plt.subplots(figsize=(8,8))
ax.axline(xy1=(0, 0), slope=1, alpha=0.5)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
sns.scatterplot(data=total[total["subcommunity"]=="c_0_sc_503"], y="snps", x="indels", ax=ax, hue="clade", alpha=0.8)
plt.ylabel("snps", fontsize=12)
plt.xlabel("indels", fontsize=12)
ax.set_title("Parsimonous root-to-tip number of snps vs indels", fontsize=14)
plt.savefig(f"snpsvsindels/c0_sc503.png")
print(total['snps'].corr(total['indels']))

fig, ax = plt.subplots(figsize=(8,8))
ax.axline(xy1=(0, 0), slope=1, alpha=0.5)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
sns.scatterplot(data=total, y="snps", x="indels", ax=ax, hue="subcommunity",size="core_size", alpha=0.8,palette=colours)
plt.ylabel("snps", fontsize=12)
plt.xlabel("indels", fontsize=12)
ax.set_title("Parsimonous root-to-tip number of snps vs indels", fontsize=14)
plt.savefig(f"snpsvsindels/total.png")




print(scores)
fig, ax = plt.subplots(figsize=(8,8))
ax.axline(xy1=(0, 0), slope=1, alpha=0.5)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
sns.displot(data=scores, kde=True)
plt.savefig(f"snpsvsindels/dis.png")
'''
rates = pd.DataFrame(avg_indel_rate)
rates.to_csv("avg_indel_rtt_rates.tsv", sep="\t")