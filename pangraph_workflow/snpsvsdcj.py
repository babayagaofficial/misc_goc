from ete3 import Tree
import pandas as pd
import os
import glob
import numpy as np
import pypangraph as pp
import matplotlib.pyplot as plt
import seaborn as sns

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

def chromosomal_snps(st,clade,subcomm,snp_tree):
    clade_num = clade.replace("cl","")
    path = f"/home/daria/Documents/projects/ABC/host_presence/phylofactor/{st}/{subcomm}/clades/{clade_num}.txt"
    members= []
    to_hash = {}
    with open(path, "r") as f:
        for line in f:
            member = line.strip()
            members.append(member.replace('#','_'))
            to_hash[member.replace('#','_')] = member
    mrca = snp_tree.get_common_ancestor(members)
    chr_snps = {}
    for member in members:
        chr_snps[to_hash[member]]=snp_tree.get_distance(member,mrca)
    chr_snps = pd.Series(chr_snps, name="chr_snps")
    return chr_snps

def plot_fcn(rtt_rel,name,type="clade"):
    fig, ax = plt.subplots(figsize=(8,8))
    ax.axline(xy1=(0, 0), slope=1, alpha=0.8)
    sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
    if type == "total":
        sns.scatterplot(data=rtt_rel, y="snps", x="dcj-indel", ax=ax, hue="subcommunity",size="core_size", alpha=0.8, palette=colours)
    elif type == "subcomm":
        sns.scatterplot(data=rtt_rel, y="snps", x="dcj-indel", ax=ax, hue="clade", alpha=0.8)
    elif type == "st":
        sns.scatterplot(data=rtt_rel, y="snps", x="dcj-indel", ax=ax, hue="st",size="core_size", alpha=0.8)
    else:
        sns.scatterplot(data=rtt_rel, y="snps", x="dcj-indel", ax=ax)
    plt.ylabel("snps", fontsize=12)
    plt.xlabel("dcj-indel", fontsize=12)
    ax.set_title("Parsimonous root-to-tip SNPs vs DCJ-Indel", fontsize=14)
    plt.savefig(f"snpsvsdcj/{name}.png")

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

cluster_path = "/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj"
clusters = [os.path.basename(el).replace('_edists.txt','') for el in glob.glob(f"{cluster_path}/*_edists.txt")]
clusters = sorted(clusters)
total = pd.DataFrame(columns=["dcj-indel", "snps"])

for cluster in clusters:

    snp_tree = Tree("/home/daria/Documents/projects/ABC/host_presence/host_tree/norm_phylogeny.nwk")

    meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

    pangraph = f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json"
    graph = pp.Pangraph.from_json(pangraph)
    core_aln = graph.core_genome_alignment()
    A = np.array(core_aln)
    core_size = len(A[0])/1000

    split = cluster.split("_")
    name = '_'.join(['c',split[3],'sc',split[5]])
    clade = '_'.join([split[0],split[1]])
    st = split[0]

    edistpath = f"/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj/{cluster}_edists.txt"
    dists = edgewise_distances(edistpath)
    outpath = f"/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/rtt_dcj/{cluster}_rtt_dcj.tsv"
    tree_labelled = Tree(f"/home/daria/Documents/projects/ABC/pangraph_workflow/fitch/relabelled_trees/{cluster}.nw", format=8)
    indel_rtt = root_to_tip_scores(tree_labelled,dists,outpath)

    snp_per_site = pd.read_csv(f"/home/daria/Documents/projects/ABC/pangraph_workflow/fitch/scores/filtered/{cluster}_rtt_snp_scores.tsv", sep="\t", index_col=0)

    time = {}
    filepath = f"/home/daria/Documents/projects/ABC/clades/trees/{cluster}.nw"
    tree = Tree(filepath, format=1)
    root = tree.get_tree_root()
    for leaf in tree.iter_leaves():
        time[leaf.name] = tree.get_distance(root,leaf)

    time = pd.Series(time, name="time")

    chr_snps = chromosomal_snps(split[0],split[1],'_'.join(['community',split[3],'subcommunity',split[5]]),snp_tree)
    chr_snps = chr_snps[time.index]

    indel_rtt = pd.Series(indel_rtt)
    snp_rtt = snp_per_site.sum(axis=1, numeric_only=True)

    rtt = pd.concat([indel_rtt,snp_rtt], axis=1)
    rtt = pd.concat([rtt,time], axis=1)
    rtt = pd.concat([rtt,chr_snps], axis=1)
    rtt.columns = ["dcj-indel", "snps", "time", "chr_snps"]
    rtt.fillna(0, inplace=True)
    rtt_rel_indels = rtt["dcj-indel"]/rtt["time"]
    rtt_rel_indels.name = "dcj-indel"
    rtt_rel_snps = rtt["snps"]/rtt["time"]
    rtt_rel_snps.name = "snps"
    rtt_rel = pd.DataFrame([rtt_rel_indels,rtt_rel_snps]).transpose()
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([name for el in rtt_rel.index], columns=["subcommunity"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([core_size for el in rtt_rel.index], columns=["core_size"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([clade for el in rtt_rel.index], columns=["clade"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([st for el in rtt_rel.index], columns=["st"],index=rtt_rel.index)],axis=1)

    #plot_fcn(rtt_rel, cluster)

    total = pd.concat([total,rtt_rel])

'''
plot_fcn(total, "total", "total")
plot_fcn(total[(total["snps"]<0.5)&(total["dcj-indel"]<0.5)], "total_zoom_in", "total")
plot_fcn(total[total["subcommunity"]=="c_0_sc_503"], "c_0_sc_503","subcomm")
plot_fcn(total[total["subcommunity"]=="c_0_sc_501"], "c_0_sc_501","subcomm")
'''

plot_fcn(total, "st_wise", "st")
plot_fcn(total[(total["snps"]<0.5)&(total["dcj-indel"]<0.5)], "st_zoom_in", "st")