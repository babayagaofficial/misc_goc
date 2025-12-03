from parse_gap import *
from compare_rates import *
import math


colours = {
    "c_0_sc_503": "#ff28ba",
    "c_0_sc_496": "#1c007d",
    "c_0_sc_502": "#ba10c2",
    "c_0_sc_500": "#00650c",
    "c_0_sc_493": "#49caff",
    "c_0_sc_494": "#ae0069",
    "c_0_sc_501": "#35926d",
    "c_0_sc_499": "#20aa00",
    "c_0_sc_498": "#f7b69a",
    "c_0_sc_485": "#ffaeeb",
    "c_0_sc_490": "#1c5951",
    "c_0_sc_487": "#416dff",
    "c_0_sc_489": "#590000",
    "c_0_sc_486": '#ca2800',
    "c_0_sc_474": '#aeff0c',
    "c_0_sc_453": '#ff316d',
    "c_0_sc_481": '#510039',
    "c_0_sc_470": '#0096a6',
    "c_0_sc_478": '#65008e',
    "c_0_sc_495": '#0431ff',
    "c_0_sc_482": '#31e7ce'
}


def plot_fcn(rtt_rel,rate,name,type="clade"):
    fig, ax = plt.subplots(figsize=(10,8))
    #ax.axline(xy1=(0, 0), slope=1, alpha=0.8)
    sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
    if type == "total":
        sns.scatterplot(data=rtt_rel, y="SNPs", x=rate, ax=ax, hue="subcommunity",size="core_size", alpha=0.8, palette=colours)
    elif type == "subcomm":
        ax.axline(xy1=(0, 0), slope=1, alpha=0.8)
        sns.scatterplot(data=rtt_rel, y="SNPs", x=rate, ax=ax, hue="clade", alpha=0.8)
    elif type == "st":
        sns.scatterplot(data=rtt_rel, y="SNPs", x=rate, ax=ax, hue="st",size="core_size", alpha=0.8)
    else:
        sns.scatterplot(data=rtt_rel, y="SNPs", x=rate, ax=ax)
    plt.ylabel("SNPs", fontsize=12)
    plt.xlabel(rate, fontsize=12)
    plt.axis('square')
    if not type == "clade":
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set_title(f"Parsimonous root-to-tip SNPs vs {rate}", fontsize=14)
    plt.savefig(f"{rate}/{name}.png")


def read_in(cluster):
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
    dists=None
    if cluster in ["st73_cl497_community_0_subcommunity_501", "st95_cl477_community_0_subcommunity_500"]:
        root_to_tip = {}
        for leaf in tree.iter_leaves():
            root_to_tip[leaf.name] = 0
        dcj_rtt = pd.Series(root_to_tip)
    else:
        tree_labelled = Tree(f"/home/daria/Documents/projects/ABC/pangraph_workflow/fitch/relabelled_trees/{cluster}.nw", format=8)
        edistpath = f"/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj/{cluster}_edists.txt"
        dists = edgewise_distances(edistpath)
        outpath = f"/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/rtt_dcj/{cluster}_rtt_dcj.tsv"
        dcj_rtt = pd.Series(root_to_tip_scores(tree_labelled,dists,outpath))

    indel_per_block = pd.read_csv(f"fitch/scores/{cluster}_rtt_scores.tsv", sep="\t", index_col=0)
    snp_per_site = pd.read_csv(f"fitch/scores/filtered/{cluster}_rtt_snp_scores.tsv", sep="\t", index_col=0)

    time = {}
    root = tree.get_tree_root()
    length = {}
    divide_by = {}
    plasmid = {}
    metadata = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")
    bl_count = graph.to_blockcount_df()
    chr_to_plasmid = {chr.name:plasmid for chr in tree.iter_leaves() for plasmid in bl_count.columns if chr.name in plasmid}
    for leaf in tree.iter_leaves():
        time[leaf.name] = tree.get_distance(root,leaf)
        length[leaf.name] = metadata[metadata["Plasmid_ID"]==chr_to_plasmid[leaf.name]]["Length"].values[0]
        divide_by[leaf.name] = time[leaf.name]*length[leaf.name]
        plasmid[leaf.name] = chr_to_plasmid[leaf.name]

    time = pd.Series(time, name="time")
    length = pd.Series(length, name="length")
    divide_by = pd.Series(divide_by, name="divide_by")
    plasmid = pd.Series(plasmid, name="plasmid")

    indel_rtt = indel_per_block.sum(axis=1, numeric_only=True)
    snp_rtt = snp_per_site.sum(axis=1, numeric_only=True)

    rtt = pd.concat([dcj_rtt,indel_rtt], axis=1)
    rtt = pd.concat([rtt,snp_rtt], axis=1)
    rtt = pd.concat([rtt,time], axis=1)
    rtt = pd.concat([rtt,length], axis=1)
    rtt = pd.concat([rtt,divide_by], axis=1)
    rtt.columns = ["DCJ-Indel", "Fitch", "SNPs", "time","length","divide_by"]
    rtt.fillna(0, inplace=True)
    rtt_rel_indels = rtt["Fitch"]/rtt["divide_by"]
    rtt_rel_indels.name = "Fitch"
    rtt_rel_snps = rtt["SNPs"]/rtt["divide_by"]
    rtt_rel_snps.name = "SNPs"
    rtt_rel_dcj = rtt["DCJ-Indel"]/rtt["divide_by"]
    rtt_rel_dcj.name = "DCJ-Indel"
    rtt_rel = pd.DataFrame([rtt_rel_dcj,rtt_rel_indels,rtt_rel_snps]).transpose()
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([name for el in rtt_rel.index], columns=["subcommunity"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([core_size for el in rtt_rel.index], columns=["core_size"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,pd.DataFrame([clade for el in rtt_rel.index], columns=["clade"],index=rtt_rel.index)],axis=1)
    rtt_rel = pd.concat([rtt_rel,plasmid],axis=1)
    rtt = pd.concat([rtt,plasmid],axis=1)
    return rtt_rel, rtt, dists

parse_gap()

cluster_path = "/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj"
clusters = [os.path.basename(el).replace('_edists.txt','') for el in glob.glob(f"{cluster_path}/*_edists.txt")]
clusters = sorted(clusters)

snps = []
indels = []

total = pd.DataFrame(columns=["DCJ-Indel", "Fitch", "SNPs"])
lengths = pd.DataFrame(columns=["DCJ-Indel", "Fitch", "SNPs"])
comparison = {"clade":[],"proportion":[],"similar":[]}
edists = {}

for cluster in clusters:

    rtt_rel,rtt,dists = read_in(cluster)

    more_dcj = (rtt_rel["DCJ-Indel"]>rtt_rel["SNPs"])

    dcj_snp_similar = (abs(rtt_rel["DCJ-Indel"]-rtt_rel["SNPs"])<0.01)

    comparison["proportion"].append(len(more_dcj[more_dcj==True])/len(more_dcj))
    comparison["clade"].append(cluster)
    comparison["similar"].append(len(dcj_snp_similar[dcj_snp_similar==True])/len(dcj_snp_similar))

    if not cluster in ["st73_cl497_community_0_subcommunity_501", "st95_cl477_community_0_subcommunity_500"]:
        for key in dists.keys():
            if key[0]=="I":
                edists[f"{key}_{cluster}"] = dists[key]
            else:
                edists[key] = dists[key]
    plot_fcn(rtt_rel, "DCJ-Indel", cluster)

    total = pd.concat([total,rtt_rel])
    lengths = pd.concat([lengths,rtt])

print(total[["Fitch", "DCJ-Indel", "SNPs"]].corr())

total.to_csv("rates_per_tip.tsv", sep="\t")
lengths.to_csv("lengths_per_tip.tsv", sep="\t")


fig, ax = plt.subplots(figsize=(10,8))
ax.axline(xy1=(0, 0), slope=1, alpha=0.8)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
print(total)
sns.scatterplot(data=total, y="Fitch", x="DCJ-Indel", ax=ax, hue="subcommunity", alpha=0.8, palette=colours)
plt.axis('square')
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
ax.set_title(f"Parsimonous root-to-tip per year per base rates for Fitch indels vs DCJ-Indel", fontsize=14)
plt.savefig("fitchvsdcj.png")




df = total.melt(id_vars=['subcommunity','core_size','clade'], value_vars=['DCJ-Indel','Fitch','SNPs'], var_name='type', value_name='rates')
df = df[['subcommunity','type','rates']]
types = sorted(colours.keys())
categories = ['DCJ-Indel','Fitch','SNPs']

# Plotting setup
fig, ax = plt.subplots(figsize=(10, 8))

# Widths and positions
n_categories = len(categories)
n_types = len(types)
box_height = 0.2
y = np.arange(len(types))


# Plot each category box within each type group
for i, t in enumerate(types):
    for j, c in enumerate(categories):
        values = df[(df['subcommunity'] == t) & (df['type'] == c)]['rates']
        ypos = y[i] + (j - (len(categories) - 1) / 2) * box_height
        ax.boxplot(values, positions=[ypos], widths=box_height,
                    vert=False,  # horizontal orientation
                    patch_artist=True,
                    boxprops=dict(facecolor=colours[t], color='black'),
                    medianprops=dict(color='black'),
                    whiskerprops=dict(color='black'),
                    capprops=dict(color='black'),
                    flierprops=dict(markerfacecolor=colours[t], marker='o',
                                    markersize=5, linestyle='none', alpha=0.6))

# Labeling
ax.set_yticks(y)
ax.set_yticklabels(types)
ax.set_ylabel('subcommunity')
ax.set_xlabel('rates')
ax.set_title('Distributions of SNP, Fitch indel, and DCJ-Indel per year per base rates per subcommunity')

# Custom legend
handles = [plt.Line2D([0], [0], color=colours[t], lw=10) for t in types]
ax.legend(handles, types, title='subcommunity',
          loc='upper center', bbox_to_anchor=(0.5, -0.12),
          ncol=math.ceil(len(types)/4))
plt.tight_layout()
plt.subplots_adjust(bottom=0.3)

plt.savefig("boxplots_corrected.png")


fig, ax = plt.subplots(figsize=(10,8))
sns.histplot(data=edists, discrete=True)
plt.xlabel("child-parent DCJ-Indel distance", fontsize=14)
plt.ylabel("count", fontsize=14)
ax.set_title("Distribution of DCJ-Indel distances between parent and child nodes", fontsize=18)
plt.savefig(f"edist_hist.png")


comparison = pd.DataFrame(data=comparison)
fig, ax = plt.subplots(figsize=(10,8))
sns.histplot(data=comparison, x="similar", ax=ax)
plt.xlabel("proportion", fontsize=10)
ax.set_title("Histogram of the proportion of plasmids which have similar\nDCJ-Indel rates and SNP rates per clade", fontsize=14)
plt.savefig(f"similar_dcj_snps.png")

dcj = lengths["DCJ-Indel"].rename("DCJ-Indel")
snps = lengths["SNPs"].rename("SNPs")
fitch = lengths["Fitch"].rename("Fitch")
branch_lengths = pd.concat([dcj,snps, fitch], keys=["DCJ-Indel", "SNPs", "Fitch"]).reset_index(level=0)
branch_lengths.columns = ["type", "branch_length"]
fig, ax = plt.subplots(figsize=(10,8))
sns.histplot(data=branch_lengths, x="branch_length", hue="type", element="step", discrete=True, ax=ax)
plt.xlabel("root-to-tip distances", fontsize=14)
plt.ylabel("count", fontsize=14)
ax.set_title("Distribution of root-to-tip distances for SNPs, Fitch indels, and DCJ-Indel", fontsize=18)
plt.savefig(f"DCJ-Indel/branch_lengths_snps.png")
#fig, ax = plt.subplots(figsize=(10,8))

sns.histplot(data=lengths, x="Fitch", discrete=True,ax=ax)
plt.savefig(f"Fitch/branch_lengths.png")


plot_fcn_compare(total, "total", colours, type="total")
plot_fcn(total, "DCJ-Indel", "total", "total")
plot_fcn(total, "Fitch", "total", "total")
sc_503 = total[total["subcommunity"]=="c_0_sc_503"]
plot_fcn(sc_503, "DCJ-Indel", "c_0_sc_503", "subcomm")
plot_fcn(sc_503, "Fitch", "c_0_sc_503", "subcomm")
sc_501 = total[total["subcommunity"]=="c_0_sc_501"]
plot_fcn(sc_501, "DCJ-Indel", "c_0_sc_501", "subcomm")
plot_fcn(sc_501, "Fitch", "c_0_sc_501", "subcomm")

rate="DCJ-Indel"
fig, ax = plt.subplots(figsize=(10,8))
#ax.axline(xy1=(0, 0), slope=1, alpha=0.8)
#ax[1].axline(xy1=(0, 0), slope=1, alpha=0.8)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
sns.scatterplot(data=total, y="SNPs", x=rate, ax=ax, hue="subcommunity",size="core_size", alpha=0.8, palette=colours)
#sns.scatterplot(data=total[(total["SNPs"]<0.5)&(total["DCJ-Indel"]<0.5)], y="SNPs", x=rate, ax=ax[1], hue="subcommunity",size="core_size", alpha=0.8, palette=colours)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

axins = ax.inset_axes(
    [0.4, 0.4, 0.6, 0.6],
    xlim=(-0.0000001, 0.5*0.00001), ylim=(-0.0000001, 0.5*0.00001), xticklabels=[], yticklabels=[])
sub = sns.scatterplot(data=total[(total["SNPs"]<0.5*0.00001)&(total["DCJ-Indel"]<0.5*0.00001)], y="SNPs", x=rate, ax=axins, hue="subcommunity",size="core_size", alpha=0.8, palette=colours, legend=False)
sub.set(xlabel=None)
sub.set(ylabel=None)

#ax.set_aspect('equal', adjustable='box')
plt.axis('square')
ax.indicate_inset_zoom(axins, edgecolor="black")
plt.ylabel("SNPs", fontsize=12)
plt.xlabel(rate, fontsize=12)
fig.suptitle(f"Parsimonous root-to-tip per year per base rates for SNPs vs {rate}", fontsize=14)
plt.savefig(f"{rate}/zoom.png")
