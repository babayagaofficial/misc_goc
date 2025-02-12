from ete3 import Tree
import seaborn as sns
import pandas as pd
import itertools
import matplotlib.pyplot as plt


filepath = '/home/daria/Documents/projects/ABC/host_presence/host_tree/norm_phylogeny.nwk'

t = Tree(filepath)


def get_dists():
    dcj={}
    with open("no_really_all_plasmids_distances.tsv", "r") as f:
        next(f)
        for line in f:
            plasmid_1, plasmid_2, dist = line.strip().split('\t')
            dcj[(plasmid_1,plasmid_2)] = int(dist)
            dcj[(plasmid_2,plasmid_1)] = int(dist)
    return dcj

def get_snps(snps, chr_1, chr_2):
    if (chr_1,chr_2) in snps.keys():
        return snps[(chr_1,chr_2)]
    elif (chr_2, chr_1) in snps.keys():
        return snps[(chr_2,chr_1)]
    else:
        if chr_1==chr_2:
            snps[(chr_1,chr_2)] = 0
        else:
            snps[(chr_1,chr_2)]=t.get_distance(chr_1,chr_2)
        return snps[(chr_1,chr_2)]

          

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")

meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

cluster_specs = pd.read_csv("/home/daria/Documents/projects/ABC/goc/cluster_specs.tsv", sep="\t")

big_subcomms = [subcomm for subcomm in list(set(typing["type"].to_list())) if len(typing[typing["type"]==subcomm])>45]
big_subcomms = sorted(big_subcomms)

dcj = get_dists()

snps = dict()

correlations = {}

for subcomm in big_subcomms:
    hosts = cluster_specs[cluster_specs["type"]==subcomm]["host"].values[0].split(",")
    plasmids = typing[typing["type"]==subcomm]["plasmid"]
    for host in hosts:
        com_dcj = []
        com_snp = []
        for plasmid_1, plasmid_2 in itertools.combinations(plasmids,2):
            if ((meta["Plasmid_ID"]==plasmid_1) & (meta["ST"] == int(host))).any() and ((meta["Plasmid_ID"]==plasmid_2) & (meta["ST"] == int(host))).any():
                chr_1 = meta[meta["Plasmid_ID"]==plasmid_1]["Genome_ID"].values[0].replace("#","_")
                chr_2 = meta[meta["Plasmid_ID"]==plasmid_2]["Genome_ID"].values[0].replace("#","_")
                com_dcj.append(dcj[(plasmid_1,plasmid_2)])
                com_snp.append(get_snps(snps,chr_1,chr_2))
        dists = pd.DataFrame(data={'chr_dists':com_snp,'DCJ':com_dcj})
        correlations[(subcomm,host)] = dists['chr_dists'].corr(dists['DCJ'])
        if not pd.isna(correlations[(subcomm,host)]):
            fig, ax = plt.subplots(figsize=(8,8))
            sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
            sns.scatterplot(data=dists, x="DCJ", y="chr_dists", ax=ax, alpha=1)
            plt.savefig(f"correlations/{subcomm}_{host}.png")
            plt.close(fig)
        
        print(correlations[(subcomm,host)])

with open("correlations/correlation_per_st.tsv", "w") as f:
    for key in correlations.keys():
        f.write(f"{key[0]}\t{key[1]}\t{correlations[key]}\n")