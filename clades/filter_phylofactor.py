import glob
import pandas as pd
import math
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, faces

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")


with open("/home/daria/Documents/projects/ABC/host_presence/exclude_in_503.txt") as f:
    to_exclude = f.read().split("\n")

          

clades = {"subcomm":[], "st":[], "clade":[], "plasmids":[], "genomes":[]}
for file in glob.glob("/home/daria/Documents/projects/ABC/host_presence/phylofactor/st*/community_*_subcommunity_*/clades/*.txt"):
    clade_plasmids = []
    clade_genomes = []
    with open(file, "r") as f:
        genomes = f.read().split("\n")
    name_split = file.split("/")
    subcom = name_split[-3]
    st = name_split[-4]
    st_num = st.replace('st',"")
    clade = name_split[-1].replace(".txt", "")
    plasmids = typing[typing["type"]==subcom]["plasmid"].to_list()
    for plasmid in plasmids:
        genome = meta[meta["Plasmid_ID"]==plasmid]["Genome_ID"].values[0]
        if genome in genomes:
            if not plasmid in to_exclude:
                clade_plasmids.append(plasmid)
                clade_genomes.append(genome)

    if len(clade_plasmids)>=4 and len(clade_plasmids)/len(genomes)>=0.5:
        clades["subcomm"].append(subcom)
        clades["st"].append(st_num)
        clades["clade"].append(clade)
        clades["plasmids"].append(set(clade_plasmids))
        clades["genomes"].append(clade_genomes)

clades_df = pd.DataFrame(clades)
for i in range(len(clades["subcomm"])):
    plasmids = clades["plasmids"][i]
    subcom = clades["subcomm"][i]
    st = clades["st"][i]
    clade = clades["clade"][i]
    blah = list(set(clades_df[(clades_df["subcomm"]==subcom) & (clades_df["st"]==st)]["clade"].values))
    blah.remove(clade)
    is_contained = False
    for other_clade in blah:
        members = set(clades_df[clades_df["clade"]==other_clade]["plasmids"].values[0])
        if plasmids.issubset(members):
            is_contained = True
    if not is_contained:

        tree = Tree(f"/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_dated_{st}.nwk")
        tree.prune(clades["genomes"][i], preserve_branch_length=True)

        tree.write(format=1, outfile=f"trees/st{st}_cl{clade}_{subcom}.nw")

        filename = f"lists/st{st}_cl{clade}_{subcom}.txt"
        leaves = tree.get_leaf_names()
        f = lambda plasmid: leaves.index(meta[meta["Plasmid_ID"]==plasmid]["Genome_ID"].values[0])
        plasmids = list(plasmids)
        plasmids.sort(key=f)
        with open(filename, "w") as f:
            for plasmid in plasmids:
                f.write(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna\n")
    



