import pandas as pd
import glob
import os
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np


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

metadata = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")
cluster_path = "/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj"
clusters = [os.path.basename(el).replace('_edists.txt','') for el in glob.glob(f"{cluster_path}/*_edists.txt")]
clusters = sorted(clusters)
is_els = {}
multi_is_els = {}
bacteriocins = {}
plasmids = []
rates = pd.read_csv("/home/daria/Documents/projects/ABC/pangraph_workflow/rates_per_tip.tsv", sep="\t", index_col="plasmid")
lengths = pd.read_csv("/home/daria/Documents/projects/ABC/pangraph_workflow/lengths_per_tip.tsv", sep="\t", index_col="plasmid")

'''
for cluster in clusters:
    clade_plasmids = [el.split("/")[-1] for el in glob.glob(f"isescan/{cluster}/*")]
    plasmids = plasmids + clade_plasmids
    for plasmid in clade_plasmids:
        try:
            csv = pd.read_csv(f"isescan/{cluster}/{plasmid}/fastas/{plasmid}.fna.csv")
            total = csv["ncopy4is"].sum()
            multi_count = csv[csv["ncopy4is"]>1]["ncopy4is"].sum()
        except:
            total = 0
            multi_count = 0
        is_els[plasmid] = total
        multi_is_els[plasmid] = multi_count

for cluster in clusters:
    genes = pd.read_csv(f"/home/daria/Documents/projects/ABC/pangraph_workflow/ggcaller/{cluster}/gene_presence_absence.csv")
    clade_plasmids = [el for el in list(genes.columns) if el not in ["Gene", "Non-unique Gene name", "Annotation"]]
    plasmids = plasmids + clade_plasmids
    split = cluster.split("_")
    subcomm = '_'.join(split[2:])
    #genes = pd.read_csv(f"/home/daria/Documents/projects/ABC/goc/ggcaller/ggcaller/{subcomm}/gene_presence_absence.csv")
    for row in genes.index:
        annotation = str(genes.loc[row,"Annotation"])
        name = genes.loc[row,"Gene"]
        if "cvaC" in annotation or "cvac" in annotation:
            if "Colicin" in annotation:
                row = genes[genes["Gene"]==name].copy()
                row.dropna(how='all', axis=1, inplace=True)
                present_in = list(row.columns)
                present_in = [el for el in present_in if el in clade_plasmids]
                for plasmid in present_in:
                    if plasmid in bacteriocins.keys():
                        bacteriocins[plasmid].append("microcin V")
                    else:
                        bacteriocins[plasmid]=["microcin V"]
        elif "cia" in annotation:
            if "Colicin" in annotation:
                row = genes[genes["Gene"]==name].copy()
                row.dropna(how='all', axis=1, inplace=True)
                present_in = list(row.columns)
                present_in = [el for el in present_in if el in clade_plasmids]
                for plasmid in present_in:
                    if plasmid in bacteriocins.keys():
                        bacteriocins[plasmid].append("colicin Ia")
                    else:
                        bacteriocins[plasmid]=["colicin Ia"]

col_carrying = {}
for cluster in clusters:
    clade_plasmids = [el.split("/")[-1] for el in glob.glob(f"isescan/{cluster}/*")]
    plasmids = plasmids + clade_plasmids
    for plasmid in clade_plasmids:
        host = metadata[metadata["Plasmid_ID"]==plasmid]["Genome_ID"].values[0]
        other_types = list(metadata[(metadata["Genome_ID"]==host) & (metadata["Plasmid_ID"]!=plasmid)]["Plasmidfinder"].values)
        for inctype in other_types:
            if "Col" in inctype:
                if plasmid in col_carrying.keys():
                    col_carrying[plasmid].append(True)
                else:
                    col_carrying[plasmid] = [True]


for cluster in clusters:
    genes = pd.read_csv(f"/home/daria/Documents/projects/ABC/pangraph_workflow/ggcaller/{cluster}/gene_presence_absence.csv")
    clade_plasmids = [el for el in list(genes.columns) if el not in ["Gene", "Non-unique Gene name", "Annotation"]]
    plasmids = plasmids + clade_plasmids
    split = cluster.split("_")
    subcomm = '_'.join(split[2:])
    #genes = pd.read_csv(f"/home/daria/Documents/projects/ABC/goc/ggcaller/ggcaller/{subcomm}/gene_presence_absence.csv")
    for row in genes.index:
        annotation = str(genes.loc[row,"Annotation"])
        name = genes.loc[row,"Gene"]
        if "Insertion element" in annotation or "insertion element" in annotation or "Insertion sequence" in annotation or "insertion sequence" in annotation:
            row = genes[genes["Gene"]==name].copy()
            row.dropna(how='all', axis=1, inplace=True)
            present_in = list(row.columns)
            present_in = [el for el in present_in if el in clade_plasmids]
            for plasmid in present_in:
                if plasmid in is_els.keys():
                    is_els[plasmid].append(name)
                else:
                    is_els[plasmid] = [name]


for cluster in clusters:
    clade_plasmids = [el.split("/")[-1] for el in glob.glob(f"isescan/{cluster}/*")]
    plasmids = plasmids + clade_plasmids
    for plasmid in clade_plasmids:
        try:
            with open(f"isescan/{cluster}/{plasmid}/fastas/{plasmid}.fna.sum") as f:
                for line in f:
                    if ".fna" in line:
                        split = line.strip().split(" ")
                        content = [el for el in split if el!=""]
                        perc = content[3]
        except:
            print(cluster, plasmid)
            perc = 0
        is_els[plasmid] = float(perc)

more_dcj = (rates["DCJ-Indel"]>rates["SNPs"])
num_more_dcj = len(rates[more_dcj])
print(more_dcj)
print(len(plasmids), num_more_dcj)
'''

neighbours = {}
cooccurrence = pd.read_csv("/home/daria/Documents/projects/ABC/supp_table_2.csv")
for cluster in clusters:
    clade_plasmids = [el.split("/")[-1] for el in glob.glob(f"isescan/{cluster}/*")]
    plasmids = plasmids + clade_plasmids
    for plasmid in clade_plasmids:
        host = metadata[metadata["Plasmid_ID"]==plasmid]["Genome_ID"].values[0]
        pT = metadata[metadata["Plasmid_ID"]==plasmid]["Plasmid type"].values[0]
        other_types = list(metadata[(metadata["Genome_ID"]==host) & (metadata["Plasmid_ID"]!=plasmid)]["Plasmid type"].values)
        subset = cooccurrence[(cooccurrence["backbone_i"].isin(other_types) & (cooccurrence["backbone_j"]==pT)) | (cooccurrence["backbone_j"].isin(other_types) & (cooccurrence["backbone_i"]==pT))]
        positive = subset[(subset["phi_cat"]=="positive") & (subset["r_cat"]=="positive")]
        associated = positive[["backbone_i","backbone_j"]].values
        if len(associated)>0:
            associated=associated[0]
            associated = [el for el in associated if el!=pT]
            if len(associated)==1:
                neighbours[plasmid] = associated[0]
            else:
                neighbours[plasmid] = " and ".join(list(set(associated)))

        
        
        

print(rates.columns)
rates["co-occuring"] = ["none" for el in plasmids]
#lengths["IS"] = [0 for el in plasmids]
#lengths = pd.concat([lengths,rates["subcommunity"]],axis=1)
print(neighbours)
for key in neighbours.keys():
    if key in rates.index:
        #rates.loc[key, "IS"] = len(is_els[key])
        #lengths.loc[key,"IS"] = len(is_els[key])
        rates.loc[key, "co-occuring"] = neighbours[key]

print(rates)
sns.scatterplot(data=rates, y="DCJ-Indel", x="SNPs", hue="co-occuring", alpha=0.8)
plt.show()
