import pandas as pd
from Bio import SeqIO
import glob
import os
from math import floor
from statistics import mode, mean, median
import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt

cores = {"tool":[], "core_size_%":[], "cluster_size":[]}

ggcaller = "ggcaller"
cluster_path = "/home/daria/Documents/projects/ABC/goc/hub_density/density_04_type/lists"
clusters = [os.path.basename(el).replace('.txt','').replace("_list", '') for el in glob.glob(f"{cluster_path}/*.txt")]
metadata = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv", usecols = ["Plasmid_ID", "Length"])
core_thresh = 0.8
modal_core_sizes = []
for cluster in clusters:
    try:
        rel_core_size = []
        dir = f"{ggcaller}/{cluster}"
        genes = pd.read_csv(f"{dir}/gene_presence_absence_roary.csv")
        fastas = [el[0] for el in pd.read_csv(f"{cluster_path}/{cluster}.txt", header=None).values]
        num_isolates = len(fastas)
        cores["cluster_size"].append(num_isolates)
        num_genes = len(genes[genes["No. isolates"]>=num_isolates*core_thresh])
        avg_core_length = genes[genes["No. isolates"]>=num_isolates*core_thresh]["Avg group size nuc"].sum() 
        for fasta in fastas:
            name = os.path.basename(fasta).replace('.fna','')
            length = metadata[metadata["Plasmid_ID"]==name]["Length"].values[0]
            ids = genes[genes["No. isolates"]==num_isolates][name].to_list()
            rel_core_size.append(floor((avg_core_length/length)*100))
        print(f"{cluster}:{median(rel_core_size)}, {num_genes}")
        cores["core_size_%"].append(median(rel_core_size))
        cores["tool"].append("pling")
    except:
        pass

cores_df = pd.DataFrame.from_dict(cores)
fig, (ax, ax2) = plt.subplots(1,2)
fig.tight_layout(pad=5.0)
fig.set_figwidth(15)
sns.scatterplot(data=cores_df, x="core_size_%", y="cluster_size", ax=ax2)
ax2.set_xlim(-2.5,100)
ax.set_title("Distribution of median relative core genome size")
sns.swarmplot(data=cores_df, y="core_size_%", x="tool", zorder=0, ax=ax)
sns.boxplot(data=cores_df, y="core_size_%", x="tool", fill=False, ax=ax)
ax.set_ylim(-2.5,100)
ax2.set_title("Median relative core genome size vs cluster size")
plt.savefig(f"both_plots_{core_thresh}%.png")

