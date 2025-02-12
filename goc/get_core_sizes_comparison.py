import pandas as pd
from Bio import SeqIO
import glob
import os
from math import floor
from statistics import mode, mean, median
import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt
#import gffutils

tools = [0.8, 0.9, 0.99]

cores = {"threshold":[], "core_size_%":[], "cluster_size":[]}

ggcaller = "ggcaller"
cluster_path = "/home/daria/Documents/projects/ABC/goc/hub_density/density_04_type/lists"
clusters = [os.path.basename(el).replace('.txt','').replace("_list", '') for el in glob.glob(f"{cluster_path}/*.txt")]
metadata = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv", usecols = ["Plasmid_ID", "Length"])
for tool in tools:
    try:
        for cluster in clusters:
            rel_core_size = []
            dir = f"{ggcaller}/{cluster}"
            genes = pd.read_csv(f"{dir}/gene_presence_absence_roary.csv")
            fastas = [el[0] for el in pd.read_csv(f"{cluster_path}/{cluster}.txt", header=None).values]
            num_isolates = len(fastas)
            cores["cluster_size"].append(num_isolates)
            avg_core_length = genes[genes["No. isolates"]>=num_isolates*tool]["Avg group size nuc"].sum() #core genes are those present in 100% of the samples; take the average length over all the samples for core gene size
            for fasta in fastas:
                name = os.path.basename(fasta).replace('.fna','')
                length = metadata[metadata["Plasmid_ID"]==name]["Length"].values[0]
                ids = genes[genes["No. isolates"]==num_isolates][name].to_list()
                rel_core_size.append(floor((avg_core_length/length)*100))
            print(f"{cluster}:{median(rel_core_size)}")
            cores["core_size_%"].append(median(rel_core_size))
            cores["threshold"].append(str(tool))
    except:
        pass

cores_df = pd.DataFrame.from_dict(cores)
fig, (ax, ax2) = plt.subplots(1,2)
fig.tight_layout(pad=5.0)
fig.set_figwidth(15)
sns.scatterplot(data=cores_df, x="core_size_%", y="cluster_size", hue="threshold", ax=ax2, alpha=0.5)
ax2.set_xlim(-2.5,100)
ax.set_title("Distribution of median relative core genome size")
sns.swarmplot(data=cores_df, y="core_size_%", x="threshold", hue="threshold", zorder=0, ax=ax)
sns.boxplot(data=cores_df, y="core_size_%", x="threshold", hue="threshold", fill=False, ax=ax)
ax.set_ylim(-2.5,100)
ax2.set_title("Median relative core genome size vs cluster size")
plt.savefig("both_plots.png")
