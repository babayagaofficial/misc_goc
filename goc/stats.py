import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def clusterwise_histograms(feature, clusters, metadata, discrete):
    Path(feature).mkdir(exist_ok=True)
    for cluster in clusters:
        plasmids = [os.path.splitext(os.path.basename(el[0]))[0] for el in pd.read_csv(f"{cluster_path}/{cluster}.txt", header=None).values]
        features = metadata[metadata["Plasmid_ID"].isin(plasmids)][feature].tolist()
        sns.displot(data=[el for el in features], discrete=discrete)
        plt.savefig(f"{feature}/{cluster}.pdf")

def samples_per_year(metadata):
    genomes = metadata[["Genome_ID", "Isolation_year"]]
    genomes = genomes.drop_duplicates()
    sns.displot(data=genomes, x = "Isolation_year", discrete=True)
    plt.savefig(f"Isolation_year/samples_per_year.pdf")

def multicopy_genes(clusters):
    Path("multicopy").mkdir(exist_ok=True)
    for cluster in clusters:
        try:
            dir = f"ggcaller/ggcaller/{cluster}"
            genes = pd.read_csv(f"{dir}/gene_presence_absence_roary.csv")
            sns.histplot(data=genes[genes["Avg sequences per isolate"]>1], x = "Avg sequences per isolate")
            plt.savefig(f"multicopy/{cluster}.pdf")
            plt.close()
        except:
            pass

cluster_path = "/home/daria/Documents/projects/ABC/goc/lists"
clusters = [os.path.basename(el).replace('.txt','').replace("_list", '') for el in glob.glob(f"{cluster_path}/*.txt")]

feature = "Isolation_year"
discrete = True
metadata = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")
clusterwise_histograms(feature,clusters,metadata,discrete)

#sns.displot(data=metadata, x = feature, discrete=True)
#plt.savefig(f"{feature}/all_subcommunities.pdf")

#samples_per_year(metadata)

#sns.scatterplot(data=metadata, x="Isolation_year", y="Length")
#plt.savefig("size_over_time.pdf")

#multicopy_genes(clusters)
