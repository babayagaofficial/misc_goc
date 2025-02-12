import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

sts = list(set(meta["ST"].to_list()))
subcomms = list(set(typing["type"].to_list()))
big_subcomms = [subcomm for subcomm in subcomms if len(typing[typing["type"]==subcomm])>1]

presence_absence = {subcomm:[] for subcomm in big_subcomms}
index_sts = []

for st in sts:
    plasmids_st = meta[meta["ST"]==st]["Plasmid_ID"].to_list()
    if len(plasmids_st)>0:
        index_sts.append(st)
        num_genomes_st = len(list(set(meta[meta["ST"]==st]["Genome_ID"].to_list())))
        for subcomm in big_subcomms:
            abs_presence = len(typing[typing["plasmid"].isin(plasmids_st) & (typing["type"]==subcomm)]["plasmid"].to_list())
            rel_presence = round(abs_presence/num_genomes_st*100)
            presence_absence[subcomm].append(rel_presence)

presence_df = pd.DataFrame(data=presence_absence, index=index_sts)
presence_df.sort_index(inplace=True)
presence_df.sort_index(axis=1, inplace=True)
presence_df.to_csv("presence_per_st.tsv", sep="\t")
sns.heatmap(data=presence_df)
plt.show()
