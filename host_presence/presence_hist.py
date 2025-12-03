import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

sts = [131,95,73,69]
subcomms = list(set(typing["type"].to_list()))
big_subcomms = [subcomm for subcomm in subcomms if len(typing[typing["type"]==subcomm])>11]

subcomm_presence = {}

with open("/home/daria/Documents/projects/ABC/host_presence/exclude_in_503.txt") as f:
    to_exclude = f.read().split("\n")

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


for st in sts:
    plasmids_st = meta[meta["ST"]==st]["Plasmid_ID"].to_list()
    for plasmid in to_exclude:
        if plasmid in plasmids_st:
            plasmids_st.remove(plasmid)
    subcomm_presence[st] = typing[typing["plasmid"].isin(plasmids_st)]["type"]
    subcomm_presence[st] = subcomm_presence[st].apply(lambda x: x if x in big_subcomms else "other")
    subcomm_presence[st] = subcomm_presence[st].apply(lambda x: x.replace("community","c").replace("sub","s"))


fig, ax = plt.subplots(2,2,figsize=(15,15))
for i in range(4):
    print(subcomm_presence[sts[i]])
    if i in [0,1]:
        sns.histplot(data=subcomm_presence[sts[i]],ax=ax[0,i],stat="percent",hue="type")
        ax[0,i].set_title(sts[i])
        ax[0,i].set_ylim(0,60)
        plt.setp(ax[0,i].get_xticklabels(), rotation=30, horizontalalignment='right')
    else:
        sns.histplot(data=subcomm_presence[sts[i]],ax=ax[1,i%2],stat="percent",hue="type")
        ax[1,i%2].set_title(sts[i])
        ax[1,i%2].set_ylim(0,60)
        plt.setp(ax[1,i%2].get_xticklabels(), rotation=30, horizontalalignment='right')

plt.savefig("hist.png")