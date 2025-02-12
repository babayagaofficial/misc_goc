import pandas as pd


typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

with open("exclude_in_503.txt") as f:
    to_exclude = f.read().split("\n")

print(to_exclude)

hosts = list(set(meta["Genome_ID"].to_list()))
names = []

subcomms = list(set(typing["type"].to_list()))
big_subcomms = [subcomm for subcomm in subcomms if len(typing[typing["type"]==subcomm])>19]

presence_absence = {subcomm:[] for subcomm in big_subcomms}

for host in hosts:
    names.append(host)
    plasmids = meta[meta["Genome_ID"]==host]["Plasmid_ID"].to_list()

    for plasmid in to_exclude:
        if plasmid in plasmids:
            plasmids.remove(plasmid)

    for subcomm in big_subcomms:
        abs_presence = len(typing[typing["plasmid"].isin(plasmids) & (typing["type"]==subcomm)]["plasmid"].to_list())
        if abs_presence == 0:
            presence_absence[subcomm].append(0)
        else:
            presence_absence[subcomm].append(1)

presence_df = pd.DataFrame(data=presence_absence, index=names)
presence_df.sort_index(inplace=True)
presence_df.sort_index(axis=1, inplace=True)
presence_df.to_csv("presence_per_host.tsv", sep="\t")
