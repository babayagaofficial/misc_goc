import pandas as pd

metadata = pd.read_csv("media-1.csv")

plasmids = metadata[(metadata["Circularity_Status"]=="circular") & (metadata["Length"]>10000)]["Plasmid_ID"].values

with open("cleaner_input.txt", "w") as f:
    for plasmid in plasmids:
        f.write(f"/nfs/research/zi/daria/game_of_clones/fastas/{plasmid}.fna\n")
