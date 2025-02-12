import pandas as pd


clusters_df = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")

clusters = list(set(clusters_df["type"].values))

for cluster in clusters:
    if len(clusters_df[clusters_df["type"]==cluster])>9:
        with open(f"lists/{cluster}.txt", "w") as f:
            for name in clusters_df[clusters_df["type"]==cluster]["plasmid"].values:
                f.write(f"/nfs/research/zi/daria/game_of_clones/fastas/{name}.fna\n")
