import pandas as pd
from statistics import mean

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
mob = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

mob_per_subcomm = {"type":[], "host": [], "#plasmids":[], "avg_length":[], "rep_types":[], "relaxase_types":[], "mpf_type":[], "predicted_mobility":[]}

subcomms = list(set(typing["type"].to_list()))

for subcomm in subcomms:
    mob_per_subcomm["type"].append(subcomm)
    plasmids = typing[typing["type"]==subcomm]["plasmid"].to_list()
    hosts = [str(el) for el in list(set(mob[mob["Plasmid_ID"].isin(plasmids)]["ST"].to_list()))]
    mob_per_subcomm["host"].append(",".join(hosts))

    mob_per_subcomm["#plasmids"].append(len(plasmids))
    mob_per_subcomm["avg_length"].append(round(mean(mob[mob["Plasmid_ID"].isin(plasmids)]["Length"].to_list())))


    reps = list(set(mob[mob["Plasmid_ID"].isin(plasmids)]["MOB-typer_Replication"].to_list()))
    relaxases = list(set(mob[mob["Plasmid_ID"].isin(plasmids)]["MOB-typer_Relaxase"].to_list()))
    mpf = list(set(mob[mob["Plasmid_ID"].isin(plasmids)]["MOB-typer_MPF"].to_list()))
    mobility = list(set(mob[mob["Plasmid_ID"].isin(plasmids)]["MOB-typer_Mobility"].to_list()))
    mob_per_subcomm["rep_types"].append(",".join(reps))
    mob_per_subcomm["relaxase_types"].append(",".join(relaxases))
    mob_per_subcomm["mpf_type"].append(",".join(mpf))
    mob_per_subcomm["predicted_mobility"].append(",".join(mobility))

mob_per_subcomm_df = pd.DataFrame(data=mob_per_subcomm)
mob_per_subcomm_df.sort_values(by=["#plasmids"], inplace=True, ascending=False)
mob_per_subcomm_df.to_csv("cluster_specs.tsv", sep="\t", index=False)

