import pickle
from typing import cast
from plasnet.communities import Communities
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mode, mean, median
import itertools

def get_dists():
    dcj={}
    with open("no_really_all_plasmids_distances.tsv", "r") as f:
        next(f)
        for line in f:
            plasmid_1, plasmid_2, dist = line.strip().split('\t')
            dcj[(plasmid_1,plasmid_2)] = int(dist)
            dcj[(plasmid_2,plasmid_1)] = int(dist)
    return dcj

communities_pickle = "/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/communities.pkl"
communities = cast(Communities, Communities.load(communities_pickle))

communities = communities.get_graphs_sorted_by_size()

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
hubs = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/hub_plasmids.csv", sep="\t")["hub_plasmids"].to_list()

dists = get_dists()

boundary_dcj = {}
median_boundary_dcj = []
median_internal_dcj = []
for community in communities:
    num=community.label.split("_")[1]
    subcoms=[com for com in list(set(typing["type"])) if len(typing[typing["type"]==com])>1 and com.split("_")[1]==num]
    singletons = [typing[typing["type"]==com]["plasmid"].values[0] for com in list(set(typing["type"])) if len(typing[typing["type"]==com])==1 and com.split("_")[1]==num]
    subcoms=sorted(subcoms, key=lambda com: len(typing[typing["type"]==com]))
    all_plasmids = typing[typing["type"].isin(subcoms)]["plasmid"].to_list() + singletons

    processed_edges = []
    for com in subcoms:
        boundary_dcj[com]=[]
        plasmids = list(typing[typing["type"]==com]["plasmid"])
        if len(plasmids)>4 and com not in ["community_0_subcommunity_377","community_0_subcommunity_378","community_0_subcommunity_379","community_0_subcommunity_380",]:
            boundary = nx.edge_boundary(community, plasmids)
            for edge in boundary:
                if not edge in processed_edges or not (edge[1],edge[0]) in processed_edges:
                    processed_edges.append(edge)
                    if edge[1] not in hubs:
                        boundary_dcj[com].append(dists[edge])
            fig, ax = plt.subplots()
            if len(boundary_dcj[com])>0:
                median_boundary_dcj.append(median(boundary_dcj[com]))
                print(median(boundary_dcj[com]), mode(boundary_dcj[com]), mean(boundary_dcj[com]))
                sns.histplot(data=boundary_dcj,x=com, ax=ax, alpha=0.5, discrete=True)
            dcj= [dists[pair] for pair in itertools.combinations(plasmids,2)]
            print(median(dcj),mode(dcj),mean(dcj))
            median_internal_dcj.append(median(dcj))
            sns.histplot(data=dcj, ax=ax, alpha=0.5, discrete=True)
            plt.savefig(f"boundary/boundary_dist_{com}.png")
            plt.close(fig)
            

fig, ax = plt.subplots()
sns.histplot(data=median_boundary_dcj, ax=ax, discrete=True)
plt.savefig(f"boundary/dist_median.png")
plt.close(fig)
fig, ax = plt.subplots()
datas = pd.DataFrame({"median":median_boundary_dcj+median_internal_dcj, "location":["boundary" for i in range(len(median_boundary_dcj))]+["internal" for i in range(len(median_internal_dcj))]})
sns.swarmplot(data=datas, y="median", x="location", zorder=0, ax=ax)
sns.boxplot(data=datas, y="median", x="location", fill=False, ax=ax)
plt.savefig(f"boundary/internal_vs_boundary.png")
