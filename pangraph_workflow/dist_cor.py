#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 16:41:24 2023

@author: daria
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

path_anno = "/home/daria/Documents/projects/ABC/goc/no_really_all_plasmids_distances.tsv"
path_align = "/home/daria/Documents/projects/ABC/pangraph_workflow/matching_output/pangraph_dcj_dists/pangraph_dcj.dist"
dist_align = pd.read_csv(path_align, sep='\t')
dist_anno = pd.read_csv(path_anno, sep='\t')
values_align = dist_align["dist"].to_numpy()
values_anno = []
for row in dist_align.index:
    plasmid_1 = dist_align.loc[row,"plasmid_1"]
    plasmid_2 = dist_align.loc[row,"plasmid_2"]
    value = dist_anno[(dist_anno["plasmid_1"]==plasmid_1)&(dist_anno["plasmid_2"]==plasmid_2)]["distance"].values[0] if ((dist_anno["plasmid_1"]==plasmid_1)&(dist_anno["plasmid_2"]==plasmid_2)).any() else dist_anno[(dist_anno["plasmid_2"]==plasmid_1)&(dist_anno["plasmid_1"]==plasmid_2)]["distance"].values[0]
    values_anno.append(value)
dists = pd.DataFrame(data={'pling':values_anno,'pangraph':values_align})
fig, ax = plt.subplots(figsize=(8,8))
ax.axline(xy1=(0, 0), slope=1, alpha=0.5)
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5, "fontsize":12})
sns.scatterplot(data=dists, x="pangraph", y="pling", ax=ax, alpha=0.1)
ax.set_xticks([i*2 for i in range(30)])
ax.set_yticks([i*2 for i in range(30)])
ax.set_xlim([-1,15])
ax.set_ylim([-1,15])
plt.xlabel("pangraph", fontsize=12)
plt.ylabel("pling", fontsize=12)
ax.set_title("Comparison of DCJ-Indel distances", fontsize=14)
plt.show()
print(dists['pangraph'].corr(dists['pling']))


sns.histplot(data=dists['pangraph'], discrete=True)

#correlation resutl: 0.9392253257275215
