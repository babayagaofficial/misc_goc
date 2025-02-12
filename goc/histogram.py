import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dists = pd.read_csv("no_really_all_plasmids_distances.tsv", sep="\t")

fig, ax = plt.subplots()
sns.histplot(data=dists,x="distance", ax=ax, discrete=True)
#ax.set_xlim(0,20)

plt.savefig("no_really_dcj_histplots.png")

print(dists["distance"].mode())
print(dists["distance"].median())
print(dists["distance"].mean())
