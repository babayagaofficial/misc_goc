import pandas as pd
import numpy as np

subcommunities = "/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv"

subcommunities_df = pd.read_csv(subcommunities, sep='\t')
subcommunities_list = list(subcommunities_df['type'])
count = [subcommunity for subcommunity in set(subcommunities_list) if subcommunities_list.count(subcommunity)>4]

singletons = [subcommunities_df[subcommunities_df["type"]==com]["plasmid"].values[0] for com in set(subcommunities_list) if subcommunities_list.count(com)==1]
print(len(singletons))

assigned = 0
for com in count:
    com_count = len(subcommunities_df[subcommunities_df["type"]==com]["plasmid"].values)
    assigned = assigned + com_count

print(assigned)

#print(count)
print(len(count))
