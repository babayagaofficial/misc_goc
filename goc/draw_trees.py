from ete3 import Tree, TreeStyle, NodeStyle
import argparse
import pandas as pd
import numpy as np
import math
import glob
import os

years = pd.read_csv('/home/daria/Documents/projects/ABC/goc/media-1.csv', usecols=['Plasmid_ID','Isolation_year'])
cluster_path = "/home/daria/Documents/projects/ABC/goc/lists"
clusters = [os.path.basename(el).replace('.txt','').replace("_list", '') for el in glob.glob(f"{cluster_path}/*.txt")]

gradient = dict()
year_vals = np.array([i for i in range(2002,2018)])
k = math.ceil(len(year_vals)/2)
print(len(year_vals), k)
for i in range(k+1):
	r = hex(math.floor(i*(255/k)))
	r = r.replace('0x','')
	if len(r) == 1:
		r = '0'+r
	hex_val = '#'+r+'00ff'
	gradient[year_vals[i]] = hex_val
	b = hex(math.floor((k-i+1)*(255/k)))
	b = b.replace('0x','')
	if len(b) ==1:
		b = '0'+b
	hex_val = '#ff00'+b
	gradient[year_vals[k-1+i]] = hex_val
gradient[np.median(year_vals)]='#ff00ff'
print(gradient)

for cluster in clusters:
	try:
		#filepath = f"ggcaller/ggcaller/{cluster}/pangenome_NJ.nwk"
		filepath = f"parsnp/parsnp/{cluster}/parsnp.nw"

		t = Tree(filepath, format=0, quoted_node_names=True)

		ts = TreeStyle()
		ts.show_leaf_name = True
		ts.scale = 250
		ts.branch_vertical_margin = 5

		for n in t.traverse():
			if n.is_leaf():
				year = years[years['Plasmid_ID']==n.name]['Isolation_year']
				if year.isnull().values.any():
					nstyle = NodeStyle()
					nstyle['shape'] = 'circle'
					nstyle['size'] = 2
					nstyle['fgcolor'] = '#000000'
				else:
					year = int(year)
					nstyle = NodeStyle()
					nstyle['shape'] = 'sphere'
					nstyle['size'] = 8
					nstyle['fgcolor'] = gradient[year]
				n.set_style(nstyle)
			else:
				nstyle=NodeStyle()
				nstyle['size']=0
				n.set_style(nstyle)

		t.render(f"trees/{cluster}.pdf", tree_style=ts)
	except:
		pass
