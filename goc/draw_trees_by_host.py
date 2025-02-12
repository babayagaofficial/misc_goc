from ete3 import Tree, TreeStyle, NodeStyle, random_color, CircleFace, TextFace, faces, AttrFace, RectFace
import argparse
import pandas as pd
import numpy as np
import math
import glob
import os

years = pd.read_csv('/home/daria/Documents/projects/ABC/goc/media-1.csv', usecols=['Plasmid_ID','ST', "Isolation_year"])
clusters = [el for el in glob.glob(f"submatrices/submatrices/trees/*")]
#cluster_path = "/home/daria/Documents/projects/ABC/goc/lists"
#clusters = [os.path.basename(el).replace('.txt','').replace("_list", '') for el in glob.glob(f"{cluster_path}/*.txt")]
#clusters.remove("community_0_subcommunity_475")
#clusters.remove("community_0_subcommunity_481")

gradient = dict()
year_vals = list(set(years['Isolation_year'].values))
hosts = list(set(years['ST'].values))

k = math.ceil(len(year_vals)/2)
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


colour_codes = []
colours = {}
for i in range(256):
	for j in range(256):
		for k in range(256):
			colour_codes.append([i,j,k])
rng = np.random.default_rng(seed=42)
all_colours = rng.choice(colour_codes, len(hosts),replace=False)
for i in range(len(hosts)):
	hexcode = [hex(all_colours[i][j]).replace('0x','') for j in range(3)]
	colours[hosts[i]] = '#'+hexcode[0]+hexcode[1]+hexcode[2]

def layout(n):
	if n.is_leaf():
		year = years[years['Plasmid_ID']==n.name]['Isolation_year']
		host = years[years['Plasmid_ID']==n.name]['ST']
		year = int(year)
		host = int(host)
		hostface = RectFace(width=4,height=4, fgcolor=gradient[year], bgcolor=gradient[year])
		faces.add_face_to_node(hostface, n, 1)
		nstyle = NodeStyle()
		nstyle['shape'] = 'circle'
		nstyle['size'] = 8
		nstyle['fgcolor'] = colours[host]
		present_hosts.add(host)
		n.set_style(nstyle)
	else:
		nstyle=NodeStyle()
		nstyle['size']=0
		n.set_style(nstyle)

for cluster in clusters:
	present_hosts=set()
	print(cluster)

	#filepath = f"ggcaller/ggcaller/{cluster}/pangenome_NJ.nwk"

	#filepath = f"/home/daria/Documents/projects/ABC/goc/parsnp/parsnp/{cluster}/parsnp.nw"
	
	filepath = cluster

	t = Tree(filepath, format=0, quoted_node_names=True)

	ts = TreeStyle()
	ts.layout_fn = layout
	ts.show_leaf_name = True
	ts.scale = 250
	ts.branch_vertical_margin = 5


	for host in present_hosts:
		ts.legend.add_face(CircleFace(10, colours[host]), column=0)
		ts.legend.add_face(TextFace(f" {host}", fsize=12), column=1)
	name = cluster.replace("submatrices/submatrices/trees","").replace(".tree","")
	t.render(f"dcj_trees/{name}.pdf", tree_style=ts)
