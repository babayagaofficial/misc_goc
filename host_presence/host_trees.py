from ete3 import Tree, TreeStyle, NodeStyle, random_color, CircleFace, TextFace, faces
import pandas as pd
import numpy as np

#filepath = '/home/daria/Documents/projects/ABC/host_presence/host_tree/norm_phylogeny.nwk'

st = 69
filepath = f'/home/daria/Documents/projects/ABC/host_presence/host_tree/ST{st}_100000000_1000_arc_BD.tre'

meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv", usecols=["Genome_ID", "ST"])

hosts = list(set(meta[meta["ST"]==st]["Genome_ID"].values))
hosts.remove('30134_6#317') #this is not in the dated tree for ST69 for some reason

'''
for i in range(len(hosts)):
	
	hosts[i] = hosts[i].replace('#','_')
'''

presence = pd.read_csv("/home/daria/Documents/projects/ABC/host_presence/presence_per_host.tsv", sep='\t', true_values=['1'], false_values=['0'], index_col = 0)
presence = presence.loc[hosts,:]

subcomms = [subcom for subcom in presence.columns if presence[subcom].any()]


t = Tree(filepath, format=1)

t.prune(hosts)


colour_codes = []
colours = {}
for i in range(256):
    for j in range(256):
        for k in range(256):
            colour_codes.append([i,j,k])            
rng = np.random.default_rng()
all_colours = rng.choice(colour_codes, len(subcomms),replace=False)
for i in range(len(subcomms)):
    hexcode = [hex(all_colours[i][j]).replace('0x','') for j in range(3)]
    colours[subcomms[i]] = '#'+hexcode[0]+hexcode[1]+hexcode[2]

def layout(n):
    if n.is_leaf():
        for community in [subcom for subcom in presence.columns if presence.loc[n.name,subcom]]:
            plasmid = CircleFace(radius=1.5, color=colours[community], style="circle")
            faces.add_face_to_node(plasmid, n, 1)
        nstyle = NodeStyle()
        nstyle['size'] = 0
        nstyle['hz_line_width'] = 1
        nstyle['vt_line_width'] = 1
        n.set_style(nstyle)
    else:
        nstyle=NodeStyle()
        nstyle['size']=0
        nstyle['hz_line_width'] = 1
        nstyle['vt_line_width'] = 1
        n.set_style(nstyle)


ts = TreeStyle()
ts.layout_fn = layout
ts.show_leaf_name = False
#ts.scale = 5000
ts.branch_vertical_margin = 2.5

for host in sorted(colours.keys()):
	ts.legend.add_face(CircleFace(2, colours[host]), column=0)
	ts.legend.add_face(TextFace(f" {host}", fsize=10), column=1)

t.render(f'{st}_dated.pdf', tree_style=ts)
