import glob
import pandas as pd
import math
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, faces

typing = pd.read_csv("/home/daria/Documents/projects/ABC/goc/typing/type_thresh_4/objects/typing.tsv", sep="\t")
meta = pd.read_csv("/home/daria/Documents/projects/ABC/goc/media-1.csv")

with open("/home/daria/Documents/projects/ABC/host_presence/exclude_in_503.txt") as f:
    to_exclude = f.read().split("\n")

gradient = dict()
year_vals = list(set(meta['Isolation_year'].values))

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


def layout(n):
    if n.is_leaf():
        year = meta[meta['Genome_ID']==n.name]['Isolation_year'].values[0]
        textface = TextFace(year, fgcolor=gradient[year])
        faces.add_face_to_node(textface,n,1)
        nstyle = NodeStyle()
        nstyle['shape'] = 'circle'
        nstyle['size'] = 2
        nstyle['fgcolor'] = gradient[year]
        n.set_style(nstyle)
    else:
        nstyle=NodeStyle()
        nstyle['size']=0
        n.set_style(nstyle)
          


for file in glob.glob("/home/daria/Documents/projects/ABC/host_presence/treeseg/st*/community_*_subcommunity_*/clades/*.txt"):
    clade_plasmids = []
    clade_genomes = []
    with open(file, "r") as f:
        genomes = f.read().split("\n")
    name_split = file.split("/")
    subcom = name_split[-3]
    st = name_split[-4]
    clade = name_split[-1].replace(".txt", "")
    plasmids = typing[typing["type"]==subcom]["plasmid"].to_list()
    for plasmid in plasmids:
        genome = meta[meta["Plasmid_ID"]==plasmid]["Genome_ID"].values[0]
        if genome in genomes:
            if not plasmid in to_exclude:
                clade_plasmids.append(plasmid)
                clade_genomes.append(genome)
                
    st_num = st.replace('st',"")
    tree = Tree(f"/home/daria/Documents/projects/ABC/host_presence/host_tree/pruned_dated_{st_num}.nwk")
    try:
        tree.prune(clade_genomes, preserve_branch_length=True)
    except:
          pass
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = True
    ts.branch_vertical_margin = 5
    tree.render(f'trees/{st}_cl{clade}_{subcom}.pdf', tree_style=ts)

    filename = f"lists/{st}_cl{clade}_{subcom}.txt"
    leaves = tree.get_leaf_names()
    f = lambda plasmid: leaves.index(meta[meta["Plasmid_ID"]==plasmid]["Genome_ID"].values[0])
    clade_plasmids.sort(key=f)
    with open(filename, "w") as f:
        for plasmid in clade_plasmids:
            f.write(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna\n")
    



