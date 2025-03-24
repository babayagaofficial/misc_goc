import pypangraph as pp
import itertools
import subprocess
import os
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
import json

def get_adj(neighbour, block):
    if (neighbour[1] and block[1]) or (not neighbour[1] and not block[1]):
        adj = (neighbour[0], 't')
    else:
        adj = (neighbour[0], 'h')
    return adj

'''
graph = pp.Pangraph.from_json("/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/st131_cl416_community_0_subcommunity_503/pangraph.json")

block_stats = graph.to_blockstats_df()
bl_count = graph.to_blockcount_df()
duplicates = list(block_stats[(block_stats["duplicated"]==True) & (block_stats["len"]>200)].index)
core_adjacencies = {}
true_adjacencies = {}

path_dict = graph.to_path_dictionary()

for duplicate in duplicates:
    dup_counts = bl_count.loc[duplicate].transpose()
    genomes = list(dup_counts[dup_counts>0].index)
    aln = graph.blocks[duplicate].to_biopython_alignment()
    core_adjacencies[duplicate]={}
    true_adjacencies[duplicate]={}
    for genome in genomes:
        path = path_dict[genome]
        j=0
        for i in range(len(path)-1):
            if path[i][0]==duplicate:
                if path[i][1]:
                    ladj = get_adj(path[i-1],path[i])
                    radj = get_adj(path[i+1], path[i])
                else:
                    ladj = get_adj(path[i+1],path[i])
                    radj = get_adj(path[i-1], path[i])
                true_adjacencies[duplicate][(genome,j)] = [ladj, radj]
                k = i-1
                not_found = True
                while not_found and k>-1:
                    if block_stats.loc[path[k][0],"core"]==1:
                        not_found = False
                        if path[i][1]:
                            ladj = get_adj(path[k], path[i])
                        else:
                            radj = get_adj(path[k], path[i])
                    k = k -1
                k = i+1
                not_found = True
                while not_found and k<len(path):
                    if block_stats.loc[path[k][0],"core"]==1:
                        not_found = False
                        if path[i][1]:
                            radj = get_adj(path[k], path[i])
                        else:
                            ladj = get_adj(path[k], path[i])
                    k = k +1
                core_adjacencies[duplicate][(genome,j)] = [ladj, radj]
                j+=1
        if path[len(path)-1][0]==duplicate:
                if path[len(path)-1][1]:
                    ladj = get_adj(path[len(path)-2],path[len(path)-1])
                    radj = get_adj(path[0], path[len(path)-1])
                else:
                    ladj = get_adj(path[0],path[len(path)-1])
                    radj = get_adj(path[len(path)-2], path[0])
                true_adjacencies[duplicate][(genome,j)] = [ladj, radj]
    print(duplicate, true_adjacencies[duplicate])
'''

def get_list(cluster_path, cluster):
    files = []
    with open(f"{cluster_path}/{cluster}.txt") as f:
        for line in f:
            files.append(line.strip())
    return [os.path.basename(el).replace('.fna','') for el in files]

def unimog_to_ilp(unimog, lp, genome1, genome2):
    try:
        subprocess.run(f"dingII generate {unimog} -mm --writeilp {lp} -p {genome1} {genome2}", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    
def dcj_dist(unimog, matching, solution, genome1, genome2):
    try:
        dcj_out = subprocess.run(f"dingII parsesol {unimog} --solgur {solution} --matching {matching} -p {genome1} {genome2}", shell=True, check=True, capture_output=True, text=True).stdout
        dist = int(dcj_out.strip().split(" ")[2])
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    return dist

def ilp_gurobi(lp, solution):
    try:
        subprocess.run(f"gurobi_cl ResultFile={solution} Threads=1 {lp} >/dev/null", shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode())
        print(e)
        raise e
    
def pangraph_dcj(cluster, genomes):
    with open(f"{cluster}.dist", "w") as f:
        f.write("plasmid_1\tplasmid_2\tdist\n")
        for genome_1, genome_2 in tqdm(itertools.combinations(genomes,2)):
            unimog = f"{cluster}.unimog"
            lp = f"{genome_1}~{genome_2}.lp"
            solution = f"{genome_1}~{genome_2}.sol"
            matching = f"matching/{genome_1}~{genome_2}.unimog"
            unimog_to_ilp(unimog,lp,genome_1,genome_2)
            ilp_gurobi(lp, solution)
            dist = dcj_dist(unimog, matching, solution, genome_1, genome_2)
            f.write(f"{genome_1}\t{genome_2}\t{dist}\n")
    subprocess.run("rm *.sol *.lp",shell=True)

cluster_path = "/home/daria/Documents/projects/ABC/clades/lists"
cluster = "st69_cl375_community_0_subcommunity_499"
genomes = get_list(cluster_path, cluster)

pangraph_dcj(cluster,genomes)

G = nx.Graph()

graph = pp.Pangraph.from_json(f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json")
path_dict = graph.to_path_dictionary()

block_stats = graph.to_blockstats_df()
bl_count = graph.to_blockcount_df()
duplicates = list(block_stats[(block_stats["duplicated"]==True) & (block_stats["len"]>200)].index)
print(duplicates)

'''
for i in genomes:
    G.add_edges_from([(f"{i}_{j}",f"{i}_{j+1}") for j in range(len(path_dict[i])-1)])
    G.add_edge(f"{i}_{len(path_dict[i])-1}",f"{i}_0")
'''

for genome_1, genome_2 in itertools.combinations(genomes,2):
    with open(f"matching/{genome_1}~{genome_2}.unimog", "r") as f:
        lines = f.readlines()
        unimog_1 = lines[1].strip(")\n").split(" ")
        unimog_1.remove("")
        unimog_2 = lines[4].strip(")\n").split(" ")
        for i in range(len(unimog_1)):
            block_id = unimog_1[i].replace("+","").replace("-","").split("_")[0]
            if int(block_id) in duplicates:
                position=None
                G.add_node(f"{genome_1}_{i}")
                G.nodes[f"{genome_1}_{i}"]["duplicate"] = block_id
                try:
                    position = unimog_2.index(unimog_1[i])
                except:
                    try:
                        if unimog_1[0]=="+":
                            position = unimog_2.index(unimog_1[i.replace("+","-")])
                        else:
                            position = unimog_2.index(unimog_1[i.replace("-","+")])
                    except:
                        pass
                if position:
                    G.add_edge(f"{genome_1}_{i}", f"{genome_2}_{position}")
                    G.nodes[f"{genome_1}_{i}"]["duplicate"] = block_id
                    G.nodes[f"{genome_2}_{position}"]["duplicate"] = block_id



#nx.draw(G)
    
#plt.show()


json_dict = nx.cytoscape_data(G)
with open(f"{cluster}.json", "w") as file:
    json.dump(json_dict, file)

'''
for duplicate in duplicates:
    print(duplicate, block_stats.loc[duplicate, "n_strains"], block_stats.loc[duplicate, "count"])
'''
matchings = []
for c in nx.connected_components(G):
    for comm in nx.community.asyn_lpa_communities(G.subgraph(c)):
        genomes = [genome.split("_")[0]+"_"+genome.split("_")[1]+"_"+genome.split("_")[2] for genome in comm]
        conflict_genomes = [genome for genome in list(set(genomes)) if genomes.count(genome)>1]
        conflicts = {}
        final_comm = comm
        removed = []
        for genome in conflict_genomes:
            conflicts[genome] = [node for node in comm if node.split("_")[0]+"_"+node.split("_")[1]+"_"+node.split("_")[2]==genome]
            to_remove = min(conflicts[genome], key=lambda x: G.subgraph(comm).degree[x])
            final_comm.remove(to_remove)
            removed.append(to_remove)
        matchings.append(comm)
        for dup in removed:
            matchings.append({dup})

print(matchings)
                
