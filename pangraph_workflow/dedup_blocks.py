import pypangraph as pp
import itertools
import subprocess
import os
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
import glob
import json

def get_adj(neighbour, block):
    if (neighbour[1] and block[1]) or (not neighbour[1] and not block[1]):
        adj = (neighbour[0], 't')
    else:
        adj = (neighbour[0], 'h')
    return adj

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
            unimog = f"unimogs/{cluster}.unimog"
            lp = f"{genome_1}~{genome_2}.lp"
            solution = f"{genome_1}~{genome_2}.sol"
            matching = f"matching_output/matching/{genome_1}~{genome_2}.unimog"
            unimog_to_ilp(unimog,lp,genome_1,genome_2)
            ilp_gurobi(lp, solution)
            dist = dcj_dist(unimog, matching, solution, genome_1, genome_2)
            f.write(f"{genome_1}\t{genome_2}\t{dist}\n")
    subprocess.run("rm *.sol *.lp",shell=True)

def get_occurance(counts, block_id):
    if block_id in counts.keys():
        counts[block_id]+=1
    else:
        counts[block_id]=1

cluster_path = "/home/daria/Documents/projects/ABC/clades/all_lists"
clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]

for cluster in clusters:
    genomes = get_list(cluster_path, cluster)

    if not os.path.exists(f"{cluster}.dist"):
        pangraph_dcj(cluster,genomes)

    G = nx.Graph()

    graph = pp.Pangraph.from_json(f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json")
    path_dict = graph.to_path_dictionary()

    block_stats = graph.to_blockstats_df()
    bl_count = graph.to_blockcount_df()
    duplicates = list(block_stats[(block_stats["duplicated"]==True) & (block_stats["len"]>200)].index)
    print(duplicates)



    for genome_1, genome_2 in itertools.combinations(genomes,2):
        with open(f"matching_output/matching/{genome_1}~{genome_2}.unimog", "r") as f:
            lines = f.readlines()
            unimog_1 = lines[1].strip(")\n").split(" ")
            unimog_1.remove("")
            unimog_2 = lines[4].strip(")\n").split(" ")
            counts_1 = {}
            counts_2 = {}
            for i in range(len(unimog_1)):
                block_id = unimog_1[i].replace("+","").replace("-","").split("_")[0]
                if int(block_id) in duplicates:
                    position=None
                    get_occurance(counts_1,block_id)
                    G.add_node(f"{genome_1}_{block_id}_{counts_1[block_id]}")
                    G.nodes[f"{genome_1}_{block_id}_{counts_1[block_id]}"]["duplicate"] = block_id
                    try:
                        position = unimog_2.index(unimog_1[i])
                        get_occurance(counts_2,block_id)
                    except:
                        try:
                            if unimog_1[0]=="+":
                                position = unimog_2.index(unimog_1[i.replace("+","-")])
                                get_occurance(counts_2,block_id)
                            else:
                                position = unimog_2.index(unimog_1[i.replace("-","+")])
                                get_occurance(counts_2,block_id)
                        except:
                            pass
                    if position:
                        G.add_edge(f"{genome_1}_{block_id}_{counts_1[block_id]}", f"{genome_2}_{block_id}_{counts_2[block_id]}")

    #json_dict = nx.cytoscape_data(G)
    #with open(f"{cluster}.json", "w") as file:
    #    json.dump(json_dict, file)

    matchings = []
    removed= []
    j=0
    with open(f"matching_output/matching_comms/{cluster}_matching.tsv", "w") as f:
        f.write("genome\tmarker\tmatching\n")
        with open(f"matching_output/matching_comms/{cluster}_conflicts.tsv", "w") as g:
            g.write("genome\tmarkers\n")
            for c in nx.connected_components(G):
                for comm in nx.community.asyn_lpa_communities(G.subgraph(c)):
                    j = j+1
                    genomes = [genome.split("_")[0]+"_"+genome.split("_")[1]+"_"+genome.split("_")[2] for genome in comm]
                    conflict_genomes = [genome for genome in list(set(genomes)) if genomes.count(genome)>1]
                    conflicts = {}
                    removed = []
                    for genome in conflict_genomes:
                        conflicts[genome] = [node for node in comm if node.split("_")[0]+"_"+node.split("_")[1]+"_"+node.split("_")[2]==genome]
                        to_remove = min(conflicts[genome], key=lambda x: G.subgraph(comm).degree[x])
                        removed.append(to_remove)
                    #matchings.append(comm)

                    for node in comm:
                        split = node.split("_")
                        genome = split[0]+"_"+split[1]+"_"+split[2]
                        marker = split[3]+"_"+split[4]
                        line = f"{genome}\t{marker}\t{j}\n"
                        f.write(line)

                    
                    for genome in conflict_genomes:
                        g.write(f"{genome}\t{conflicts[genome]}\n")


    #print(matchings)
                    
    print(removed)
