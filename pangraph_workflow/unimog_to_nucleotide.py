import pypangraph as pp
from Bio.Seq import Seq
import glob
import subprocess
import os

def read_in_unimog(unimog_path, cluster):
    genomes = {}
    with open(unimog_path,"r") as f:
        for line in f:
            if line[0]==">":
                entry = line.strip().replace(">","")
            else:
                genome = line.strip(")\n").split(" ")
                clean_genome = [el.split("_")[0] for el in genome if el!='']
                if entry[0]=='I':
                    name = f"{cluster}_{entry}"
                    if name in genomes.keys():
                        genomes[f"{name}_2"] = clean_genome
                    else:
                        genomes[name] = clean_genome
                else:
                    genomes[entry] = clean_genome
    return genomes

cluster_path = "/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/res_adj"
clusters = [os.path.basename(el).replace('_edists.txt','') for el in glob.glob(f"{cluster_path}/*_edists.txt")]

for cluster in clusters:

    unimog_path = f"/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/unimogs/{cluster}.unimog"

    try:
        subprocess.check_call(f"python /home/daria/Documents/projects/ABC/spp_dcj_experimental/scripts/parse_solution.py /home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/ilp/{cluster}_pei.txt /home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/full_adj/{cluster}.txt /home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/ilp/{cluster}_ids.txt /home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/sol/{cluster}.sol -t -fmb /home/daria/Documents/projects/ABC/pangraph_workflow/fitch/ranges/{cluster}_range.tsv --write-unimog {unimog_path}", shell=True)
    
        graph = pp.Pangraph.from_json(f"/home/daria/Documents/projects/ABC/pangraph_workflow/pangraph/{cluster}/pangraph.json")
        genomes = read_in_unimog(unimog_path, cluster)

        seqs = {}
        for entry in genomes.keys():
            seqs[entry] = Seq("")
            for el in genomes[entry]:
                if el[0] == '-':
                    id = int(el[1:])
                    block = graph.blocks[id]
                    seq = Seq(block.consensus())
                    seqs[entry] = seqs[entry] + seq.reverse_complement()
                else:
                    id = int(el)
                    block = graph.blocks[id]
                    seq = Seq(block.consensus())
                    seqs[entry] = seqs[entry] + seq

        os.mkdir(f"recon_fastas/{cluster}")
        for entry in seqs.keys():
            with open(f"recon_fastas/{cluster}/{entry}.fasta", "w") as fasta:
                fasta.write(f">{entry}\n{str(seqs[entry])}\n")
    except subprocess.CalledProcessError as e:
        pass


        
            