import glob
import os

cluster_path = "/home/daria/Documents/projects/ABC/clades/lists" #path to directory with a text file for each cluster, which contains a list of fasta paths to memebers of the cluster
out = "pangraph_2"
clusters = [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]
clusters.remove("st69_cl461_community_0_subcommunity_499")
clusters.remove("st95_cl461_community_0_subcommunity_502")


def get_list(cluster):
    files = []
    with open(f"{cluster_path}/{cluster}.txt") as f:
        for line in f:
            files.append(line.strip())
    return files

def get_gff_fasta(cluster, gff_id):
    gene_fasta = f"unmapped_genes/{cluster}/{gff_id}.fna"
    with open(f"{ggcaller_path}/{cluster}/gene_calls.ffn") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == gff_id:
                with open(gene_fasta, "w") as output:
                    new_record = SeqRecord(record.seq, id=gff_id)
                    SeqIO.write(new_record, output, "fasta")
    return gene_fasta


rule all:
    input:
        [f"{out}/{cluster}/pangraph.gfa" for cluster in clusters]

rule ggcaller:
    input:
        fasta_list = lambda wildcards: f"{cluster_path}/{wildcards.cluster}.txt"
    output:
        ann_dir = directory(f"ggcaller/{{cluster}}")
    conda:
        "ggc_env"
    resources:
        mem_mb=lambda wildcards, attempt: 40000*attempt
    threads: 8
    shell:
        "ggcaller --refs {input.fasta_list} --out {output.ann_dir} --repeat --threads {threads}"


rule start:
    input:
        ann_dir = f"ggcaller/{{cluster}}",
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        same_start = directory("fastas/{cluster}")
    params:
        cluster = lambda wildcards: wildcards.cluster
    script: "change_start.py"


rule pangraph:
    input:
        "fastas/{cluster}"
    output:
        gfa = f"{out}/{{cluster}}/pangraph.gfa",
	    json = f"{out}/{{cluster}}/pangraph.json"
    params:
        fastas = lambda wildcards: [file.replace(".fna", "_start_moved.fna").replace("/home/daria/Documents/projects/ABC/goc/fastas", f"fastas/{wildcards.cluster}") for file in get_list(wildcards.cluster)]
    #conda:
    #    "mmseqs2"
    resources:
        mem_mb=lambda wildcards, attempt: 40000*attempt
    threads: 8
    shell:
        """
        pangraph build --circular -k minimap2 -s 20 -b 5 -a 150 --len 200 {params.fastas} > {output.json}
        pangraph export gfa --output {output.gfa} --minimum-length 200 {output.json}
        """
