import gffutils
from pymummer import coords_file, alignment, nucmer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from pathlib import Path
import subprocess
import pandas as pd

ggcaller = snakemake.input.ann_dir
fastas = snakemake.input.fastas
cluster = snakemake.params.cluster
os.mkdir(snakemake.output.same_start)
gff_dir = f"{ggcaller}/GFF"
genomes = [os.path.splitext(os.path.basename(el))[0] for el in fastas]
determine_strand = lambda match: '+' if match.on_same_strand() else '-'
dir = f"{ggcaller}"
genes = pd.read_csv(f"{dir}/gene_presence_absence_roary.csv")
num_isolates = len(genomes)
core = genes[genes["No. isolates"]==num_isolates]
gene = core.loc[0,"Gene"]
db = {}
for genome in genomes:
    db = gffutils.create_db(f"{gff_dir}/{genome}.gff",  ":memory:")
    start = float('inf')
    for id in core[core["Gene"]==gene][genome].values[0].split(";"):
        number = id.split("_")[2]
        gff_id = f"{genome}_{number}"
        try:
            strand = db[gff_id].strand
            start = min(start, db[gff_id].start)
            if strand == '-' and start == db[gff_id].start:
                end = db[gff_id].end
        except gffutils.FeatureNotFoundError: #some genes map across contigs and then aren't in the gff file
            fasta = Path(f"/home/daria/Documents/projects/ABC/goc/fastas/{genome}.fna")
            if not fasta.is_file():
                try:
                    subprocess.run(f"scp daria@codon-slurm-login.ebi.ac.uk:/nfs/research/zi/daria/game_of_clones/fastas/{genome}.fna /home/daria/Documents/projects/ABC/goc/fastas", shell=True, check=True, capture_output=True)
                except subprocess.CalledProcessError as e:
                    print(e.stderr.decode())
                    print(e)
                    raise e
            output_dir = Path(f"unmapped_genes/{cluster}")
            output_dir.mkdir(parents=True, exist_ok=True)
            gene_fasta = get_gff_fasta(cluster, gff_id)
            results = f"unmapped_genes/{cluster}/{gff_id}.nucmer"
            runner = nucmer.Runner(fasta, gene_fasta, results)
            runner.run()
            matches = sorted([match for match in coords_file.reader(results)], key=lambda x: x.qry_start, reverse = True)
            if matches[0].qry_length - 50 <= matches[0].hit_length_qry+matches[-1].hit_length_qry <= matches[0].qry_length + 50:
                strand = determine_strand(matches[0])
                start = matches[0].ref_start
            elif len(matches)==1 and matches[0].qry_length - 50 <= matches[0].hit_length_qry <= matches[0].qry_length + 50:
                strand = determine_strand(matches[0])
                start = matches[0].ref_start
            else:
                raise Exception("Unaccounted for case. Investigate.")
    if strand == '+':
        seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{genome}.fna", "fasta")
        rev = SeqRecord(seq_in.seq[start:]+seq_in.seq[:start], genome, "", "")
        SeqIO.write(rev, f"fastas/{cluster}/{genome}_start_moved.fna", "fasta")
    elif strand == '-':
        seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{genome}.fna", "fasta")
        move_seq = seq_in.seq[end:]+seq_in.seq[:end]
        rev = SeqRecord(move_seq.reverse_complement(), genome, "", "")
        SeqIO.write(rev, f"fastas/{cluster}/{genome}_start_moved.fna", "fasta")
