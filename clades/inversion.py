import pymummer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

'''
plasmid = "30224_1#353_3"
seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna", "fasta")
inv = SeqRecord(seq_in.seq[15000:31223], plasmid, "", "")
print(len(inv.seq))
SeqIO.write(inv, "inversion.fna", "fasta")

list_file = "/home/daria/Documents/projects/ABC/clades/lists/st131_cl416_community_0_subcommunity_503.txt"
'''
'''
plasmid = "31663_7#276_3"
seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna", "fasta")
inv = SeqRecord(seq_in.seq[24096:], plasmid, "", "")
print(len(inv.seq))
extracted_region = "insertion_5.fna"
SeqIO.write(inv, extracted_region, "fasta")

list_file = "/home/daria/Documents/projects/ABC/clades/lists/st73_cl703_community_0_subcommunity_493.txt"
'''

plasmid = "30859_5#362_3"
seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna", "fasta")
inv = SeqRecord(seq_in.seq[4981:6232], plasmid, "", "")
print(len(inv.seq))
extracted_region = "insertion_1.fna"
SeqIO.write(inv, extracted_region, "fasta")

list_file = "/home/daria/Documents/projects/ABC/clades/lists/st73_cl703_community_0_subcommunity_493.txt"

with open(list_file, "r") as f:
    files = f.read().split("\n")

orientation = {}
for file in files:
    plasmid = file.split("/")[-1].replace(".fna","")
    nucmer_file = f"{plasmid}_inv.nucmer"
    matches = []
    n = pymummer.nucmer.Runner(
        extracted_region,
        file,
        nucmer_file,
        min_id=80,
        breaklen=500,
        maxmatch=True,
        simplify=True
    )
    n.run()

    matches = [x for x in pymummer.coords_file.reader(nucmer_file)]
    for match in matches:
        if match.ref_length+1000 >= match.hit_length_ref >= match.ref_length-2000:
            if match.on_same_strand():
                orientation[plasmid] = '+'
            else:
                orientation[plasmid] = '-'

f = lambda item: files.index(f"/home/daria/Documents/projects/ABC/goc/fastas/{item[0]}.fna")
orientation = dict(sorted(orientation.items(), key=f))
with open("inversion_st73_cl703_c0_sc493.tsv", "w") as f:
    for key in orientation.keys():
        f.write(f"{key}\t{orientation[key]}\n")
