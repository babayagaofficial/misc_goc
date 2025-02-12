from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

plasmid = "30348_2#115_2"
seq_in = SeqIO.read(f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}.fna", "fasta")
rev = SeqRecord(seq_in.seq.reverse_complement(), plasmid, "", "")
SeqIO.write(rev, f"/home/daria/Documents/projects/ABC/goc/fastas/{plasmid}_rev.fna", "fasta")