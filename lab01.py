import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
nucleotides = ['A', 'T', 'C', 'G']
weights = [0.1, 0.3, 0.2, 0.4]
length = 1000

def random_seq(nucleotides, weights, n):
	random_set = random.choices(nucleotides, weights, k=n)
	return ''.join(random_set)

seq = random_seq(nucleotides, weights, length)
print(seq)
record = SeqRecord(
	id='lab01_randomseq_weight',
	seq=Seq(seq)
)
SeqIO.write(record, "lab01_randomseq_weight.fasta", "fasta")
