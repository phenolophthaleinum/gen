import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

transition_matrix = {
'A': {'A':0.6, 'C':0.1, 'T':0.2, 'G':0.1},
'C': {'A':0.1, 'C':0.5, 'T':0.1, 'G':0.3},
'T': {'A':0.4, 'C':0.05, 'T':0.5, 'G':0.05},
'G': {'A':0.05, 'C':0.2, 'T':0.05, 'G':0.7}
}
length = 1000
def random_seq(transition_matrix, n):
	seq = [random.choice(list(transition_matrix.keys()))]  ## init. uniform dist.
	for i in range(n):
		key = seq[i]
		nucleotides = transition_matrix[key]
		keys = list(nucleotides.keys())   
		values = list(nucleotides.values())
		nuc_out= random.choices(keys, values)
		seq.extend(nuc_out)
	return ''.join(seq)

seq = random_seq(transition_matrix, length-1)

record = SeqRecord(
	id='lab01_markov',
	seq=Seq(seq)
)
SeqIO.write(record, "lab01_markov.fasta", "fasta")
