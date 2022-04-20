import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

transition_matrix = {
	'A': {'A':0.95, 'B':0.05},
	'B': {'A':0.1, 'B':0.9}
}

emission_matrix = {
	'A': {'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25},
	'B': {'A':0.1, 'C':0.4, 'T':0.1, 'G':0.4}
}
length = 1000
def random_seq(transition_matrix, emission_matrix, n):
	state = random.choice(list(transition_matrix.keys()))  ## init. uniform dist.
	seq = []
	hidden = [state]
	for i in range(n):
		nucleotides = emission_matrix[state]
		keys = list(nucleotides.keys())   
		values = list(nucleotides.values())
		nuc_out= random.choices(keys, values)
		seq.extend(nuc_out)
		state = random.choices(
				list(transition_matrix[state].keys()),
				list(transition_matrix[state].values())
			)[0]
		hidden.extend(state)
	return [''.join(seq), ''.join(hidden)]

seq = random_seq(transition_matrix, emission_matrix, length)
print(seq)
record = SeqRecord(
	id='lab01_hmm',
	seq=Seq(seq[0])
)
SeqIO.write(record, "lab01_hmm.fasta", "fasta")

