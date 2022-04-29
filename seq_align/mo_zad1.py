import random
import sys


length = int(sys.argv[1])
chroms = 2


transition_matrix = {
	'A': {'A':0.999, 'B':0.001},
	'B': {'A':0.001, 'B':0.999}
}

emission_matrix = {
	'A': {
		'A':{'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25},
		'C':{'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25},
		'T':{'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25},
		'G':{'A':0.25, 'C':0.25, 'T':0.25, 'G':0.25},
		},
	'B': {
		'A':{'A':0.1, 'C':0.1, 'T':0.7, 'G':0.1},
		'C':{'A':0.4, 'C':0.1, 'T':0.4, 'G':0.1},
		'T':{'A':0.7, 'C':0.1, 'T':0.1, 'G':0.1},
		'G':{'A':0.4, 'C':0.1, 'T':0.4, 'G':0.1}
		}
}


def random_seq(transition_matrix, emission_matrix, n):
	state = random.choice(list(transition_matrix.keys()))  ## init. uniform dist.
	seq = [random.choice(list(emission_matrix[state].keys()))]
	hidden = [state]
	for i in range(n-1):
		prev_n = seq[i]
		nucleotides = emission_matrix[state][prev_n]
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


outseq = []
for chrom in range(0,chroms):

	seq_oh = random_seq(transition_matrix, emission_matrix,  length)
	seq = seq_oh[0]
	seq = seq.lower()

	seq_h = seq_oh[1]
	seq_h = seq_h.lower()
	outseq.append(seq)

out_name = sys.argv[2]

with open(out_name, "w") as outfasta:
	for i,x in enumerate(outseq):
		outfasta.write(f">{i}\n{x}\n")

