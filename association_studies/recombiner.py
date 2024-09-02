import random
import argparse




parser = argparse.ArgumentParser(description='Simulate simple populations. Recombination + SNPs only, no selection. ')
parser.add_argument('-l', help = 'Length of the sequence [1000]', default = 1000, type =int)
parser.add_argument('-m', help = 'Mutation rate [0.01]', default = 0.01, type =float)
parser.add_argument('-n', help = 'Initial population size [10]', default = 10, type = int)
parser.add_argument('-g', help = 'Number of generations [100]', default = 10, type = int)
parser.add_argument('-x', help = 'Number of offspring (per subpopulation) [100]', default = 100, type = int)
parser.add_argument('-b', help = 'Bottleneck probability [0.1]', default = 0.1, type = float)
parser.add_argument('-f', help = 'Genomes per bottleneck event [5]', default = 5, type = int)
parser.add_argument('-r', help = 'Random seed [random]', type = int)

args = parser.parse_args()

if args.r: random.seed(args.r)

length = args.l
nucleotides = ['A', 'C', 'T', 'G']


def random_seq(nucleotides, n):
	random_set = [random.choice(nucleotides) for i in range(n)]
	return ''.join(random_set)



seq = random_seq(nucleotides, length)


out_table_indexes = []
def snps(seq, prob=0.001):
	outseq = ""
	for idx, i in enumerate(seq):
		r = random.uniform(0,1)
		if r <= prob:
			outseq += random.choice(["A","C","T","G"])
			if idx not in out_table_indexes:
				out_table_indexes.append(idx)
		else:
			outseq += i
	return outseq


num_seqs = args.n

init_pop = []

for n in range(num_seqs):
	outseq = snps(seq, prob = args.m)
	outclass = str(n)*length
	init_pop.append([outseq,outclass])

#print(init_pop)


########## phenotype

causal_allele_pos = int(random.uniform(0, len(init_pop[0][0])))
out_table_indexes.append(causal_allele_pos)
causal_allele = "A"

fixed = [c for c in init_pop[0][0]]
fixed[causal_allele_pos] = causal_allele
fixed = "".join(fixed)
init_pop[0][0] = fixed
for n in range(len(init_pop)-1):
	fixed = [c for c in init_pop[n+1][0]]
	fixed[causal_allele_pos] = "C"
	fixed = "".join(fixed)
	init_pop[n+1][0] = fixed

parent1 = init_pop[0][0]
parent2 = init_pop[1][0]


def recombine(seq, seq2):	
	xo_spot = random.uniform(1, len(seq[0]))
	xo_spot = int(xo_spot)
	#print(xo_spot)

	outseq = seq[0][:xo_spot]+seq2[0][xo_spot:]
	outseq = snps(outseq)
	
	outclass = seq[1][:xo_spot]+seq2[1][xo_spot:]

	return([outseq,outclass])




def random_cross(pop, xes = args.x):   ## init pop, number of crosses
	outpop = []

	for x in range(xes):
		a = random.choice(pop)
		b = random.choice(pop)
		while a==b:
			b = random.choice(pop)

		outpop.append(recombine(a,b))

	return outpop


pop2 = [random_cross(init_pop)]
# print(pop2[0])

# for i,n in enumerate(pop2[0]):
# 	print(f"{i}. {n[1]}")



ngen = args.g

bubble_probability = args.b
no_genomes_per_bubble = args.f



for gen in range(ngen-1):  
	for bubble,genomes in enumerate(pop2):
		roll = random.uniform(0,1)
		if roll < bubble_probability:
			new_bubble = random.sample(pop2[bubble], no_genomes_per_bubble)  
			pop2.append(random_cross(new_bubble))		### should be a tree
		pop2[bubble] = random_cross(pop2[bubble])




penetrance = 0.9


import os
import shutil
import numpy as np

outdir = os.getcwd()+"/out"
if not os.path.exists(outdir):
	os.mkdir("out")
else:
	shutil.rmtree(outdir)
	os.mkdir("out")


print(len(out_table_indexes))
out_table_indexes = sorted(out_table_indexes)

snp_table = [] ### przepisz to w pandas..

with open("out/phen_cov", "w") as phen_cov:
	for n_bubble, bubble in enumerate(pop2):
		group_mean = random.uniform(1,10)
		for i,n in enumerate(bubble):
			second_chrom = random.choice(pop2[n_bubble])
			phen = 0

			rng = np.random.default_rng()
			
			
			group_effect = rng.normal(group_mean, group_mean*0.1)

			error = rng.normal(5,2)
			effect = 0
			flag = 0


			if n[0][causal_allele_pos] == causal_allele:
				flag += 1
			if second_chrom[0][causal_allele_pos] == causal_allele:
				flag += 1

			roll = random.uniform(0,1)
			if roll < penetrance:
				if flag != 0:
					effect += rng.normal(30,5)
				if flag ==2:   									### Additive
					effect += rng.normal(30,5)

			phen = error + effect + group_effect ## mozna dodac interakcje 
			

				
			print(f"{n_bubble}_{i}\t{n[0]}\t{n[1]}\t{n_bubble}\t{phen}\t{flag}")  ## with seqs
			#print(f"{n_bubble}_{i}\t{n_bubble}\t{phen}\t{flag}")
			phen_cov.write(f"{n_bubble}_{i}\t{n_bubble}\t{phen}\t{flag}\n")
			with open(f"out/{n_bubble}_{i}.fa", "w") as fasta_out:
				fasta_out.write(f">{n_bubble}_{i}_1\n{n[0]}\n>{n_bubble}_{i}_2\n{second_chrom[0]}\n")

			snp_vect = []
			for snp in out_table_indexes:
				genotype = "".join(sorted([n[0][snp],second_chrom[0][snp]])) ### posortowane, zeby stracic phase ####### dodac bledy
				snp_vect.append(genotype)
			snp_table.append([f"{n_bubble}_{i}"]+snp_vect)

with open("out/snp_table", "w") as snp_out:
	snp_out.write("\t".join([str(k) for k in out_table_indexes])+"\n")
	for snp in snp_table:
		snp_out.write("\t".join(snp)+"\n")


with open("out/causal_allele", "w") as ca:
	ca.write(f"{causal_allele_pos}\t{causal_allele}\n")

with open("out/ref_seq.fa" , "w") as fq:
	fq.write(f">ref_seq\n{seq}\n")


with open("out/parent1.fa" , "w") as fq:
	fq.write(f">parent1\n{parent1}\n")

with open("out/parent2.fa" , "w") as fq:
	fq.write(f">parent2\n{parent2}\n")

### Add admixture model

### add snp table output