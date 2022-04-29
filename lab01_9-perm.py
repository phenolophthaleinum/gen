import io
import subprocess
from collections import defaultdict

import numpy as np
from Bio import SeqIO
import pandas as pd
import plotly.express as px
from PIL import Image
import matplotlib.pyplot as plt
import random


record, = SeqIO.parse("lab01_markov.fasta", 'fasta')
seq = str(record.seq)
# seq = "CCCCAACGTACAAGATCTACCATACCGCACGTCAGACACAGGTTCCCCTACGCCCTCACTTACCGATGTGGTAAATCTAATGAACTCTCCTAGCTGACGGTGCGAACAGTAACAGCCCTGCAAAGAACCTTCGGGTACTAAGAGCAGCAGAGACGACGTGCTGGACTACACATGATGCATTGCACCTTGGCCCG"
# seq = seq.lower()


def get_kmers(seq, size):
    kmers = defaultdict(int)
    n_kmers = len(seq) - size + 1
    for i in range(n_kmers):
        kmer = seq[i:i + size]
        kmers[kmer] += 1
    print(kmers)
    kmers_sum = sum(kmers.values())
    for kmer in kmers:
        kmers[kmer] = kmers[kmer] / kmers_sum
    return kmers


kmers = get_kmers(seq, 2)
max_kmer = sorted(kmers.items(), key=lambda n: n[1], reverse=True)[0]


def bootstrap(seq, reps, kmer):
    kmer_out = {kmer: []}
    for rep in range(1, reps):
        out = ""
        for i in (range(1, len(seq))):
            idx = random.randint(0, len(seq))
            out += seq[idx:idx + 1]
        kmer_out[kmer].append(get_kmers(out, 2)[kmer])
    return kmer_out


# bs = bootstrap(seq, 10000, max_kmer[0])


def perm(seq, reps, kmer):
    kmer_out = {kmer: []}
    print(seq)
    for rep in range(1, reps):
        out = "".join(random.sample(seq, len(seq)))
        kmer_out[kmer].append(get_kmers(out, 2)[kmer])
    return kmer_out


bs = perm(seq, 1000, max_kmer[0])

print(bs)

l_ci = np.quantile(bs[max_kmer[0]], 0.025)
l_hi = np.quantile(bs[max_kmer[0]], 0.975)

print((l_ci, l_hi))
print(get_kmers(seq, 2)[max_kmer[0]])


plt.hist(bs[max_kmer[0]])
plt.axvline(x=l_ci, color="black", linestyle=':')
plt.axvline(x=l_hi, color="black", linestyle=':')
plt.axvline(x=get_kmers(seq, 2)[max_kmer[0]], color="orange", linestyle='--')
plt.show()
