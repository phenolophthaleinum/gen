import subprocess

import numpy as np
from Bio import SeqIO
import pandas as pd

filenames = [
    "lab01_hmm.fasta",
    "lab01_randomseq.fasta",
    "lab01_markov.fasta",
    "lab01_randomseq_weight.fasta",
]
cmd = 'jellyfish count -m 2 -s 100M -t 10 -C lab01_randomseq.fasta -o randomseq_2mer.jf'

for file in filenames:
    output = f"{file}_2mer_dump.fasta"
    with open(output, 'w') as f:
        subprocess.run(['jellyfish', 'count', '-m', '2', '-s', '100M', '-t', '10', '-C', f'{file}', '-o', f'{file}_2mer.jf'])
        subprocess.run(['jellyfish', 'dump', f'{file}_2mer.jf', '-c'], stdout=f)
        subprocess.run(['rm', f'{file}_2mer.jf'])
    #parsed_dict = SeqIO.to_dict(SeqIO.parse(output, 'fasta'))
    t = pd.read_table(output, header=None, names=["key", 'val'], sep=" ")
    kmers = t.set_index('key')['val'].to_dict()
    # kmers = {str(parsed_dict[e].seq): int(e) for e in parsed_dict}
    record, = SeqIO.parse(file, 'fasta')
    s = record.seq
    compo = {base: s.count(base) for base in ['A', 'C', 'G', 'T']}
    probs = []
    # no i tu jest problem - jaki jest wzor w koncu xd
    for kmer in kmers:
        probs.append(
            (kmers[kmer] / sum(list(kmers.values()))) / ((compo[kmer[0]] / len(s)) * (compo[kmer[1]] / len(s))))
    print(f"{file} probs:\n {probs}")
    print(f"{np.prod(probs)}")

