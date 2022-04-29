import io
import subprocess
from collections import defaultdict

import numpy as np
from Bio import SeqIO
import pandas as pd
import plotly.express as px
from PIL import Image


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


filenames = [
    "lab01_hmm.fasta",
    "lab01_randomseq.fasta",
    "lab01_markov.fasta",
    "lab01_randomseq_weight.fasta",
]
cmd = 'jellyfish count -m 2 -s 100M -t 10 -C lab01_randomseq.fasta -o randomseq_2mer.jf'
plot_data = []
for file in filenames:
    # output = f"{file}_2mer_dump.fasta"
    # with open(output, 'w') as f:
    #     subprocess.run(['jellyfish', 'count', '-m', '2', '-s', '100M', '-t', '10', '-C', f'{file}', '-o', f'{file}_2mer.jf'])
    #     subprocess.run(['jellyfish', 'dump', f'{file}_2mer.jf', '-c'], stdout=f)
    #     subprocess.run(['rm', f'{file}_2mer.jf'])
    # #parsed_dict = SeqIO.to_dict(SeqIO.parse(output, 'fasta'))
    # t = pd.read_table(output, header=None, names=["key", 'val'], sep=" ")
    # kmers = t.set_index('key')['val'].to_dict()
    # kmers = {str(parsed_dict[e].seq): int(e) for e in parsed_dict}
    # >>>>>>>>>>>>>>>>> old
    record, = SeqIO.parse(file, 'fasta')
    s = str(record.seq)
    kmers = get_kmers(s, 3)
    compo = {base: s.count(base) for base in ['A', 'C', 'G', 'T']}
    probs = {}
    # nie wzor byl problemem
    for kmer in kmers:
        # probs.append(
        #     (kmers[kmer] / sum(list(kmers.values()))) / ((compo[kmer[0]] / len(s)) * (compo[kmer[1]] / len(s))))
        probs[kmer] = kmers[kmer] / ((compo[kmer[0]] / len(s)) * (compo[kmer[1]] / len(s)))
    # print(f"{file} probs:\n {probs}")
    # print(f"Prod: {np.prod(list(probs.values()))}")
    print(f"Mean: {np.mean(list(probs.values()))}")
    df = pd.DataFrame(kmers, index=[0])
    df_l = pd.melt(df, value_vars=list(df))
    df_l["file"] = file
    print(df_l)
    plot_data.append(df_l)
df_plot = pd.concat(plot_data)
print(plot_data)
    # na pewno histogram?
    # fig = px.histogram(df_l, x="variable")
fig = px.bar(df_plot, x='variable', y='value', facet_row="file", color="value",
             color_continuous_scale=px.colors.sequential.deep)
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
img_bytes = fig.to_image(format="png", width=800, height=800, scale=2)
image = Image.open(io.BytesIO(img_bytes))
image.save(f"3mer_hist.png")
