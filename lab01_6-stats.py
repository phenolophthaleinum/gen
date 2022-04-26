import numpy as np
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import Bio.SeqUtils
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from PIL import Image
import io


def window(seq, size, step, func):
    length = len(seq)
    start = 0
    stop = size
    out = []
    while stop <= length:  ### right border exception
        subseq = seq[start:stop]
        out.append(func(subseq))
        start += step
        stop += step
    return out


def gc_content(seq):
    gc = 0
    for i in seq:
        if i in ["g", "c"]:
            gc += 1
    return gc / len(seq)


def gc_skew(seq):
    try:
        return (seq.count("G") - seq.count("C")) / (seq.count("G") + seq.count("C"))
    except ZeroDivisionError:
        return 1


def at_skew(seq):
    try:
        return (seq.count("A") - seq.count("T")) / (seq.count("A") + seq.count("T"))
    except ZeroDivisionError:
        return 1


def compo(seq):
    return {"A": seq.count("A"), "C": seq.count("C"), "T": seq.count("T"), "G": seq.count("G")}


def shannon(seq):
    composition = compo(seq)
    entropy = 0
    print(composition)
    for base in composition:
        p = composition[base] / float(len(seq))
        try:
            entropy_b = p * (np.math.log(p, 2))
        except ValueError:
            entropy_b = 0
        entropy += entropy_b
    return entropy * -1


filenames = [
    "lab01_hmm.fasta",
    "lab01_randomseq.fasta",
    "lab01_markov.fasta",
    "lab01_randomseq_weight.fasta",
]
records = []
data = []
data_gc = {}
data_comp = {}
data_gcskew = {}
data_atskew = {}
data_entropy = {}
for filename in filenames:
    with open(filename) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(record.id)
            s = str(record.seq)
            # print(f"Relative composition: {window(s, 100, 100, compo)}")
            # print(f"GC content: {window(s, 100, 100, Bio.SeqUtils.GC)}")
            # print(f"GC skew: {window(s, 100, 100, gc_skew)}")
            # print(f"AT skew: {window(s, 100, 100, at_skew)}")
            data_gc[filename] = window(s, 20, 5, Bio.SeqUtils.GC)
            data_gc["type"] = "gc_content"
            data.append(pd.DataFrame(data_gc))
            data_comp[filename] = window(s, 20, 5, compo)
            #data_comp["type"] = "composition"
            # data.append(data_comp)
            data_gcskew[filename] = window(s, 20, 5, gc_skew)
            data_gcskew["type"] = "gc_skew"
            data.append(pd.DataFrame(data_gcskew))
            data_atskew[filename] = window(s, 20, 5, at_skew)
            data_atskew["type"] = "at_skew"
            data.append(pd.DataFrame(data_atskew))
            data_entropy[filename] = window(s, 20, 5, shannon)
            data_entropy["type"] = "shannon_entropy"
            data.append(pd.DataFrame(data_entropy))
# print(gc_df_data)
#df = pd.DataFrame.from_dict(data, orient='columns')
df = pd.concat(data[-4:])

# df2 = pd.DataFrame.from_dict(data_comp, orient="columns")
df2 = pd.concat(pd.DataFrame(v).assign(file=k) for k, v in data_comp.items())
print(df2)
# print(gc_df)

fig = px.line(df, facet_col="type", facet_col_wrap=2)
for k in fig.layout:
    if re.search('yaxis[1-9]+', k):
        fig.layout[k].update(matches=None)
fig.update_yaxes(showticklabels=True)
fig2 = px.line(df2, facet_col="file", facet_col_wrap=2)
#fig.show()
#fig2.show()
img_bytes = fig.to_image(format="png", width=1200, height=800, scale=2)
img_bytes2 = fig2.to_image(format="png", width=1200, height=800, scale=2)
image = Image.open(io.BytesIO(img_bytes))
image2 = Image.open(io.BytesIO(img_bytes2))
image.save("stats1.png")
image2.save("stats-composit.png")
# fig.write_image("stats1.png")
# fig2.write_image("stats-composit.png")
# line_plot = df.plot.line()
# plt.show()
