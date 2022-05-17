#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import sys
import io
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-o", help="output", dest="output", type=str)

args = parser.parse_args()

def pairwise(iterable):
    """s -> (s0, s1), (s2, s3), (s4, s5), ..."""
    for a, b in zip(iterable, iterable[1:]):
        yield a, b 

WINDOW = 70_000

df = pd.read_csv(io.StringIO(sys.stdin.read()), skiprows=5, delimiter = "\t", usecols=[1, 2, 3], names=["position", "from", "to"])

ranges = list(range(0, df.position.max(), WINDOW))
counts = []
for i in range(len(ranges)-1):
    x = df.loc[(df.position >= ranges[i]) & (df.position <= ranges[i+1])]
    counts.append(len(x))

ranges = [f"{r[0]}-{r[1]}" for r in pairwise(ranges)]

fig, ax = plt.subplots(figsize=(21, 11))
ax.ticklabel_format(useOffset=False)
ax.plot(ranges, counts, c="purple")
ax.set_xticks(ranges)
plt.xticks(rotation=45)
plt.savefig(f"{args.output}.png", dpi=400)
