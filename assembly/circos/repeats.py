import pandas as pd
from Bio import SeqIO
from collections import defaultdict

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


def get_repeats(df, seqlen, size, step, chrom):
    # chroms = df[0].unique()
    # df["P1"] = df["P1"].apply(lambda x: int(x))
    # tags = ' - '.join(df[["TAG_R", "TAG_Q"]].values[0])
    start = 0
    stop = size
    out = {}
    if size >= seqlen:
        repeats = len(df.loc[(df[1] >= start) & (df[2] <= stop + 1) & (df[0] == chrom)])
        out[f"{start}-{stop}"] = repeats
        # print(out)
        return out
    while stop <= seqlen:
        repeats = len(df.loc[(df[1] >= start) & (df[2] <= stop + 1) & (df[0] == chrom)])
        out[f"{start}-{stop}"] = repeats
        start += step
        stop += step
    # print(out)
    # df_snps = pd.DataFrame(out, index=[0])
    # df_snps_melt = df_snps.melt(value_name='snp_number', var_name="pos")
    # df_snps_melt["alignment"] = tags
    return out


df = pd.read_table("ragtag.scaffold.fasta.out.gff.repeats.tbl", header=None, delimiter=" ")
print(df)

res = {}
with open("ragtag.scaffold.fasta") as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        chrom_len = len(record.seq)
        print(chrom_len)
        chrom_name = record.id
        print(chrom_name)
        res[chrom_name] = get_repeats(df, chrom_len, 50000, 25000, chrom_name)
# print(res)
    # records = list(SeqIO.parse(handle, "fasta"))
# print(records)
# print(len(records))
# res = get_repeats(df, 7500000, 50000, 25000)
# print(res['NC_039455.1_RagTag'])
data = defaultdict(list)
# print(res["NC_040193.1_RagTag"].items())
for d in res:
    pairs = res[d].items()
    for p in pairs:
        pos = p[0]
        value = p[1]
        start = pos.split("-")[0]
        end = pos.split("-")[1]
        data["chr"].append(d)
        data["start"].append(start)
        data["end"].append(end)
        data['value'].append(value)

df_final = pd.DataFrame(data)
print(df_final)
df_final.to_csv("repeats.histo", header=False, sep='\t', index=False)
