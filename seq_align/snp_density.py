import pandas as pd
import os
import plotly.express as px
from PIL import Image
from Bio import SeqIO
import io
import re


def get_snp(df, seqlen, size, step):
    start = 0
    stop = size
    out = {}
    df["P1"] = df["P1"].apply(lambda x: int(x))
    tags = ' - '.join(df[["TAG_R", "TAG_Q"]].values[0])
    while stop <= seqlen:
        snps = len(df.loc[(df['P1'] >= start) & (df['P1'] <= stop)])
        out[f"{start}-{stop}"] = snps
        start += step
        stop += step
    df_snps = pd.DataFrame(out, index=[0])
    df_snps_melt = df_snps.melt(value_name='snp_number', var_name="pos")
    df_snps_melt["alignment"] = tags
    return df_snps_melt


snpfiles = [
    'seq_align/nucmer_genus_align__filter_snp_TCI.snps',
    'seq_align/nucmer_strain_align3_filter_snp_TCI.snps'
]
genomes = [
    'seq_align/licheniformis.fna',
    'seq_align/subtilis2.fna'
]
data = []
for snpfile, genome in zip(snpfiles, genomes):
    record, = SeqIO.parse(genome, 'fasta')
    filename = snpfile.split('.')[0]
    out = f'{filename}_fix.snps'
    # os.system(f'tail -n +3 {snpfile} > {out}')
    df = pd.read_table(out, header=None,
                       names=["P1", "SUB_R", "SUB_Q", "P2", "BUFF", "DIST", "FRM1", "FRM2", "TAG_R", "TAG_Q"])[1::]
    seq_len = len(record.seq)
# df_snps = get_snp(df, 75000, 25000)
# df_snps_melt = df_snps.melt(value_name='snp_number', var_name="pos")
    data.append(get_snp(df, seq_len, 75000, 25000))
final_data = pd.concat(data)
print(final_data)
fig = px.line(final_data, x='pos', y='snp_number', facet_row="alignment", facet_row_spacing=0.2)
# for k in fig.layout:
#     if re.search('xaxis[1-9]+', k):
#         fig.layout[k].update(matches=None)
# fig.update_yaxes(showticklabels=True)
# fig.update_yaxes(matches=None)
fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True, matches=None))
# for k in fig.layout:
#     if re.search('xaxis[1-9]+', k):
#         fig.layout[k].update(matches=None)
# fig.update_yaxes(showticklabels=True)
img_bytes = fig.to_image(format="png", width=1400, height=1200, scale=2)
image = Image.open(io.BytesIO(img_bytes))
image.save("snps-test.png")
# fig.show()
