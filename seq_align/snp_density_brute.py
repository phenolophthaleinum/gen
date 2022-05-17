import pandas as pd
from Bio import SeqIO
import plotly.express as px
from PIL import Image
import io


def get_snp(df, seqlen, size, step, gen_name):
    start = 0
    stop = size
    out = {}
    df["snp_pos"] = df["snp_pos"].apply(lambda x: int(x))
    while stop <= seqlen:
        snps = len(df.loc[(df['snp_pos'] >= start) & (df['snp_pos'] <= stop)])
        out[f"{start}-{stop}"] = snps
        start += step
        stop += step
    df_snps = pd.DataFrame(out, index=[0])
    df_snps_melt = df_snps.melt(value_name='snp_number', var_name="pos")
    df_snps_melt['alignment'] = f'ref: {gen_name}'
    return df_snps_melt

snpfiles = [
    # 'seq_align/nucmer_genus_align__filter_snp_TCI.snps',
    # 'seq_align/nucmer_strain_align3_filter_snp_TCI.snps'
    "seq_align/mch/mch_genus_align2.snps",
    'seq_align/mch/mch_strain_align2.snps'
]
coords = [
    'seq_align/mch/mch_genus_align2_new.filtered',
    'seq_align/mch/mch_strain_align2_new.filtered'
]
genomes = [
    'seq_align/licheniformis.fna',
    'seq_align/subtilis2.fna'
]
data = []
for snpfile, coords, genome in zip(snpfiles, coords, genomes):
    record, = SeqIO.parse(genome, 'fasta')
    seq_len = len(record.seq)
    df = pd.read_csv(coords, skiprows=4, delimiter="\t", usecols=[0, 1],
                     names=["start", "end"])
    l_start = df['start'].to_list()
    l_end = df['end'].to_list()
    alignments = list(zip(l_start, l_end))
    df_snps = pd.read_csv(snpfile, skiprows=5, delimiter="\t", usecols=[1],
                          names=["snp_pos"])
    df_snps['snp_pos'] = df_snps['snp_pos'].apply(lambda x: int(x))
    df_list = []
    for start, end in alignments:
        df_list.append(df_snps.loc[(df_snps['snp_pos'] >= start) & (df_snps['snp_pos'] <= end)])
    df_snps_filter = pd.concat(df_list)
    data.append(get_snp(df_snps_filter, seq_len, 75000, 50000, genome))
final_data = pd.concat(data)
print(final_data)
fig = px.line(final_data, x='pos', y='snp_number', facet_row="alignment", facet_row_spacing=0.2)
fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True, matches=None))
fig.update_xaxes(tickangle=45)
img_bytes = fig.to_image(format="png", width=1400, height=1200, scale=2)
image = Image.open(io.BytesIO(img_bytes))
image.save("snps-brute.png")

# df = pd.read_csv("seq_align/mch/mch_genus_align2_new.filtered", skiprows=4, delimiter = "\t", usecols=[0, 1], names=["start", "end"])
# l_start = df['start'].to_list()
# l_end = df['end'].to_list()
# alignments = list(zip(l_start, l_end))
# df_snps = pd.read_csv("seq_align/mch/mch_genus_align2.snps", skiprows=5, delimiter = "\t", usecols=[1], names=["snp_pos"])
# df_snps['snp_pos'] = df_snps['snp_pos'].apply(lambda x: int(x))
# df_list = []
# for start, end in alignments:
#     df_list.append(df_snps.loc[(df_snps['snp_pos'] >= start) & (df_snps['snp_pos'] <= end)])
# df_snps_filter = pd.concat(df_list)
#
#
# record, = SeqIO.parse('seq_align/licheniformis.fna', 'fasta')
# seq_len = len(record.seq)
# df_plot_strain = get_snp(df_snps_filter, seq_len, 75000, 50000)
# fig = px.line(df_plot_strain, x='pos', y='snp_number')
# fig.update_xaxes(tickangle=45)
# img_bytes = fig.to_image(format="png", width=1400, height=1200, scale=2)
# image = Image.open(io.BytesIO(img_bytes))
# image.save("snps-test_genus_brute_newfilter.png")
