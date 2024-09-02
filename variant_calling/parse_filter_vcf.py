import vcf


table = {'A': 4,
         'C': 5,
         'G': 6,
         'T': 7}
vaf_tresh = 0.05


def vaf(t1_alt, t1_ref):
    return t1_alt / (t1_alt + t1_ref)


vcf_reader = vcf.Reader(open('somatic_short_pass_norm.vcf', 'r'))
filtered_recs2 = []
for record in vcf_reader:
    # check if SNV else treat as indel
    if not (hasattr(record.samples[0].data, 'TIR') and hasattr(record.samples[0].data, 'TAR')):
        alt = record.ALT
        ref = record.REF
        alt_counts = record.samples[1].data[table[str(alt[0])]][0]
        ref_counts = record.samples[0].data[table[ref[0]]][0]
        if vaf(alt_counts, ref_counts) > vaf_tresh:
            filtered_recs2.append(record)
        continue
    alt_counts = record.samples[1].data.TIR[0]
    ref_counts = record.samples[0].data.TAR[0]
    if vaf(alt_counts, ref_counts) > vaf_tresh:
        filtered_recs2.append(record)

vcf_writer = vcf.Writer(open("somatic_norm_filtered.vcf", 'w'), vcf_reader)
for record in filtered_recs2:
    vcf_writer.write_record(record)
