import sys
import math



in_fasta = sys.argv[1]
fadict = {}
with open(in_fasta) as plik:
	idx = ""
	for line in plik:
		if line.startswith(">"):
			idx = line.strip().split()[0][1:]
			fadict[idx] = ""
		else:
			fadict[idx] += line.strip()

chromname = list(fadict.keys())[0]

seq = fadict[chromname]
seq = seq.lower()
seq_l = len(seq)


def gc_content(seq):
	gc = 0
	for i in seq:
		if i in ["g","c"]:
			gc += 1
	return(gc/len(seq))


def gc_skew(seq):
	gc = {"g":0,"c":0}
	for i in seq:
		if i in gc.keys():
			gc[i] += 1
	g = gc["g"]
	c = gc["c"]
	return((c-g)/(g+c))


def shannon(seq):
	n_dic = {}
	seq_len = len(seq)
	for n in seq:
		if n in n_dic:
			n_dic[n] += 1
		else:
			n_dic[n] = 1
	se = 0
	for n in n_dic:
		pn = n_dic[n]/seq_len
		se += (pn)*(math.log(pn,2))
	return((-se))


def lzc(seq):
	seqs = []
	start = 0
	end = 1
	step = 1
	length = len(seq)
	while end <= length:
		if not seq[start:end] in seqs:
			seqs.append(seq[start:end])
			start += step
			end += 1
			step = 1

		else:
			end += 1
			step += 1
		
	return(len(seqs)/length)


def window(seq, size, step, func):
	length = len(seq)
	start = 0
	stop = size
	out = []
	while(stop <= length):        ### right border exception
		subseq = seq[start:stop]
		out.append([start, stop, func(subseq)])   # with step == size, or further aggregate
		#out.append(func(subseq))
		start += step
		stop += step
	return(out)


w_size = 500
w_step = w_size


def minmax(x, minx, maxx):
	return((x-minx)/(maxx-minx))



gc = window(seq,w_size,w_step,gc_content)
print(gc)

import pandas as pd
from matplotlib import pyplot as plt

gc_df = pd.DataFrame(gc)
# gc_df.hist()
# plt.show()


with open("gc_content.histo", "w") as out:
	for line in gc:
		val = line[2]
		if 0 < val <= 0.25:
			col = "103,201,129"
		elif 0.25 < val <= 0.35:
			col = "201,193,103"
		else:
			col = "204,116,92"

		lineout = [chromname]+line+[f"fill_color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")


shan = window(seq,w_size,w_step,shannon)



with open("shannon.histo", "w") as out:
	for line in shan:
		col = "red"
		lineout = [chromname]+line+[f"fill_color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")

lempel = window(seq,w_size,w_step,lzc)


with open("lempel.histo", "w") as out:
	for line in lempel:
		col = "blue"
		lineout = [chromname]+line+[f"color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")

print(len(seq))

skew = window(seq, w_size, w_step, gc_skew)
	


with open("gc_skew.histo", "w") as out:
	for line in skew:
		val = line[2]
		if val > 0:
			col = "red"
		else:
			col = "blue"
		lineout = [chromname]+line+[f"fill_color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")



with open("gc_skew_cumul.histo", "w") as out:
	gcs = 0
	for line in skew:
		gcs += line[2]
		col = "pink"
		lineout = [chromname]+[line[0]]+[line[1]]+[str(gcs)]+[f"color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")



# from plotnine import *
#
# skewdf = pd.DataFrame(skew, columns = ["start","stop","skew"])
# print(skewdf)
#
# skew_cuml = []
# cumul = 0
# for idx, row in skewdf.iterrows():
# 	print(row)
# 	cumul += row["skew"]
# 	skew_cuml.append(cumul)
# print(skew_cuml)
# skewdf["skew_cuml"] = skew_cuml
#
# skewplot = (ggplot(skewdf, aes(x="start", y=skew_cuml)) + geom_line())
# print(skewplot)




def kmer_freq(seq, k=20):
	length = len(seq)
	start = 0
	stop = k
	kmer_dic = {}
	while(stop <= length):        ### right border exception
		subseq = seq[start:stop]
		if subseq in kmer_dic:
			kmer_dic[subseq] += 1 
		else:
			kmer_dic[subseq] = 1
		start += 1
		stop += 1
	return kmer_dic

kmers = kmer_freq(seq)
kmers = sorted(kmers.items(), key= lambda n: n[1], reverse = True)

print(kmers[0])

max_kmer = kmers[0][0]

print(max_kmer)

def kmer_pos(seq, kmer):
	start = 0 
	stop = len(kmer)
	pos = []
	while(stop <= len(seq)):
		if seq[start:stop] == kmer:
			pos.append([start,stop])
		start += 1
		stop += 1
	return pos


kpos = kmer_pos(seq,max_kmer)

print(kpos)


with open("max_kmer.histo", "w") as out:
	col = "darkblue"
	for line in kpos:
		lineout = [chromname]+line+[f"color={col}"]
		out.write("\t".join([str(x) for x in lineout])+"\n")


def merge_overlaps(ranges): ## assume sorted
	out = []
	i = 0 
	jump = 1
	while i < len(ranges):
		start = ranges[i][0]
		stop = ranges[i][1]
		while i+jump < len(ranges):
			start_next = ranges[i+jump][0]
			stop_next = ranges[i+jump][1]
			if start_next <= stop:
				if stop_next > stop:
					stop = stop_next
				jump +=1 
			else:
				out.append([start,stop])
				jump +=1 
				break
		i += jump
		jump = 1
	out.append([start,stop])
	return(out)



kpos_merged = merge_overlaps(kpos)

print(kpos_merged)


with open("max_kmer_fix.links", "w") as out:
	col = "blue"
	used_ids = []
	for i in range(len(kpos_merged)):
		for j in range(len(kpos_merged)):
			if (i!=j) and int(f"{i}{j}") not in used_ids:
				# print([i, j], sep=' | ')
				lineout = [f"{i}{j}",chromname]+kpos_merged[i]+[f"color={col}"]              ## uwaga - tu jest ij i ji takie samo, czyli jest z powtorzeniami...
				out.write("\t".join([str(x) for x in lineout])+"\n")
				lineout = [f"{i}{j}",chromname]+kpos_merged[j]+[f"color={col}"]
				out.write("\t".join([str(x) for x in lineout])+"\n")
				used_ids.append(int(f"{i}{j}"))
		


