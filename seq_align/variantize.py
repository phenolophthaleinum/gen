import var

import os
import sys

in_fasta = sys.argv[1]
opts = sys.argv[2]

fadict = {}
with open(in_fasta) as plik:
	idx = ""
	for line in plik:
		if line.startswith(">"):
			fadict[line.strip()] = ""
			idx = line.strip()
		else:
			fadict[idx] += line.strip()
	


for key in fadict:
	seq = fadict[key].lower()
	for opt in opts:
		if opt == "s":
			seq = var.snps(seq)
		elif opt == "i":
			seq = var.small_indels(seq)
		elif opt == "l" or opt == "I":
			seq = var.large_indels(seq)
		elif opt == "v":
			seq = var.inv(seq)
		elif opt == "t":
			seq = var.tan_dup(seq)
		elif opt == "d":
			seq = var.disp_dup(seq)
		elif opt == "n":
			seq = var.nr_trans(seq)
		elif opt == "r":
			seq = var.r_trans(seq)
		else:
			print("WRONG OPTION")
			raise

	print(f"{key}\n{seq}\n")