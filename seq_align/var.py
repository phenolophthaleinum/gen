import random
# import matplotlib.pyplot as plt
import numpy as np
import math



# print(np.random.poisson(1, 1000)+1)
# plt.hist(np.random.poisson(2, 1000)+1, bins=100)
# plt.show()

#rng = np.random.default_rng()
# sample = rng.negative_binomial(1.5,0.1,1000)
# plt.hist(sample)
# plt.show()



def snps(seq, prob=0.001):
	outseq = ""
	for i in seq:
		r = random.uniform(0,1)
		if r <= prob:
			outseq += random.choice(["A","C","T","G"])
		else:
			outseq += i
	return outseq



def small_indels(seq, prob=0.0001):
	rng = np.random.default_rng()
	outseq = ""
	i=0
	while i < len(seq): 
		nuc = seq[i]
		r = random.uniform(0,1)
		if r <= prob:
			v_type = random.choice(["i","d"])
			size = int(np.ceil(rng.gamma(0.5,50))) 
			if v_type == "i":
				outseq += nuc 
				for n in range(0,size):
					outseq += random.choice(["A","C","T","G"])
			else:
				for n in range(0,size):
					i += 1
		else:
			outseq += nuc
		i += 1
	return outseq



def large_indels(seq, prob=0.0001):
	rng = np.random.default_rng()
	outseq = ""
	i=0
	while i < len(seq): 
		nuc = seq[i]
		r = random.uniform(0,1)
		if r <= prob:
			v_type = random.choice(["i","d"])
			size = int(np.ceil(rng.gamma(2,300))) 
			if v_type == "i":
				outseq += nuc 
				for n in range(0,size):
					outseq += random.choice(["A","C","T","G"])
			else:
				for n in range(0,size):
					i += 1
		else:
			outseq += nuc
		i += 1
	return outseq


def inv(seq, prob=0.0001):
	rng = np.random.default_rng()
	outseq = ""
	i=0
	while i < len(seq): 
		nuc = seq[i]
		r = random.uniform(0,1)
		if r <= prob:
			size = int(np.ceil(rng.gamma(2,300))) 
			outseq += seq[i:i+size][::-1].upper()
			i += size	### can be modified to allow for overlapping events (but biologically makes no sense?)
		else:
			outseq += nuc
			i += 1
	return outseq


def tan_dup(seq, prob=0.0001):
	rng = np.random.default_rng()
	outseq = ""
	i=0
	n_dup = rng.poisson(8)+1
	while i < len(seq): 
		nuc = seq[i]
		r = random.uniform(0,1)
		if r <= prob:
			size = int(np.ceil(rng.gamma(0.3,1000))) 
			outseq += seq[i:i+size].upper()
			for n in range(n_dup):
				outseq += seq[i:i+size].upper()
			i += size	### can be modified to allow for overlapping events (but biologically makes no sense?)
		else:
			outseq += nuc
			i += 1
	return outseq



def disp_dup(seq, prob=0.0001):
	rng = np.random.default_rng()
	i=0
	events = math.floor(len(seq)*prob)
	for event in range(events):
		outseq = ""
		seed = int(random.uniform(0,len(seq)))
		target = int(random.uniform(0,len(seq)))
		size = int(np.ceil(rng.gamma(3,200))) 
		outseq += seq[:target]
		outseq += seq[seed: seed+size].upper()
		outseq += seq[target:]
		seq = outseq
	return outseq



def nr_trans(seq, prob=0.0001):   ### non-reciprocal
	rng = np.random.default_rng() 
	i=0
	events = math.floor(len(seq)*prob)
	if events == 0 :
		return seq
	for event in range(events):
		outseq = ""
		size = int(np.ceil(rng.gamma(3,1000)))
		seed = int(random.uniform(0,len(seq)))
		target = int(random.uniform(0,len(seq)))
		if size > abs(seed-target):
			size = abs(seed-target)
		if size > len(seq)-seed:
			size = len(seq)-seed
		if target < seed:
			outseq += seq[:target]
			outseq += seq[seed: seed+size].upper()
			outseq += seq[target:seed]
			outseq += seq[seed+size:]
		elif seed < target:
			outseq += seq[:seed]
			outseq += seq[seed+size:target]
			outseq += seq[seed: seed+size].upper()
			outseq += seq[target:]
			
		else:
			outseq = seq
		seq = outseq
	return outseq


def r_trans(seq, prob=0.0001): ## reciprocal              
	rng = np.random.default_rng() 
	i=0
	events = math.floor(len(seq)*prob)
	if events == 0 :
		return seq
	for event in range(events):
		outseq = ""
		seed = int(random.uniform(0,len(seq)))
		target = int(random.uniform(0,len(seq)))
		s_size = int(np.ceil(rng.gamma(3,1000)))           ### can be constrained so they dont overlap
		t_size = int(np.ceil(rng.gamma(3,1000)))
		if target < seed:
			outseq += seq[:target]
			outseq += seq[seed: seed+s_size].upper()
			outseq += seq[target+t_size:seed]
			outseq += seq[target:target+t_size].upper()
			outseq += seq[seed+s_size:]
		elif seed < target:
			outseq += seq[:seed]
			outseq += seq[target:target+t_size].upper()
			outseq += seq[seed+s_size:target]
			outseq += seq[seed: seed+s_size].upper()
			outseq += seq[target+t_size:]
			
		else:
			outseq = seq
		seq = outseq
	return outseq