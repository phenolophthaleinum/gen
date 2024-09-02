import os
import threading
from timeit import default_timer as timer
from Bio import SeqIO
from io import StringIO

ret = []


def read_parallel(file, chunk):
    ret.append(file.read(chunk))


start = timer()
# X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta
# C:/Users/Maciej/vscode-projects/v-tests/seq.fasta
f_size = os.path.getsize("X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta")
part_num = 64
parts = [f_size // part_num + (1 if x < f_size % part_num else 0) for x in range(part_num)]
# print(parts)
threads = list()
with open('X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta', 'r') as f:
    for part in parts:
        # t = threading.Thread(target=f.read, args=(part,))
        t = threading.Thread(target=read_parallel, args=(f, part))
        threads.append(t)
        t.start()
for index, thread in enumerate(threads):
    thread.join()
    # print(f"thread {index} done")
txt = ''.join(ret)
# with StringIO(txt) as fasta_io:
#     records = SeqIO.parse(fasta_io, "fasta")
#     print(len(list(records)))
# recs = txt.split('>')[1::]
recs = ['>' + x for x in txt.split('>')][1::]
biorecs = [SeqIO.read(StringIO(e), 'fasta') for e in recs]
end = timer()
runtime = end - start
print(f"Done in {runtime:.6f}")
print(len(ret))
print(txt[0:20])
print(len(biorecs))
start = timer()
with open('X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta', 'r') as f2:
    txt = f2.read()
end = timer()
runtime = end - start
print(f"Done in {runtime:.6f}")
# print(ret)
# print(ret[0][-5:])
