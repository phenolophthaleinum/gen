from ctypes import *
import os


so_file = "C:/Users/Maciej/vscode-projects/v-tests/test_py.so"
if os.name == 'nt':
    so_file = "C:/Users/Maciej/vscode-projects/v-tests/test_py.dll"
v_lib = CDLL(so_file)
v_lib.v_parallel_read.argtypes = c_wchar_p, c_int
print(v_lib.__dict__)
print(v_lib.v_parallel_read.argtypes)
# assert len(v_lib.v_parallel_read(c_wchar_p("C:/Users/Maciej/vscode-projects/v-tests/seq.fasta"), c_int(4))) > 0, "Error reading file in parallel from V."
v_lib.v_parallel_read("C:/Users/Maciej/vscode-projects/v-tests/seq.fasta", 4)
