import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import Seq
import sys

in_fasta = sys.argv[1]
out_fasta = sys.argv[2]
REF = sys.argv[3]
REF_DICT = sys.argv[4]

L_ref = 29903

all_fasta = [r for r in SeqIO.parse(in_fasta, 'fasta')]

ref = [r for r in all_fasta if REF in r.id][0].seq

################ SET UP DICTIONARY ##################
ref_positions = list(range(sum(np.array(ref)!='-')))
alignment_positions = [i for i,base in enumerate(np.array(ref)) if base!='-']
ref_dict = np.load(REF_DICT).item()

rev_ref_dict = {}
for k in ref_dict.keys():
    rev_ref_dict[ref_dict[k]] = k
    
rev_ref_dict_backwards = {}
for k in list(ref_dict.keys())[::-1]:
    rev_ref_dict_backwards[ref_dict[k]] = k

    
for ii in range(1,len(all_fasta)):
    s = list(all_fasta[ii].seq)
    
    # Mask sites as described in virological.org.
    for i in [187, 1059, 2094, 3037, 3130, 6990, 8022, 10323, 10741, 11074, 13408,
              14786, 19684, 20148, 21137, 24034, 24378, 25563, 26144, 26461, 26681,
              28077, 28826, 28854, 29700, 4050, 13402, 11083, 15324, 21575]:
         s[rev_ref_dict[i-1]] = 'N'
            
    # Mask beginning 55 sites.
    for i in range(rev_ref_dict[55]):
        s[i] = 'N'
    
    # Mask end 100 sites.
    for i in range(rev_ref_dict_backwards[L_ref-100],len(s)):
        s[i] = 'N'
        
    all_fasta[ii].seq = Seq.Seq(''.join(s))
all_fasta[0].id = REF
SeqIO.write(all_fasta, out_fasta, 'fasta')