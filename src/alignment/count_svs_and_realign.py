# File at /scratch/groups/dpwall/personal/briannac/sars_cov2_microdeletions/realign_msa.py

import numpy as np
from collections import Counter
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio import Seq
import matplotlib
import numpy as np
import sys

ALIGNMENT = sys.argv[1]
REALIGNMENT=sys.argv[2]
SV_COUNTER = sys.argv[3]
SV_DICT = sys.argv[4]
REF_DICT = sys.argv[5]
REF = sys.argv[6]

insertions_counter = defaultdict(int)
deletions_counter = defaultdict(int)
insertions_dict = defaultdict(list)
deletions_dict = defaultdict(list)

all_fasta = [r for r in SeqIO.parse(ALIGNMENT,'fasta')]
ref = [r for r in all_fasta if REF in r.id][0].seq

################ SET UP DICTIONARY ##################
ref_positions = list(range(sum(np.array(ref)!='-')))
alignment_positions = [i for i,base in enumerate(np.array(ref)) if base!='-']
ref_dict = {}
prev_j = 0
j = 0
for i,base in enumerate(np.array(ref)):
    if base=='-':
        ref_dict[i] = prev_j
    else: 
        ref_dict[i] = j
        prev_j = j
        j = j + 1
np.save(REF_DICT, ref_dict)

##################### INSERTIONS #####################

insertions = np.array([i for i,base in enumerate(np.array(ref)) if base=='-'])
for f in all_fasta[1:]:
    seq = np.concatenate([['.'], np.array(f.seq), ['.']])
    seq = np.array(['-'] + ['-' for _ in ref] + ['-'])
    seq[insertions+1] = np.array(f.seq)[insertions]
    #seq = np.array(np.array(all_fasta[1].seq)[insertions])
    #if set(seq)!={'-'}: continue
    stops = np.where((seq[1:]=='-') & (seq[:-1]!='-'))[0]
    starts =  np.where((seq[1:]!='-') & (seq[:-1]=='-'))[0]
    for start, stop in zip(starts, stops):
        #if start==0: continue
        #if stop==len(f.seq):continue
        idx = (ref_dict[start+1], ''.join(seq[(start+1):(stop+1)]), start+1)
        insertions_counter[idx] = insertions_counter[idx] + 1
        insertions_dict[idx] = insertions_dict[idx] + [f.id]

ref_no_del = ''.join(np.array(ref)[np.array(ref)!='-'])
insertions_flanking = {}

for k in insertions_counter.keys():
    start = k[0]
    insert = k[1]
    flanks = str(ref_no_del[:start])+insert+str(ref_no_del[(start):])
    if flanks in insertions_flanking.keys():
        insertions_flanking[flanks] = insertions_flanking[flanks] + [k]
    else:
        insertions_flanking[flanks] = [k] 
        
synon_insertions = [d for d in insertions_flanking.values() if len(d)>1]

for d in synon_insertions:
    true_insertion = d[0]
    for false_insertion in d[1:]:
        insertions_dict[true_insertion]=list(set(insertions_dict[false_insertion]+insertions_dict[true_insertion]))
        insertions_counter[true_insertion]=len(insertions_dict[true_insertion])
        ids = [i for i,f in enumerate(all_fasta) if f.id in insertions_dict[false_insertion]]
        for i in ids:
            orig_seq = np.array(all_fasta[i].seq)
            orig_seq[false_insertion[2]:(false_insertion[2]+len(false_insertion[1]))] = ref[false_insertion[2]:(false_insertion[2]+len(false_insertion[1]))]
            orig_seq[true_insertion[2]:(true_insertion[2]+len(true_insertion[1]))] = '-'
            all_fasta[i].seq = Seq.Seq(''.join(orig_seq))
        insertions_counter.pop(false_insertion)
        insertions_dict.pop(false_insertion)

new_insertions_dict = dict()
new_insertions_counter = dict()

for d in insertions_dict.keys():
    new_insertions_dict[(d[0]+1, len(d[1]), 'I')] = insertions_dict[d]
    new_insertions_counter[(d[0]+1, len(d[1]), 'I')] = insertions_counter[d]

##################### DELETIONS #####################

for f in all_fasta[1:]:
    seq = np.array([s for i,s in enumerate(f.seq) if i not in insertions])
    if '-' not in seq: continue
    seq = np.concatenate([['.'], seq, ['.']])
    starts = np.where((seq[1:]=='-') & (seq[:-1]!='-'))[0]
    stops =  np.where((seq[1:]!='-') & (seq[:-1]=='-'))[0]
    for start, stop in zip(starts, stops):
        if start==0: continue
        if stop==len(f.seq):continue
        idx = (start, stop-start)
        deletions_counter[idx] = deletions_counter[idx] + 1
        deletions_dict[idx] = deletions_dict[idx]+ [f.id]

deletions_flanking = {}
for k in deletions_counter.keys():
    start = k[0]
    stop = k[0] + k[1]
    flanks = (str(all_fasta[0].seq[:start])+str(all_fasta[0].seq[(stop):]), k[1])
    if flanks in deletions_flanking.keys():
        deletions_flanking[flanks] = deletions_flanking[flanks] + [k]
    else:
        deletions_flanking[flanks] = [k] 
        
synon_deletions = [d for d in deletions_flanking.values() if len(d)>1]
for d in synon_deletions:
    true_deletion = d[0]
    for false_deletion in d[1:]:
        deletions_dict[true_deletion]=list(set(deletions_dict[false_deletion]+deletions_dict[true_deletion]))
        deletions_counter[true_deletion]=len(deletions_dict[true_deletion])
        ids = [i for i,f in enumerate(all_fasta) if f.id in deletions_dict[false_deletion]]
        for i in ids:
            orig_seq = np.array(all_fasta[i].seq)
            orig_seq[false_deletion[0]:(false_deletion[0]+false_deletion[1])] = ref[false_deletion[0]:(false_deletion[0]+false_deletion[1])]
            orig_seq[true_deletion[0]:(true_deletion[0]+true_deletion[1])] = '-'
            all_fasta[i].seq = Seq.Seq(''.join(orig_seq))
        deletions_counter.pop(false_deletion)
        deletions_dict.pop(false_deletion)

new_deletions_dict = dict()
new_deletions_counter = dict()

for d in deletions_dict.keys():
    new_deletions_dict[(d[0]+1, d[1], 'D')] = deletions_dict[d]
    new_deletions_counter[(d[0]+1, d[1], 'D')] = deletions_counter[d]
    

    
# Combine and save realigned files.
new_deletions_counter.update(new_insertions_counter)
sv_counter = new_deletions_counter
new_deletions_dict.update(new_insertions_dict)
sv_dict = new_deletions_dict

SeqIO.write(all_fasta, REALIGNMENT,'fasta')
np.save(SV_COUNTER, sv_counter)
np.save(SV_DICT, sv_dict)