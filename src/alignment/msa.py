import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import pysam 
import sys
import copy

DELETIONS_TABLE = sys.argv[1]  #'results/alignment/gisaid_hcov-19_test.svs_table.tsv'
MSA_FILE = sys.argv[2]  #'results/alignment/gisaid_hcov-19_test.realigned.fasta'
REFERENCE = sys.argv[3]  #'NC_045512.2'
REF_DICT = sys.argv[4]  #'results/alignment/gisaid_hcov-19_test.alignment_to_ref_pos_dict.npy'
MSA_FILEPATH = sys.argv[5]  # 'figs/msa/gisaid_hcov-19_test/'

def PlotMSA(alignment, start, stop, ref_seq, ref_dict, max_seqs=10, title='msa', fig_file=None):
    
    ref = copy.deepcopy(ref_seq)
    seqs = [ref] + alignment[:max_seqs]
    seq_names = [s.id.split('|')[1].replace('/2020', '').replace('hCoV-19/', '') for s in seqs]
    color_dict={'A':1, 'T':2, 'C':3, 'G':4, 'N':5, '-':6}
    t = pd.DataFrame(np.array([[np.array([color_dict[s] for s in a.seq[start:stop]] + [1,2,3,4,5,6])] for a in seqs])[:,0,:])
    labels = pd.DataFrame(np.array([[np.array([s for s in a.seq[start:stop]] +['','','','','',''])] for a in seqs])[:,0,:])
    t.index = seq_names
    t.columns = [ref_dict[a+1] for a in list(range(start, stop))] + ['A', 'T', 'C', 'G', 'N', '-']
    labels.index = seq_names
    labels.columns = [ref_dict[a+1] for a in list(range(start, stop))] + ['A', 'T', 'C', 'G', 'N', '-']
    f, ax = plt.subplots(figsize=(20,len(seqs)/2))
    ax.set_aspect('equal')
    msa = sns.heatmap(t, annot=labels, fmt='', ax=ax, cmap=ListedColormap(['red', 'yellow', 'blue', 'green', 'grey', 'white']), cbar=False,
                      xticklabels=round(len(t.columns)/10))
    plt.xticks(rotation=0)
    ax.set_xlim([0,len(t.columns)-6])
    ax.axhline(1, *ax.get_xlim(), linewidth=5, color='k')
    plt.title(title)
    plt.tight_layout()
    if fig_file:
        plt.savefig(fig_file)
    plt.show()
    
   
    
deletions = pd.read_csv(DELETIONS_TABLE, sep='\t')
alignment = [r for r in SeqIO.parse(MSA_FILE, 'fasta')]
alignment = [r for r in alignment if '|' in r.id]
ref_dict = np.load(REF_DICT).item()
ref_seq = [a for a in alignment if REFERENCE in a.id][0]
rev_ref_dict = dict()
for k in ref_dict.keys():
    rev_ref_dict[ref_dict[k]] = k

    
for d in deletions.iloc:
    pad = 20
    if d['Start Pos']<50: continue
    if d['Start Pos']>29500: continue
    START = rev_ref_dict[int(d['Start Pos'])]-pad
    STOP = rev_ref_dict[min(29903,int(d['Start Pos']) + int(d['Length (bp)']))]+pad
    samples = d['Possible Samples'].split(', ')
    alignment_filt = [a for a in alignment if a.id.split('|')[1].replace('hCoV-19/', '').replace('/2020', '') in samples]
    PlotMSA(alignment_filt, START, STOP, ref_seq, ref_dict, max_seqs=10, title='Multiple Sequence Alignment',
            fig_file=MSA_FILEPATH + '/%i_%s_%i_bp.jpg' % (d['Start Pos'], d['Structural Variant'], d['Length (bp)']))