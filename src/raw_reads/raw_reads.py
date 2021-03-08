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

RAW_READS_TABLE = sys.argv[1] #'results/alignment/gisaid_hcov-19_test.raw_reads.tsv'
REFERENCE = sys.argv[2] #'data/reference_genomes/NC_045512.2.fasta'
FIG_PATH=sys.argv[3] # 'figs/raw_reads/gisaid_hcov-19_test'

def PlotReads(bam_file, chrom, start, stop, reference_file, pad=20, max_reads=20, title='raw reads', fig_file=None):
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        reads = [read for read in samfile.fetch(chrom, start=start, stop=STOP)]
    
    ref = [r for r in SeqIO.parse(reference_file, 'fasta')][0].seq

    # Set up data frame
    pos = list(range(start-pad, stop+pad+1))
    df = pd.DataFrame(np.zeros((min(len(reads)+1, max_reads+1),len(pos))))
    labels = pd.DataFrame([['' for _ in range(min(len(reads), max_reads)+1)] for _ in range(len(pos))]).transpose()
    labels.columns = [p+1 for p in pos]
    df.columns = [p+1 for p in pos]


    color_dict={'A':1, 'T':2, 'C':3, 'G':4, None:5, 'N':6}
    np.random.shuffle(reads)
    for i,read in (enumerate(reads)):
        prev_i=0
        query=read.seq
        if i>=(max_reads-1): break
        #[(q_i, ref_i) for (q_i, ref_i) in read.get_aligned_pairs()]
        for (q_i, r_i) in [(q_i, r_i) for (q_i, r_i) in read.get_aligned_pairs() if not ((r_i and (r_i>(stop+pad))) or (r_i and (r_i<(start-pad))))]:
            if type(r_i)!=int:
                if prev_i==0: continue
                c_name = str(prev_i) + '.'
                if not (c_name in df.columns):
                    idx=int(np.where(df.columns==prev_i)[0][0]+1)
                    df.insert(copy.copy(idx), c_name, 5.0)
                    labels.insert(idx, c_name, ['' for i in range(len(labels))])
                    prev_i = c_name 
                df.loc[i+1][c_name] = color_dict[query[q_i]]                
                labels.loc[i+1][c_name] = query[q_i]
                #continue
            else: 
                prev_i = r_i
                if type(q_i)!=int:
                    df.loc[i+1][r_i+1]=5
                    labels.loc[i+1][r_i+1] = '-'
                elif q_i>=len(query):
                    break
                    df.loc[i+1][r_i+1] = color_dict[query[q_i]]
                    labels.loc[i+1][r_i+1] = query[q_i]
                else:
                    df.loc[i+1][r_i+1] = color_dict[query[q_i]]
                    labels.loc[i+1][r_i+1] = query[q_i]
    labels.loc[0] = [ref[c-1] if type(c)==int else '-' for c in df.columns]
    labels.columns = [c if type(c)==int else ' ' for c in labels.columns]
    df.columns = [c if type(c)==int else ' ' for c in df.columns]
    df.loc[0] = [color_dict[ref[c-1]] if type(c)==int else 5 for c in df.columns]
    df['A']=1
    df['T']=2
    df['C']=3
    df['G']=4
    df['-']=5
    df['N']=6

    labels['A']=''
    labels['T']=''
    labels['C']=''
    labels['G']=''
    labels['-']=''
    labels['N']=''

    f, ax = plt.subplots(figsize=(20,max_reads/2))
    ax.set_aspect('equal')

    msa = sns.heatmap(df, annot=labels, fmt='', ax=ax, cmap=ListedColormap(['white', 'red', 'yellow', 'blue', 'green', 'white', 'grey']), cbar=False,
                      xticklabels=round(len(df.columns)/10))
    plt.xticks(rotation=0)
    ax.set_xlim([0,len(df.columns)-6.1])
    ax.axhline(1, *ax.get_xlim(), linewidth=5, color='k')
    ax.axhline(1, *ax.get_xlim(), linewidth=4, color='k')
    plt.yticks([])
    plt.title(title)
    plt.tight_layout()
    if fig_file:
        plt.savefig(fig_file)
    plt.show()

    
raw_reads = pd.read_csv(RAW_READS_TABLE, sep='\t')


for r in raw_reads.iloc:
    BAM_FILE = 'results/raw_reads/' + str(r['pos']) + '_' + r['short_seq_name'] + '.bam'
    print(BAM_FILE)
    START = r['pos']
    STOP = r['pos']+3
    if START>29500: continue
    if START<100: continue
    PAD = 50
    CHROM = 'NC_045512.2'
    MAX_READS=30
    PlotReads(BAM_FILE, CHROM, START, STOP, REFERENCE, pad=PAD, max_reads=MAX_READS,
              title=r['short_seq_name'] + ', variant at ' + str(r['pos']), fig_file=FIG_PATH + '/' + str(r['pos']) + '_' + r['short_seq_name'] + '.jpg')