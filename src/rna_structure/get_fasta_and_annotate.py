import pandas as pd
from Bio import SeqIO
from Bio import Seq
import sys

START = int(sys.argv[1])
STOP = int(sys.argv[2])
FASTA_OUT = sys.argv[3] #'results/rna_structure/300_500.fasta' 
REF = sys.argv[4] #'data/reference_genomes/NC_045512.2.fasta' 
SVS_TABLE = sys.argv[5] #'results/alignment/gisaid_hcov-19_2020_06_03_22.svs_table.tsv'
ANNOTATION_FILE = sys.argv[6] #'results/rna_structure/300_500.annotations.txt' 
CONSTRAINT_FILE = sys.argv[7]

svs = pd.read_csv(SVS_TABLE, sep='\t')
seq = [i for i in SeqIO.parse(REF, 'fasta')]
seq[0].seq = seq[0].seq[START:STOP]
SeqIO.write(seq, FASTA_OUT, 'fasta')
a = ''
for s in svs[(svs['Start Pos']>START) & (svs['Start Pos']<STOP)].iloc:
    if s['Structural Variant']=='D':
        a = a + '%i %i 20 GREEN omark '% (s['Start Pos']-0*START, s['Start Pos']+s['Length (bp)']-1-0*START)
    else:
        a = a + '%i %i 20 .99 .6 .2 omark '% (s['Start Pos']-0*START-1, s['Start Pos']-0*START)      
    a = a + '%i -2.0 -2.0 (%i) Label ' % (s['Start Pos']-0*START-1, s['Start Pos'])
a=a[:-1]
a = a.replace(' "', '"')
with open(ANNOTATION_FILE, 'w') as f:
    f.write(a)

with open(CONSTRAINT_FILE, 'w') as f:
    if START<20:
        f.write('P 1 0 3\n')
    else:
        f.write('P 1 0 3\n')
    f.write('P %i 0 5' % (len(seq[0].seq)-5))