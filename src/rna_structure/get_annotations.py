import pandas as pd
import sys

SV_FILE = sys.argv[1] #'results/alignment/gisaid_hcov-19_test.svs_table.tsv'
ANNOTATION_FILE = sys.argv[2] #'results/rna_structure/gisaid_hcov-19_test_annotations.txt'

svs = pd.read_csv(SV_FILE, sep='\t')
a = ''
for s in svs.iloc:
    if s['Structural Variant']=='D':
        if ((s['Start Pos']+s['Length (bp)']-1)>28900) or (s['Start Pos']<100):
            a = a + '%i %i 14 RED omark ' % (s['Start Pos'], s['Start Pos']+s['Length (bp)']-1)
        else:
            a = a + '%i %i 14 BLUE omark '% (s['Start Pos'], s['Start Pos']+s['Length (bp)']-1)
    else:
        if s['Start Pos']==29903:
            a = a + '%i %i 30 RED omark '% (s['Start Pos']-1, s['Start Pos'])      
        elif (s['Start Pos']>28900) or (s['Start Pos']<100):
            a = a #+ '%i %i 30 RED omark '% (s['Start Pos'], s['Start Pos']+1)
        else:
            a = a + '%i %i 30 GREEN omark '% (s['Start Pos'], s['Start Pos']+1)      
    a = a + '%i -2.0 -2.0 (%i) Label ' % (s['Start Pos'], s['Start Pos'])
a=a[:-1]
a = a.replace(' "', '"')
with open(ANNOTATION_FILE, 'w') as f:
    f.write(a)