import pandas as pd
import numpy as np
import sys
import glob

DELETIONS_TABLE = sys.argv[1] #'results/alignment/gisaid_hcov-19_test.deletions_table.tsv'
BED_FILE = sys.argv[2] #'results/UCSC_genome_browser_annotations/deletions_test.bed'
dfs = []
for f in glob.glob('results/raw_reads/*.cov.vcf'):
    try:
        df = pd.read_csv(f, sep='\t', comment='#', header=None)
    except: continue
    dfs = dfs + [df]
vcfs = pd.concat(dfs)
vcfs.columns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL','FILTER', 'INFO', 'EXTRA1', 'EXTRA2']

deletions = pd.read_csv(DELETIONS_TABLE, sep='\t')
deletions = deletions[deletions['Start Pos']>50]
deletions = deletions[deletions['Start Pos']<29900]

bed = pd.DataFrame()
bed['chrom'] = ['NC_045512v2' for _ in deletions.iterrows()]
bed['start'] = deletions['Start Pos'].values
bed['end'] = deletions['Start Pos'].values + deletions['Length (bp)'].values
bed['name'] = [str(j) + str(i) + 'b' for i,j in zip(deletions['Length (bp)'], deletions['Structural Variant'])]
deletions.loc[[s > 29900 for s in deletions['Start Pos']],'Available Samples'] = np.nan
s = np.array([300 for _ in range(len(bed))])

for i in range(len(bed)):
    
    pos = deletions['Start Pos'].iloc[i]
    length = deletions['Length (bp)'].iloc[i]
    df_deletion = vcfs[(vcfs.POS==(pos-1)).values & np.array([len(r)==(length+1) for r in vcfs.REF])]
    depths = np.array([float(i.split(';')[0].replace('DP=', '')) for i in df_deletion['INFO']])
    freqs = np.array([float(i.split(';')[1].replace('AF=', '')) for i in df_deletion['INFO']])
    validated = sum((depths > 10) & (freqs > .6))>0
    
    if deletions['# Sequences'].iloc[i]>1:
        s[i] = 600
        pos = deletions['Start Pos'].iloc[i]
        length = deletions['Length (bp)'].iloc[i]
        
        
    if validated:
        s[i] = 1000

bed['strand'] = '+'
bed['score'] = s

bed.to_csv(BED_FILE, index=None, sep='\t', header=None)
with open(BED_FILE, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write('track name=Deletions description="Deletions" useScore=1 visibility=dense' + '\n' + content)