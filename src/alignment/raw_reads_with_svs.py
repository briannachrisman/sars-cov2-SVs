import numpy as np
from Bio import SeqIO
import pandas as pd
import sys

DELETIONS_DICT = sys.argv[1]
DELETIONS_TABLE = sys.argv[2]
NCBI_READS_TABLE = sys.argv[3]
RAW_READS_TABLE = sys.argv[4]

deletions_dict = np.load(DELETIONS_DICT, allow_pickle=True).item()
ncbi_reads_table = pd.read_csv(NCBI_READS_TABLE, sep=',')

ids = {((i+'|').split('|')[1]) for i in np.concatenate([d for d in deletions_dict.values()])}
ids = {i.replace('hCoV-19/', '').replace('/2020', '') for i in ids}
raw_reads_names = ncbi_reads_table['Sample Name'].values
ids_dict = {}
raw_read_names = ncbi_reads_table['Sample Name'].values
for i in ids:
    ids_dict[i] = i.replace('hCoV-19/', '').replace('/2020', '').replace('-', '').replace('_', '')
    ids_dict[i] = ids_dict[i].replace('/', '') 
    ids_dict[i] = ids_dict[i].replace('AustraliaNSW', 'AustraliaVIC')
    ids_dict[i] = ids_dict[i].replace('AustraliaQLDID', 'AustraliaVIC')
    if 'CAMB' in ids_dict[i]:
        ids_dict[i] = ids_dict[i].replace('England', 'COGUK')
    if 'Denmark' in ids_dict[i] :
        if 'HH' in ids_dict[i] :
            ids_dict[i] = ids_dict[i].replace('Denmark', 'DK')
        if 'SSI' in ids_dict[i] :
            ids_dict[i] = ids_dict[i].replace('Denmark', 'DK')
    ids_dict[i] = ids_dict[i].replace('Scotland', 'COGUK')
    ids_dict[i] = ids_dict[i].replace('ALAB', '')
    ids_dict[i] = ids_dict[i].replace('COGUK', 'England')
    ids_dict[i] = ids_dict[i].replace('USAVADCLS', 'VADCLS')
    ids_dict[i] = ids_dict[i].replace('USAWAUW', 'WAUW')
    ids_dict[i] = ids_dict[i].replace('USA', '')

seq_names = ncbi_reads_table['Sample Name']
raw_reads_dict = {}
for i, name in enumerate(seq_names):
    raw_reads_dict[i] = name.replace('hCoV-19/', '').replace('/2020', '').replace('-', '').replace('_', '')
    raw_reads_dict[i] = raw_reads_dict[i].replace('/', '') 
    if 'CAMB' in raw_reads_dict[i]:
        raw_reads_dict[i] = raw_reads_dict[i].replace('England', 'COGUK')
    if 'Denmark' in raw_reads_dict[i] :
        if 'HH' in raw_reads_dict[i] :
            raw_reads_dict[i] = raw_reads_dict[i].replace('Denmark', 'DK')
        if 'SSI' in raw_reads_dict in raw_reads_dict[i] :
            raw_reads_dict[i] = raw_reads_dict[i].replace('Denmark', 'DK')
    raw_reads_dict[i] = raw_reads_dict[i].replace('ALAB', '')
    raw_reads_dict[i] = raw_reads_dict[i].replace('Scotland', 'COGUK')

    #if 'ScotlandCVR' in raw_reads_dict[i]:
    #    raw_reads_dict[i] = raw_reads_dict[i].replace('Scotland', 'COGUK')
    raw_reads_dict[i] = raw_reads_dict[i].replace('COGUK', 'England')
    raw_reads_dict[i] = raw_reads_dict[i].replace('USAVADCLS', 'VADCLS')
    raw_reads_dict[i] = raw_reads_dict[i].replace('USAWAUW', 'WAUW')
    raw_reads_dict[i] = raw_reads_dict[i].replace('USA', '')


n_seqs = [len(d) for d in deletions_dict.values()]
locs = [d[0] for d in deletions_dict.keys()]
lengths = [d[1] for d in deletions_dict.keys()]
types = [d[2] for d in deletions_dict.keys()]
countries=[', '.join(list({idx.split('|')[1].split('/')[1] for idx in deletions_dict[k]})) for k in deletions_dict.keys()]
available_samples = []
possible_samples = []
for k in deletions_dict.keys():
    seq_names=[d.split('|')[1] for d in deletions_dict[k]]
    seq_names = [s.replace('hCoV-19/', '').replace('/2020', '') for s in seq_names]
    possible_samples = possible_samples + [', '.join(
        [s.replace('hCoV-19/', '').replace('/2020', '') for s in seq_names])]
    available_samples = available_samples + [', '.join(
        [s.replace('hCoV-19/', '').replace('/2020', '') for s in seq_names if ids_dict[s] in raw_reads_dict.values()])]
deletions_df = pd.DataFrame([locs, lengths, types, n_seqs, countries, available_samples, possible_samples]).transpose()
deletions_df.columns=['Start Pos', 'Length (bp)', 'Structural Variant', '# Sequences', 'Countries', 'Available Samples', 'Possible Samples']
deletions_df = deletions_df.sort_values(['Start Pos'])


t = []
for row in deletions_df.iloc:
    if len(row['Available Samples'])==0: continue
    for a in row['Available Samples'].split(', '):
        ncbi_row = ncbi_reads_table.iloc[[r==ids_dict[a] for r in raw_reads_dict.values()]].iloc[0]
        t = t + [[ncbi_row['Run'], ncbi_row['LibraryLayout'], ncbi_row['Assay Type'], ncbi_row['Platform'], 'NCBI', a.replace('hCoV-19/', '').replace('/2020', '').replace('/', '_'), row['Start Pos'], a]]
        
raw_reads_df = pd.DataFrame(t)
raw_reads_df.columns = ['Run', 'LibraryLayout', 'Assay Type', 'Platform', 'database', 'short_seq_name', 'pos', 'full_seq_name']

deletions_df.to_csv(DELETIONS_TABLE, sep='\t', index=None)
raw_reads_df.to_csv(RAW_READS_TABLE, sep='\t', index=None)