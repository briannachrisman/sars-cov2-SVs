import matplotlib
matplotlib.use('Agg')
import json
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import sys

NODES_FILE = sys.argv[1] #'results/tree/gisaid_hcov-19_test.nodes.json'
DELETIONS_COUNTER = sys.argv[2] #'results/alignment/gisaid_hcov-19_test.deletions_counter.npy'
REF_DICT =  sys.argv[3]
HISTOGRAM_FILE = sys.argv[4]

with open(NODES_FILE) as json_file:
    data = json.load(json_file)

ref_dict = np.load(REF_DICT).item()

mutations = Counter(np.concatenate([np.array(data['nodes'][i]['muts']) for i in data['nodes']])).most_common()
svs = [(i,j) for i,j in mutations if '-' in i]


d_pos = [ref_dict[int(d[1:-1])] for d,i in svs if d[0] if d[-1]=='-']
d_n_times = [i for d,i in svs if d[-1]=='-']
i_pos = [ref_dict[int(d[1:-1])] for d,i in svs if d[0] if d[0]=='-']
i_n_times = [i for d,i in svs if d[0]=='-']


svs_counter = np.load(SVS_COUNTER).item()
print(len(svs_counter), 'unique SVs.')

f,ax =plt.subplots(figsize=(10,6), nrows=2, sharey='row', sharex='row')
ax[0].plot([d[0] for d in svs_counter.keys() if d[2]=='D'],
         [svs_counter[d] for d in svs_counter.keys() if d[2]=='D'], '.', color='orange')
ax[0].plot(d_pos, d_n_times, 'bx')
ax[0].hist(np.concatenate([np.array([d[0] for _ in range(svs_counter[d])]) for d in svs_counter.keys() if d[2]=='D']),
           bins=np.linspace(0,30000,100), alpha=.2, color='orange')
ax[0].set_yscale('log')
ax[0].set_xlim(0,30000)
ax[0].legend(['Sequences', 'Lineages'])
ax[0].set_title('Deletions')


ax[1].plot([d[0] for d in svs_counter.keys() if d[2]=='I'],
         [svs_counter[d] for d in svs_counter.keys() if d[2]=='I'], '.', color='orange')
ax[1].plot(i_pos, i_n_times, 'bx')
ax[1].hist(np.concatenate([np.array([d[0] for _ in range(svs_counter[d])]) for d in svs_counter.keys() if d[2]=='I']), 
           bins=np.linspace(0,30000,100), alpha=.3, color='orange')
ax[1].set_yscale('log')
ax[1].set_xlim(0,30000)
ax[1].set_title('Insertions')
plt.xlabel('Start Position')
plt.ylabel('# With Structural Variant')
plt.tight_layout()
plt.show()
f.savefig(HISTOGRAM_FILE, bbox_inches='tight')