from Bio import Phylo
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import networkx as nx
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import sys

TREE_FILE = sys.argv[1] #'results/tree/gisaid_hcov-19_test.tree.nwk'
NODES_FILE = sys.argv[2] #'results/tree/gisaid_hcov-19_test.nodes.json'
DELETIONS_COUNTER = sys.argv[3] #'results/alignment/gisaid_hcov-19_test.deletions_counter.npy'
TREE_DIR = sys.argv[4] #'figs/tree/'

with open(NODES_FILE) as json_file:
    nodes = json.load(json_file)
REF = 'NC_045512.2'
tree = Phylo.read(TREE_FILE, "newick")
tree.root_with_outgroup(REF)
tree.prune(REF)
tree.root.name = REF
tree_nx = Phylo.to_networkx(tree)


nodelist = [n for n in tree_nx.nodes() if 'NODE' not in n.name]
ROOT = [n for n in tree_nx.nodes()][0].name
pos =  graphviz_layout(tree_nx, prog='/share/software/user/open/graphviz/2.40.1/bin/twopi', root=[n for n in tree_nx.nodes()][0])


for loc, length in [i[0] for i in Counter(np.load(DELETIONS_COUNTER).item()).most_common()]:
    node_color = ['blue' for n in nodelist]
    node_size = [1 for _ in nodelist]
    for i,n in enumerate(nodelist):
        if n.name:
            if ROOT == n.name:
                node_color[i] = 'blue'
                node_size[i] = 10
            elif 'NODE' in n.name:
                node_color[i]='white'
                node_size[i] = 0
            elif '-' in nodes['nodes'][n.name]['sequence'][(loc-1):(loc+length)]:
                node_color[i] = 'red'
                node_size[i] = 20

    edge_color = ['blue' for _ in tree_nx.edges()]
    for i,e in enumerate(tree_nx.edges()):
        if '-' in nodes['nodes'][e[1].name]['sequence'][(loc-1):(loc+length)]:
            edge_color[i] = 'red'
    plt.figure(figsize=(5, 5))
    nx.draw(tree_nx, pos=pos, node_size=node_size, alpha=0.5, node_color = node_color, edge_color=edge_color, nodelist=nodelist, with_labels=False, arrows=False)
    plt.axis('equal')
    plt.title('%.0f-bp deletion at %.0f ' % (length, loc))
    plt.savefig(TREE_DIR + ('/%.0f_%.0fbp.jpg' % (length, loc)), bbox_inches='tight')
    plt.show()