#!/usr/bin/env python

from Bio import SeqIO
from collections import Counter
import sys

fasta = sys.argv[1]
filt_fasta = sys.argv[2]

seqs = [r for r in SeqIO.parse(fasta, "fasta")]

# Filter to valid sequences, and sequencs with <20 N's.
seqs_nucleotides = [s for s in seqs if len(set(s.seq.upper()).union(set(['A', 'C', 'T', 'G', 'N', '-'])))==6]
seqs_updated_nonmissing = [s for s in seqs_nucleotides if Counter(list(s.seq.upper())[100:-100])['N']<20]

# Save filtered sequences.
SeqIO.write(seqs_updated_nonmissing, filt_fasta, "fasta")
