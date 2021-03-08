import numpy as np
import pandas as pd
import pysam 
import sys


bam_file = sys.argv[1]
results_out = sys.argv[2]
chrom= 'chrSCV'
deletion_thresh=1
with pysam.AlignmentFile(bam_file, "rb") as samfile:
    n_reads = samfile.count()
starts_all = np.array([np.nan for _ in range(n_reads)])
ends_all = np.array([np.nan for _ in range(n_reads)])
print("n reads: ", n_reads)

with pysam.AlignmentFile(bam_file, "rb") as samfile:
    for i, read in enumerate(samfile.fetch(chrom)):
        seq = np.array([r[1] for r in read.aligned_pairs if (r[0] and r[1])])
        starts = seq[:-1]
        ends = seq[1:]
        idx = np.argmax(ends-starts)
        #if (ends[idx]-starts[idx])>deletion_thresh:
        starts_all[i] = starts[idx]
        ends_all[i] = ends[idx]
                

                

good_idx = ~np.isnan(starts_all)
starts_all = starts_all[good_idx]
ends_all = ends_all[good_idx]
df = pd.DataFrame()
df['starts'] = starts_all
df['ends'] = ends_all
df.to_csv(results_out)