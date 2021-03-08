# Filter Raw Reads
import pandas as pd
import sys

RAW_READS_METADATA = sys.argv[1]
RAW_READS_METADATA_FILT = sys.argv[2]

raw_read_metadata = pd.read_csv(
    RAW_READS_METADATA)
raw_read_metadata = raw_read_metadata[raw_read_metadata.Bases.values<1e9]
raw_read_metadata.index = [2+i for i in range(len(raw_read_metadata))]
raw_read_metadata[['Run','Platform', 'LibraryLayout', 'Assay Type']].to_csv(RAW_READS_METADATA_FILT, sep='\t')


