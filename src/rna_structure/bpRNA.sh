

cd /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations
cp results/rna_structure/mxfold/mxfold.ct intermediate_files/rna_structure/mxfold.bpseq
cp /oak/stanford/groups/dpwall/computeEnvironments/bpRNA/bpRNA.pl intermediate_files/rna_structure/bpRNA.pl
cd intermediate_files/rna_structure
perl bpRNA.pl mxfold.bpseq
cp mxfold.st  /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/results/rna_structure/bpRNA/mxfold.st


cd /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations
awk -v OFS="\t" '$1=$1' results/rna_structure/rna_fold/SARS-CoV-2_RNAfold_folded.ct  | tail -n +2 | cut -f1,2,5  > intermediate_files/rna_structure/rnafold.bpseq
cp /oak/stanford/groups/dpwall/computeEnvironments/bpRNA/bpRNA.pl intermediate_files/rna_structure/bpRNA.pl
cd intermediate_files/rna_structure
perl bpRNA.pl rnafold.bpseq
cp rnafold.st  /home/groups/dpwall/briannac/SARS-CoV-2_structural_variations/results/rna_structure/bpRNA/rnafold.st
