#!/bin/bash
#SBATCH --job-name=tree
#SBATCH -p dpwall
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/tree.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/tree.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15

# sh file is in $MY_SCRATCH/SARS-CoV-2_Microdeletions/src/tree/tree.sh #

ml python/3.6.1


#####################################
### CHANGE DEPENDING ON DATASET #####
VERSION_NAME=gisaid_hcov-19_2020_06_03_22
######################################


######### Set up file names #########
PARENT_PATH=/scratch/groups/dpwall/personal/briannac/SARS-CoV-2_Microdeletions
cd $PARENT_PATH
RESULTS=results/tree

REF=data/reference_genomes/NC_045512.2.fasta
ALIGNMENT=results/alignment/${VERSION_NAME}.realigned.fasta
ALIGNMENT_DEDUP=$RESULTS/${VERSION_NAME}.dedup.fasta
ALIGNMENT_MASK=$RESULTS/${VERSION_NAME}.masked.fasta

TREE=$RESULTS/${VERSION_NAME}.tree.tmp.nwk
TREE_REFINED=$RESULTS/${VERSION_NAME}.tree.nwk
NODES=$RESULTS/${VERSION_NAME}.nodes.tmp.json
NODES_REFINED=$RESULTS/${VERSION_NAME}.nodes.json
REF_DICT=results/alignment/${VERSION_NAME}.alignment_to_ref_pos_dict.npy


# Remove duplicate sequences and renumber duplicate IDs.
echo "Removing duplicates..."
python3 src/tree/remove_duplicate_seqs.py ${ALIGNMENT} NC_045512.2 > ${ALIGNMENT_DEDUP}


############################
### Mask weird sites #######
############################
echo "Masking sites..."
python3 -u src/tree/mask.py ${ALIGNMENT_DEDUP} ${ALIGNMENT_MASK} NC_045512.2 ${REF_DICT}


############################
######### Make Tree ########
############################
AUGUR_RECURSION_LIMIT=10000
echo "Computing Tree...."
augur tree \
    --alignment ${ALIGNMENT_MASK} \
    --output ${TREE} \
    --nthreads 10 \
    --substitution-model HKY
echo "Computed Tree"


############################
######## Refine Tree #######
############################
echo "Refining Tree...."
augur refine \
    --tree ${TREE} \
    --alignment ${ALIGNMENT_MASK} \
    --output-tree ${TREE_REFINED} \
    --output-node-data ${NODES} \
    --coalescent opt \
    --clock-filter-iqd 4 \
    --root=NC_045512.2
echo "Refined Tree"


#############################
## Compute Ancestral Nodes ##
#############################
echo "Computing ancestral nodes...."
augur ancestral \
    --tree ${TREE_REFINED} \
    --alignment ${ALIGNMENT_MASK} \
    --output-node-data ${NODES_REFINED} \
    --inference joint   
echo "Computed ancestral nodes"


echo "DONE SUCESS"


\rm $RESULTS/*.tmp*