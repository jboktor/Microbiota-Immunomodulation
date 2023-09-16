#!/bin/bash
#SBATCH --nodes=1   # number of nodes
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.

source /home/${USER}/.bashrc
source activate colabfold

cd /central/groups/MazmanianLab/joeB/alphafold
colabfold_search input_sequences.fasta /path/to/db_folder msas

# ./scripts/download_all_data.sh /central/groups/MazmanianLab/joeB/Downloads/alphafold_db