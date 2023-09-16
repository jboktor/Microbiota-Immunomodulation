#!/bin/bash
#SBATCH --nodes=1   # number of nodes
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.

source /home/${USER}/.bashrc
source activate colabfold

cd /central/groups/MazmanianLab/joeB/alphafold
./scripts/download_all_data.sh /central/groups/MazmanianLab/joeB/Downloads/alphafold_db