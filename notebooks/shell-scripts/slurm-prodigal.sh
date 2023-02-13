#!/bin/bash
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mail-user=jboktor@caltech.edu   # email address
#SBATCH --mail-type=FAIL   # Notify on failure.

source /home/${USER}/.bashrc
source activate pdmbsR

# Script to download any file via wget and sbatch 
# When running the command specify the URL and the output directory
# as well as the walltime (--time) and memory (--mem-per-cpu) and job name (--job-name)
# e.g.)
# sbatch \
# --time=1-00:00:00 --mem-per-cpu=100G --job-name="Download-XXX" \
# slurm_wget.sh \
# -u my-url.tar.gz -o path-to-my-output/

while getopts "i:o:g:" opt
do
    case "$opt" in 
        i ) INPUT="$OPTARG" ;;
        o ) OUTPUT="$OPTARG" ;;
        g ) GENE_COORDS="$OPTARG" ;;
    esac
done


prodigal -i "${INPUT}" -o "${GENE_COORDS}" -a "${OUTPUT}"
rm "${GENE_COORDS}"