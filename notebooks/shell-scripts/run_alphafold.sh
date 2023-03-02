#!/bin/bash
# Description: AlphaFold Singulatiry
# Author: Sanjay Kumar Srikakulam
# reworked for Singularity: Naveed Near-Ansari

usage() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-o <output_dir>   Path to a directory that will store the results."
        echo "-m <model_names>  Name of model to use <monomer|monomer_casp14|monomer_ptm|multimer>"
        echo "-f <fasta_path>   Path to a FASTA file containing one sequence"
        echo "-t <max_template_date> Maximum template release date to consider (ISO-8601 format - i.e. YYYY-MM-DD). Important if folding historical test sets"
        echo "Optional Parameters:"
        echo "-b <benchmark>    Run multiple JAX model evaluations to obtain a timing that excludes the compilation time, which should be more indicative of the time required for inferencing many proteins (default: 'False')"
        echo "-d <data_dir>     Path to directory of supporting data"
        echo "-g <use_gpu>      Enable NVIDIA runtime to run with GPUs (default: True)"
        echo "-a <gpu_devices>  Comma separated list of devices to pass to 'CUDA_VISIBLE_DEVICES' (default: 0)"
        echo "-n <number>       How many predictions (each with a different random seed) will be generated per model"
        echo "-p <preset>       Choose preset model configuration - no ensembling (full_dbs) or 8 model ensemblings (casp14) (default: 'full_dbs')"
        echo ""
        exit 1
}

while getopts ":d:o:m:n:f:t:a:p:g:b" i; do
        case "${i}" in
        d)
                data_dir=$OPTARG
        ;;
        o)
                output_dir=$OPTARG
        ;;
        m)
                model_names=$OPTARG
        ;;
        n)
                num_mult=$OPTARG
        ;;
        f)
                fasta_path=$OPTARG
        ;;
        t)
                max_template_date=$OPTARG
        ;;
        b)
                benchmark=true
        ;;
        g)
                use_gpu=$OPTARG
        ;;
        a)
                gpu_devices=$OPTARG
        ;;
        p)
                preset=$OPTARG
        ;;
        esac
done

# Parse input and set defaults
if [[ "$output_dir" == "" || "$model_names" == "" || "$fasta_path" == "" || "$max_template_date" == "" ]] ; then
    usage
fi


if [[ "$data_dir" == "" || "$data_dir" == "None"  ]] ; then
    data_dir=/central/software/alphafold/data/
fi

if [[ "$benchmark" == "" ]] ; then
    benchmark=false
fi

if [[ "$use_gpu" == "" ]] ; then
    use_gpu=true
fi

if [[ "$gpu_devices" == "" ]] ; then
    gpu_devices=0
fi

if [[ "$preset" == "" ]] ; then
    preset="full_dbs"
fi

if [[ "$preset" != "full_dbs" && "$preset" != "casp14" ]] ; then
    echo "Unknown preset! Using default ('full_dbs')"
    preset="full_dbs"
fi

if [[ "$model_names" == "multimer" ]] ; then
   add_dbs=" --pdb_seqres_database_path=/data/pdb_seqres/pdb_seqres.txt  --uniprot_database_path=/data/uniprot/uniprot.fasta "
fi
if [[ "$model_names" == "monomer" || "$model_names" == "monomer_casp14" || "$model_names" == "monomer_ptm" ]] ; then
    add_dbs=" --pdb70_database_path=/data/pdb70/pdb70 "
fi

if [[ "$num_mult" == "" || "$num_mult" == "None" ]] ; then
   num_mult=5
fi


# Export ENVIRONMENT variables and set CUDA devices for use
#export CUDA_VISIBLE_DEVICES=-1
#export RELAX_GPU="--nouse_gpu_relax"
if [[ "$use_gpu" == true ]] ; then
    export CUDA_VISIBLE_DEVICES=0
    export RELAX_GPU="--use_gpu_relax"

    if [[ "$gpu_devices" ]] ; then
        export CUDA_VISIBLE_DEVICES=$gpu_devices
    fi
fi

export TF_FORCE_UNIFIED_MEMORY='1'
export XLA_PYTHON_CLIENT_MEM_FRACTION='4.0'

# Path and user config (change me if required)
data_dir=/central/software/alphafold/data/
bfd_database_path="/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
mgnify_database_path="/data/mgnify/mgy_clusters.fa"
template_mmcif_dir="data/pdb_mmcif/mmcif_files"
obsolete_pdbs_path="/data/pdb_mmcif/obsolete.dat"
pdb70_database_path="/data/pdb70/"
uniclust30_database_path="/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
uniref90_database_path="/data/uniref90/uniref90.fasta"


CMD="
singularity run --nv \
 -B /central/software/alphafold/data:/data \
 -B .:/etc \
 -B /central:/central \
 --pwd  /app/alphafold /central/software/alphafold/2.2.0/container/alphafold_2.2.0.sif \
 --fasta_paths=$fasta_path  \
 --uniref90_database_path=/data/uniref90/uniref90.fasta  \
 --data_dir=/data \
 --mgnify_database_path=/data/mgnify/mgy_clusters.fa   \
 --bfd_database_path=/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
 --uniclust30_database_path=/data/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
 --template_mmcif_dir=/data/pdb_mmcif/mmcif_files  \
 --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \
 --max_template_date=$max_template_date \
 --output_dir=$output_dir  \
 --model_preset=$model_names \
 --db_preset=$preset \
 $add_dbs  \
 --num_multimer_predictions_per_model=$num_mult \
 --use_gpu_relax \
 --benchmark=$benchmark"

echo $add_num_mult
echo $CMD
$CMD
