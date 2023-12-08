import os
import requests
import time
from random import random
import time
from argparse import ArgumentParser, Namespace, FileType

parser = ArgumentParser()
parser.add_argument('--pdb_path', type=str, default='data/dummy_data/1a0q_protein.pdb', help='Path to the protein .pdb file')
parser.add_argument('--ligand', type=str, default='COc(cc1)ccc1C#N', help='Either a SMILES string or the path to a molecule file that rdkit can read')
parser.add_argument('--tmp_dir', type=str, default='mydock', help='Name of the results folder')
parser.add_argument('--out_dir', type=str, default='results/user_inference', help='Directory where the outputs will be written to')
parser.add_argument('--out_name', type=str, default='mydock', help='Name of the results folder')
parser.add_argument('--gpu_id', type=str, default='0', help='set GPU to run on, (0-7)')
args = parser.parse_args()

# Set the current GPU 
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id)

# creating pair specific tmp dir for fasta and ESM embeddings
pair_tmp_dir = os.path.join(args.tmp_dir, args.out_name)
esm_embeddings_path = os.path.join(pair_tmp_dir, "esm2_output")
esm_prepped_fa = os.path.join(pair_tmp_dir, f'{args.out_name}.fa')

os.makedirs(args.out_dir, exist_ok=True)
os.makedirs(args.tmp_dir, exist_ok=True)
os.makedirs(pair_tmp_dir, exist_ok=True)
os.makedirs(esm_embeddings_path, exist_ok=True)

# The line `os.chdir("/mnt/nvme0/jbok/docking/DiffDock")` changes the current working directory to the
# specified path. In this case, it changes the current working directory to
# "/mnt/nvme0/jbok/docking/DiffDock".
os.chdir("/mnt/nvme0/jbok/docking/DiffDock")



# Prepping fasta for ESM embedding
command_prep_fa = f"python datasets/esm_embedding_preparation.py " \
        f"--protein_path {args.pdb_path} " \
        f"--out_file {esm_prepped_fa}"

# Execute the command
os.system(command_prep_fa)

os.environ['HOME'] = "esm/model_weights"
os.environ['PYTHONPATH'] = os.environ.get('PYTHONPATH', '') + ':/mnt/nvme0/jbok/docking/DiffDock/esm'

# Embed Fasta file
command_esm_embed = f"python esm/scripts/extract.py esm2_t33_650M_UR50D " \
        f"{esm_prepped_fa} " \
        f"{esm_embeddings_path} " \
        f"--repr_layers 33 " \
        f"--include per_tok " \
        f"--truncation_seq_length 30000"

# Execute the command
os.system(command_esm_embed)

# Running DiffDock
command_diffdock = f"python -m inference " \
        f"--protein_path {args.pdb_path} " \
        f"--ligand \"{args.ligand}\" " \
        f"--esm_embeddings_path {esm_embeddings_path} " \
        f"--out_dir {args.out_dir} " \
        f"--out_name {args.out_name} " \
        f"--inference_steps 20 " \
        f"--samples_per_complex 10 " \
        f"--batch_size 6"

# Execute the command
os.system(command_diffdock)
