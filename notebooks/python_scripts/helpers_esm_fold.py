from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
import torch
from tqdm import tqdm

# Utility functions are from: https://github.com/huggingface/notebooks/blob/main/examples/protein_folding.ipynb 

def convert_outputs_to_pdb(outputs):
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs["atom37_atom_exists"]
    pdbs = []
    for i in range(outputs["aatype"].shape[0]):
        aa = outputs["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs["residue_index"][i] + 1
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs["plddt"][i],
            chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs

def esm_fold_wrapper(sequence_list, output_paths):
    # Point to tokenizer and model
    tokenizer = AutoTokenizer.from_pretrained(
        "facebook/esmfold_v1",
        cache_dir = "/central/groups/MazmanianLab/joeB/cache")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", 
        cache_dir = "/central/groups/MazmanianLab/joeB/cache",
        low_cpu_mem_usage=True)
    
    # GPU optional params
    # torch.backends.cuda.matmul.allow_tf32 = True
    model = model.cuda()
    model.esm = model.esm.half()
    tokenized_input = tokenized_input.cuda()
    
    # tokenize sequences
    seqs_tokenized = tokenizer(sequence_list, padding=False, add_special_tokens=False)['input_ids']
    outputs = []
    
    # loop through sequences and process
    with torch.no_grad():
        for input_ids in tqdm(seqs_tokenized):
            input_ids = torch.tensor(input_ids, device='cuda').unsqueeze(0)
            output = model(input_ids)
            outputs.append({key: val.cpu() for key, val in output.items()})
    
    pdb_list = [convert_outputs_to_pdb(output) for output in outputs]

    for pdb_path, pdb in zip(output_paths, pdb_list):
        with open(pdb_path, "w") as f:
            f.write("".join(pdb))

def esm_fold_wrapper_cpu(sequence_list, output_paths):
    # Point to tokenizer and model
    tokenizer = AutoTokenizer.from_pretrained(
        "facebook/esmfold_v1",
        cache_dir = "/central/groups/MazmanianLab/joeB/cache")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", 
        cache_dir = "/central/groups/MazmanianLab/joeB/cache",
        low_cpu_mem_usage=True)

    # tokenize sequences
    seqs_tokenized = tokenizer(sequence_list, padding=False, add_special_tokens=False)['input_ids']
    outputs = []
    
    # loop through sequences and process
    with torch.no_grad():
        for input_ids in tqdm(seqs_tokenized):
            input_ids = torch.tensor(input_ids).unsqueeze(0)
            output = model(input_ids)
            outputs.append({key: val.cpu() for key, val in output.items()})
    
    pdb_list = [convert_outputs_to_pdb(output) for output in outputs]

    for pdb_path, pdb in zip(output_paths, pdb_list):
        with open(pdb_path, "w") as f:
            f.write("".join(pdb))
    
