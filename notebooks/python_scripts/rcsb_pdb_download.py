import os
import requests

def rcsb_seq_search(seq):
    params = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                "operator": "exact_match",
                "value": "Homo sapiens",
                "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name"
                }
            },
            {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.95,
                "sequence_type": "protein",
                "value": seq
                }
            }
            ]
        },
        "return_type": "entry"
    }
    endpoint_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(endpoint_url, json=params)
    if not response.ok:
        response.raise_for_status()    
    response_list = response.json()['result_set']
    return response_list
    
def rcsb_download_pdb(pdb_id, output_dir):
    # Download the PDB file of the hit protein
    pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response_pdb = requests.get(pdb_url)
    if not response_pdb.ok:
        response_pdb.raise_for_status()
    # Create a directory to save the PDB files
    os.makedirs(output_dir, exist_ok=True)
    # Write the PDB file to disk
    pdb_file = os.path.join(output_dir, f'{pdb_id}.pdb')
    with open(pdb_file, 'w') as f:
        f.write(response_pdb.text)
    print(f'PDB file saved to {pdb_file}')
