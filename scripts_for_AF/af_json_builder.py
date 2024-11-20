import json
import requests
import os

def retrieve_sequences_from_fasta(pdb_id):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
    response = requests.get(url)
    fasta_lines = response.text.split('\n')
    
    sequences = []
    current_sequence = ""
    for line in fasta_lines:
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
            current_sequence = ""
        else:
            current_sequence += line.strip()
    
    if current_sequence:
        sequences.append(current_sequence)
    
    return sequences

def generate_alphafold_job(pdb_id):
    # Retrieve sequences from RCSB
    sequences = retrieve_sequences_from_fasta(pdb_id)
    
    # Construct the AlphaFold job
    sequences_dicts = [
        {"proteinChain": {"sequence": sequence, "count": 1}}
        for sequence in sequences
    ]
    
    job = {
        "name": pdb_id,
        "modelSeeds": [],
        "sequences": sequences_dicts
    }
    
    return job

# Example usage
pdb_ids = ['7M3N', '6LGW', '6U12', '8GV7', '6W4S', '7WO5', '6Z3P', '7C01', '7LVW', '7VNG', '6WJ1', '6YIO', '7Z4T', '7CQC', '6JB8', '6OTC', '6TYS', '6ZTR', '7KF1', '6HER', '6OEJ', '6PZW', '7T77', '8GZ5', '7SOC', '7L5J', '7CJ2', '6XZF', '6VO1', '8GV6', '7T73', '6QFC', '6XC2', '6PXH', '6X97', '6OFI', '6W52', '7ZF9', '6PZ8', '7DUO', '7TFO', '7ZR7', '7SJN', '7SHY', '6P50', '7WRV', '7KET', '6ZER', '6VN0', '6ZDH', '8F8X', '6QB6', '6Q0O', '6U54', '6ZDG', '8B7W', '7MMN', '7SU1', '5ZUF', '7SJO', '7OM4', '7Z2M', '6ORN', '6HHD', '7SWN', '7MPG', '7LO6', '6Z3Q', '7EW5', '7SBG', '6ZFO', '7KF0', '6YLA', '6SV2', '7T25', '7S0E', '7TYV', '6PZZ', '7NX3', '6UUH', '6QD7', '7UED', '6WIZ', '7SU0', '7VYR', '7SGM', '6QFA', '6P4B', '7JVA', '7KFW', '7X7O', '7NP1', '6JHT', '7A5S']

# Split PDB IDs into groups of 20
pdb_groups = [pdb_ids[i:i + 20] for i in range(0, len(pdb_ids), 20)]

# Generate JSON files for each group
for idx, pdb_group in enumerate(pdb_groups, start=1):
    alphafold_jobs = [generate_alphafold_job(pdb_id) for pdb_id in pdb_group]
    
    # Save each group to a separate JSON file
    output_file = f"af_abag_jobs{idx}.json"
    with open(output_file, "w") as f:
        json.dump(alphafold_jobs, f, indent=2)

print(f"AlphaFold jobs saved to {os.path.join(os.getcwd(), output_file)}")
