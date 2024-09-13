import argparse
import glob
from Bio.PDB import PDBParser, Superimposer, PDBIO

# Set up argument parsing
parser = argparse.ArgumentParser(description="Align multiple PDB files and save the aligned structures into a single PDB file.")
parser.add_argument('-i', '--input', required=True, help="Input wildcard pattern for PDB files (e.g., '*_00*.pdb')")
parser.add_argument('-o', '--output', required=True, help="Output file name for the aligned PDB structures (e.g., 'aligned.pdb')")
args = parser.parse_args()

# Use glob to match the input pattern and get the list of PDB files
pdb_files = sorted(glob.glob(args.input))
if not pdb_files:
    raise ValueError("No PDB files found matching the pattern.")

# Initialize PDB parser
pdb_parser = PDBParser()

# Parse the first structure to use as reference
ref_structure = pdb_parser.get_structure('reference', pdb_files[0])
ref_model = ref_structure[0]  # Assuming only one model per PDB file

# Get atoms from the reference structure (assuming we're aligning CA atoms)
ref_atoms = [atom for atom in ref_model.get_atoms() if atom.name == 'CA']

# Initialize Superimposer object
super_imposer = Superimposer()

# Initialize PDBIO for saving
pdb_io = PDBIO()

model_counter =1

# Align and save
with open(args.output, "w") as outfile:
    for pdb_file in pdb_files:
        # Parse structure
        structure = pdb_parser.get_structure('target', pdb_file)
        model = structure[0]  # Assuming only one model per PDB file

        if pdb_file == pdb_files[0]:
            # For the first file (reference), no need to align
            target_atoms = ref_atoms
        else:
            # Get target atoms for alignment
            target_atoms = [atom for atom in model.get_atoms() if atom.name == 'CA']

            # Align to reference
            super_imposer.set_atoms(ref_atoms, target_atoms)
            super_imposer.apply(model.get_atoms())

        # Write MODEL record
        outfile.write(f"MODEL        {model_counter}\n")

        # Save aligned structure
        pdb_io.set_structure(structure)
        pdb_io.save(outfile)

        # Write ENDMDL record
        outfile.write("ENDMDL\n")

        model_counter += 1
        
print(f"Aligned structures have been saved to {args.output}")
