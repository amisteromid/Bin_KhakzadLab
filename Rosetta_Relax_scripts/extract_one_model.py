from Bio.PDB import MMCIFParser, PDBIO
from sys import argv

pdb_id=argv[1]

# Replace 'your_file.cif' with the path to your CIF file
cif_file_path = f"{pdb_id}"
# Define the output PDB file name
output_pdb_path = f'{pdb_id[:4]}.pdb'

# Initialize the parser
parser = MMCIFParser()

# Parse the structure from the CIF file
structure = parser.get_structure('NMR_Structure', cif_file_path)

# Get the 5th model (index 4) from the structure
model_5 = structure[4]

# Initialize PDBIO and set the structure to model_5
io = PDBIO()
io.set_structure(model_5)

# Save the 5th model as a PDB file
io.save(output_pdb_path)

print(f"The 5th model has been saved as '{output_pdb_path}'")
