import freesasa
freesasa.setVerbosity(freesasa.silent)
import gemmi
import tempfile
import numpy as np
'''
def calculate_sasa(model):
    # Create a temporary file to write the model PDB
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_pdb:
        # Create a new structure and add the model to this structure
        structure = gemmi.Structure()
        structure.add_model(model)  # Corrected line without adopt=True
        
        # Write the structure (with the model) to a temporary PDB file
        structure.write_pdb(temp_pdb.name)
        temp_pdb.flush()
        fs_structure = freesasa.Structure(temp_pdb.name, options={'hetatm': True})
        result = freesasa.calc(fs_structure)
        
    # Generate a list of SASA values for each atom in the model
    sasa_values = [result.atomArea(i) for i in range(fs_structure.nAtoms())]

    return sasa_values
def calculate_sasa_unbound(model):
    sasa_values_unbound = []
    atom_to_chain_map = []

    # Iterate over each chain to calculate its SASA independently
    for chain in model:
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_pdb:
            # Create a new structure
            structure = gemmi.Structure()
            
            # Create a new model with a specified name
            temp_model = gemmi.Model('temp')
            # Clone the current chain and add it to the temp_model
            cloned_chain = chain.clone()
            temp_model.add_chain(cloned_chain)
            # Add the temp_model to the structure
            structure.add_model(temp_model, pos=-1)  # pos=-1 adds the model at the end
            
            # Write the single-chain structure to a temp PDB file
            structure.write_pdb(temp_pdb.name)
            temp_pdb.flush()
            
            # Calculate SASA for the single chain
            fs_structure = freesasa.Structure(temp_pdb.name, options={'hetatm': True})
            result = freesasa.calc(fs_structure)

            # Record SASA values and the chain each atom belongs to
            for i in range(fs_structure.nAtoms()):
                sasa_values_unbound.append(result.atomArea(i))
                atom_to_chain_map.append(chain.name)

    return sasa_values_unbound


doc = gemmi.read_pdb("6wa1.pdb", max_line_length=80)
sasa_dic1={}
sasa_dic2={}
for mid, model in enumerate(doc):
	sasa_values = calculate_sasa_unbound(model)
	sasa_dic1[mid]=np.array(sasa_values)
sasa_mean1 = np.mean(np.stack(list(sasa_dic1.values())), axis=0)
for mid, model in enumerate(doc):
	sasa_values = calculate_sasa(model)
	sasa_dic2[mid]=np.array(sasa_values)
sasa_mean2 = np.mean(np.stack(list(sasa_dic2.values())), axis=0)

counter=0
for i in range(len(sasa_mean1)):
	if sasa_mean1[i] > 1: print (i,end='+')
	if abs(sasa_mean1[i] - sasa_mean2[i])>0.5:
		counter += 1
		#print (sasa_mean1[i],end='\n')
		#print (sasa_mean1[i],sasa_mean2[i], end='\n')
		print (i, end='+')
print ()
print (counter)
print (len(sasa_mean1)==len(sasa_mean2))
'''




structure = freesasa.Structure("2apn.pdb",options={'hetatm': False,'join-models': False})
result = freesasa.calc(structure)
# Print the number of atoms, including those in ligands
#print(structure.nAtoms())
residues_sasa = result.residueAreas()
#sasa_list = [sasa['Total'] for residue_number, sasa in residues_sasa['A'].items()]
if len(residues_sasa) > 1:
	raise ValueError("Structure contains more than one chain. Please specify the chain of interest.")
first_chain_id = next(iter(residues_sasa))
sasa_values = [residue_area.total for residue_area in residues_sasa[first_chain_id].values()]
print (sasa_values)
#print (residues_sasa['A']['1'].total)

'''
# Threshold for considering an atom as exposed
threshold_sasa = 5  # This is just an example; adjust as needed
surface = []
# Iterate through all atoms in the structure
for atom_index in range(structure.nAtoms()):
    # Retrieve the SASA for the current atom
    atom_sasa = result.atomArea(atom_index)
    
    # Check if the atom's SASA exceeds the threshold
    if atom_sasa > threshold_sasa:
        # Generate PyMOL command to color this atom; assuming PyMOL 'id' starts at 1
        atom_id = atom_index + 1
        print (atom_id,end='+')
        #surface.append(atom_id)
        #print(f"color blue, id {atom_id}")
'''


'''
# Define a threshold for considering an atom/residue as exposed (surface) or buried (core)
# This is a simplistic approach; you might need a more sophisticated method for precise analysis
threshold_sasa = 15.0 # Å², example threshold, adjust based on your criteria

# Iterate through residues and determine if they are surface or core based on their SASA
for chain in residues_sasa:
    print (f"color blue, chain {chain} and resi", end= ' ')
    for residue, area in residues_sasa[chain].items():
        label = "surface" if area.total > threshold_sasa else "core"
        if label == "surface":
            print (residue,end='+')
        print(f"Residue {residue} in chain {chain} is {label}: SASA = {area.total:.2f} Å²")
    print ()


# Iterate through residues and determine if they are surface or core based on their SASA
for chain in residues_sasa:
   for residue, area in residues_sasa[chain].items():
        label = "surface" if area.total > threshold_sasa else "core"
        print(f"Residue {residue} in chain {chain} is {label}: SASA = {area.total:.2f} Å²")
'''
