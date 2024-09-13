from tmalign import mmalign
from sys import argv
from pymol import cmd
import os

# python3 run_pairwise_mmalign.py Relax 6wa1 /home/omokhtar/Desktop/rosetta
# python3 run_pairwise_mmalign.py NMR 6wa1 ./
# python3 run_pairwise_mmalign.py merged_relax 6wa1 /home/omokhtar/Desktop/rosetta/final_relaxed/


def NMR_var(pdb, is_available=False):
    cmd.delete("all")
    tmscores = []
    if is_available != False:
        pdb=pdb.upper()
        cmd.load(f"/home/omokhtar/Desktop/rosetta/final_relaxed/{pdb}_relaxed.pdb")
        cmd.split_states(f"{pdb}_relaxed")
    else:
        cmd.fetch(pdb)
        cmd.split_states(pdb)
    for i,n in enumerate(cmd.get_names()):
	    for ii,nn in enumerate(cmd.get_names()):
		    if ii>i:
			    s = mmalign(n,nn,exe='/home/omokhtar/MMalign')
			    tmscores.append(s)

    return tmscores
    
def Relax_var(pdb_id, directory="./"):
    cmd.delete("all")
    tmscores = []
    
    conformers = []
    # Check if the directory exists
    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
        return conformers
     
    for filename in os.listdir(directory):
        # Check if the filename starts with the pdb_id and ends with '.pdb'
        if filename.startswith(pdb_id) and filename.lower().endswith('.pdb'):
            conformers.append(os.path.join(directory, filename))
    for each in conformers:
        cmd.load(each)
    for i,n in enumerate(cmd.get_names()):
	    for ii,nn in enumerate(cmd.get_names()):
		    if ii>i:
			    s = mmalign(n,nn,exe='/home/omokhtar/MMalign')
			    tmscores.append(s)

    return tmscores

    
if __name__ == "__main__":
	stype = argv[1]
	pdb = argv[2]
	directory = argv[3]
	
	if stype == "Relax":
		print (pdb)
		print (Relax_var(pdb,directory=directory))
	elif stype == "NMR":
		print (pdb)
		print (NMR_var(pdb))
	elif stype == "merged_relax":
		print (pdb+"_relax")
		print (NMR_var(pdb,is_available=True))
		print (pdb+"_NMR")
		print (NMR_var(pdb))
