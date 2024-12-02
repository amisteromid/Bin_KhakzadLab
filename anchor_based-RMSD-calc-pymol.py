from pymol import cmd

cmd.reinitialize()
cmd.load("2ajw.pdb")
cmd.load("pred1.pdb")

def align_and_calc_rmsd(mol1="2ajw", mol2="pred1", selection_og="resi 17-22", selection_pred="resi 17-22"):
    cmd.select("anchor_og", f"{mol1} and not ({selection_og})")
    cmd.select("anchor_pred", f"{mol2} and not ({selection_pred})")
    
    align_result = cmd.align("anchor_pred", "anchor_og")
    
    cmd.select("linker_og", f"{mol1} and {selection_og} and backbone")
    cmd.select("linker_pred", f"{mol2} and {selection_pred} and backbone")
    
    rmsd = cmd.rms_cur("linker_pred", "linker_og", matchmaker=4)
    
    print(f"Alignment RMSD (non-linker regions): {align_result[0]:.2f}")
    print(f"Number of atoms aligned (non-linker regions): {align_result[1]}")
    print(f"RMSD of linker region: {rmsd:.2f}")
    
    cmd.delete("anchor_og")
    cmd.delete("anchor_pred")
    cmd.delete("linker_og")
    cmd.delete("linker_pred")
    cmd.delete("pred_copy")
    
    return rmsd

align_and_calc_rmsd("2ajw", "pred1", "resi 17-22", "resi 17-22")
