from pymol.cgo import *
from pymol import cmd
from chempy import cpv
import itertools

# Setup - clean and display CA atoms
cmd.remove("not chain A")
cmd.remove("not polymer")
cmd.hide("cartoon", "all")
cmd.show("sphere", "name CA")
cmd.color("lightblue", "name CA")
cmd.set("sphere_scale", 0.5)

def calculate_pseudo_cb(chain, resi, bond_length=1.54):
    """
    Calculate pseudo CB position for glycine using the formula:
    CB' = CA - bond_length * ||(N-CA) + (C-CA)|| / ((N-CA) + (C-CA))
    
    Parameters:
    - chain: chain identifier
    - resi: residue number
    - bond_length: CA-CB bond length in Angstroms (typical ~1.54Å)
    
    Returns:
    - pseudo_cb_coords: [x, y, z] coordinates or None if calculation fails
    """
    # Get atom coordinates
    ca_coords = cmd.get_atom_coords(f"chain {chain} and resi {resi} and name CA")
    n_coords = cmd.get_atom_coords(f"chain {chain} and resi {resi} and name N")
    c_coords = cmd.get_atom_coords(f"chain {chain} and resi {resi} and name C")
    
    # Check if all atoms exist
    if ca_coords is None or n_coords is None or c_coords is None:
        return None
    
    # Calculate vectors from CA
    n_ca_vector = cpv.sub(n_coords, ca_coords)  # N - CA
    c_ca_vector = cpv.sub(c_coords, ca_coords)  # C - CA
    
    # Sum the vectors: (N-CA) + (C-CA)
    sum_vector = cpv.add(n_ca_vector, c_ca_vector)
    
    # Check for zero vector (shouldn't happen in real structures)
    if cpv.length(sum_vector) < 0.001:
        return None
    
    # Normalize the sum vector
    normalized_sum = cpv.normalize(sum_vector)
    
    # Calculate pseudo CB position: CB' = CA - bond_length * normalized_sum
    pseudo_cb_coords = cpv.sub(ca_coords, cpv.scale(normalized_sum, bond_length))
    
    return pseudo_cb_coords

def cgo_arrow(atom1, atom2, radius=0.15, arrow_length=None, hlength=0.8, hradius=0.25, color='lightblue'):
    """
    Create a CGO arrow between two atoms with improved proportions
    
    Parameters:
    - radius: cylinder thickness
    - arrow_length: total length of arrow (None = use actual CA-CB distance)
    - hlength: arrowhead length
    - hradius: arrowhead radius
    """
    xyz1 = cmd.get_atom_coords(atom1)
    xyz2 = cmd.get_atom_coords(atom2)
    
    # Handle missing atoms gracefully
    if xyz1 is None or xyz2 is None:
        return []
    
    original_vector = cpv.sub(xyz2, xyz1)
    original_length = cpv.length(original_vector)
    
    # Skip very short vectors
    if original_length < 0.1:
        return []
    
    # Use specified arrow length or actual distance
    if arrow_length is None:
        final_length = original_length
        xyz2_final = xyz2
    else:
        final_length = arrow_length
        # Scale the vector to desired length
        normalized_vector = cpv.normalize(original_vector)
        xyz2_final = cpv.add(xyz1, cpv.scale(normalized_vector, final_length))
    
    # Auto-adjust arrowhead size based on vector length if needed
    if hlength < 0:
        hlength = min(radius * 5.0, final_length * 0.3)
    if hradius < 0:
        hradius = hlength * 0.8
    
    # Ensure arrowhead doesn't exceed vector length
    hlength = min(hlength, final_length * 0.4)
    
    # Calculate where cylinder ends and cone begins
    direction = cpv.normalize(cpv.sub(xyz2_final, xyz1))
    xyz2_shift = cpv.sub(xyz2_final, cpv.scale(direction, hlength))
    
    # Get color tuple
    color_rgb = cmd.get_color_tuple(color)
    
    obj = [
        # Cylinder shaft
        CYLINDER, *xyz1, *xyz2_shift, radius, 
        *color_rgb, *color_rgb,
        
        # Cone arrowhead
        CONE, *xyz2_shift, *xyz2_final, hradius, 0.0, 
        *color_rgb, *color_rgb, 1.0, 0.0
    ]
    
    return obj

def cgo_arrow_to_coords(xyz1, xyz2, radius=0.15, arrow_length=None, hlength=0.8, hradius=0.25, color='lightblue'):
    """
    Create a CGO arrow between two coordinate points
    
    Parameters:
    - xyz1, xyz2: coordinate tuples [x, y, z]
    - radius: cylinder thickness
    - arrow_length: total length of arrow (None = use actual distance)
    - hlength: arrowhead length
    - hradius: arrowhead radius
    """
    # Handle missing coordinates gracefully
    if xyz1 is None or xyz2 is None:
        return []
    
    original_vector = cpv.sub(xyz2, xyz1)
    original_length = cpv.length(original_vector)
    
    # Skip very short vectors
    if original_length < 0.1:
        return []
    
    # Use specified arrow length or actual distance
    if arrow_length is None:
        final_length = original_length
        xyz2_final = xyz2
    else:
        final_length = arrow_length
        # Scale the vector to desired length
        normalized_vector = cpv.normalize(original_vector)
        xyz2_final = cpv.add(xyz1, cpv.scale(normalized_vector, final_length))
    
    # Auto-adjust arrowhead size based on vector length if needed
    if hlength < 0:
        hlength = min(radius * 5.0, final_length * 0.3)
    if hradius < 0:
        hradius = hlength * 0.8
    
    # Ensure arrowhead doesn't exceed vector length
    hlength = min(hlength, final_length * 0.4)
    
    # Calculate where cylinder ends and cone begins
    direction = cpv.normalize(cpv.sub(xyz2_final, xyz1))
    xyz2_shift = cpv.sub(xyz2_final, cpv.scale(direction, hlength))
    
    # Get color tuple
    color_rgb = cmd.get_color_tuple(color)
    
    obj = [
        # Cylinder shaft
        CYLINDER, *xyz1, *xyz2_shift, radius, 
        *color_rgb, *color_rgb,
        
        # Cone arrowhead
        CONE, *xyz2_shift, *xyz2_final, hradius, 0.0, 
        *color_rgb, *color_rgb, 1.0, 0.0
    ]
    
    return obj

def ca_cb_arrows(selection="all", color="lightblue", radius=0.15, arrow_length=None, hlength=0.8, hradius=0.25, glycine_color="orange"):
    """
    Create CA->CB arrows for all residues, including pseudo CB for glycine
    
    Parameters:
    - selection: PyMOL selection string
    - color: arrow color for non-glycine residues
    - radius: cylinder radius (matches sphere scale better)
    - arrow_length: total arrow length (None = use actual CA-CB distance)
    - hlength: arrowhead length 
    - hradius: arrowhead radius
    - glycine_color: special color for glycine pseudo CB arrows
    """
    model = cmd.get_model(f"({selection}) and name CA")
    arrow_count = 0
    glycine_count = 0
    
    for a in model.atom:
        resi = a.resi
        chain = a.chain
        
        sel_ca = f"chain {chain} and resi {resi} and name CA"
        sel_cb = f"chain {chain} and resi {resi} and name CB"
        
        # Check if CB atom exists (not glycine)
        if cmd.count_atoms(sel_cb) == 1:
            # Regular residue with CB atom
            arrow_obj = cgo_arrow(sel_ca, sel_cb, 
                                radius=radius, 
                                arrow_length=arrow_length,
                                hlength=hlength, 
                                hradius=hradius, 
                                color=color)
            
            if arrow_obj:  # Only load if arrow was created successfully
                cmd.load_cgo(arrow_obj, f"arrow_{chain}_{resi}")
                arrow_count += 1
        else:
            # Check if this is glycine (no CB atom)
            resn = cmd.get_model(f"chain {chain} and resi {resi} and name CA").atom[0].resn
            if resn == "GLY":
                # Calculate pseudo CB for glycine
                ca_coords = cmd.get_atom_coords(sel_ca)
                pseudo_cb_coords = calculate_pseudo_cb(chain, resi)
                
                if ca_coords is not None and pseudo_cb_coords is not None:
                    arrow_obj = cgo_arrow_to_coords(ca_coords, pseudo_cb_coords,
                                                  radius=radius,
                                                  arrow_length=arrow_length,
                                                  hlength=hlength,
                                                  hradius=hradius,
                                                  color=glycine_color)
                    
                    if arrow_obj:
                        cmd.load_cgo(arrow_obj, f"arrow_gly_{chain}_{resi}")
                        arrow_count += 1
                        glycine_count += 1
    
    print(f"Created {arrow_count} CA->CB arrows ({glycine_count} glycine pseudo CB arrows)")

def create_ca_edges(selection="all", distance_cutoff=9.0, color="gray50", radius=0.08):
    """
    Create black edges between CA atoms within specified distance
    
    Parameters:
    - selection: PyMOL selection string
    - distance_cutoff: maximum distance in Angstrom for creating edges
    - color: edge color (default black)
    - radius: cylinder radius for edges
    """
    model = cmd.get_model(f"({selection}) and name CA")
    edge_count = 0
    
    # Get all CA atoms with their coordinates
    ca_atoms = []
    for a in model.atom:
        resi = a.resi
        chain = a.chain
        sel_ca = f"chain {chain} and resi {resi} and name CA"
        coords = cmd.get_atom_coords(sel_ca)
        if coords is not None:
            ca_atoms.append({
                'selection': sel_ca,
                'coords': coords,
                'resi': resi,
                'chain': chain
            })
    
    # Get color tuple
    color_rgb = cmd.get_color_tuple(color)
    
    # Create edges between all pairs within distance cutoff
    for i, atom1 in enumerate(ca_atoms):
        for j, atom2 in enumerate(ca_atoms[i+1:], i+1):  # Avoid duplicates
            distance = cpv.distance(atom1['coords'], atom2['coords'])
            
            if distance <= distance_cutoff:
                # Create cylinder edge
                edge_obj = [
                    CYLINDER, 
                    *atom1['coords'], *atom2['coords'], 
                    radius,
                    *color_rgb, *color_rgb
                ]
                
                cmd.load_cgo(edge_obj, f"edge_{atom1['chain']}_{atom1['resi']}_{atom2['chain']}_{atom2['resi']}")
                edge_count += 1
    
    print(f"Created {edge_count} CA-CA edges within {distance_cutoff}Å")

# Optional: Create arrows with fixed lengths
def create_uniform_arrows(length=2.0, color="lightblue", glycine_color="orange"):
    """Create all arrows with the same length regardless of actual CA-CB distance"""
    cmd.delete("arrow_*")  # Remove existing arrows
    ca_cb_arrows(arrow_length=length, color=color, glycine_color=glycine_color)
    print(f"Created uniform arrows of length {length}")

def create_scaled_arrows(scale_factor=1.5, color="lightblue", glycine_color="orange"):
    """Create arrows scaled by a factor of the actual CA-CB distance"""
    cmd.delete("arrow_*")
    # This would need a more complex implementation to scale existing distances
    # For now, just create uniform arrows
    create_uniform_arrows(length=2.0 * scale_factor, color=color, glycine_color=glycine_color)

# Optional: Alternative color schemes and styling with length control
def apply_color_scheme(scheme="default", arrow_length=None):
    """Apply different color schemes to arrows with optional length control"""
    cmd.delete("arrow_*")  # Clear existing arrows
    
    if scheme == "rainbow":
        ca_cb_arrows(color="rainbow", glycine_color="magenta", arrow_length=arrow_length)
    elif scheme == "raspberry":
        ca_cb_arrows(color="raspberry", glycine_color="raspberry", radius=0.12, hlength=0.6, hradius=0.2, arrow_length=arrow_length)
    elif scheme == "green":
        ca_cb_arrows(color="green", glycine_color="yellow", radius=0.18, hlength=1.0, hradius=0.3, arrow_length=arrow_length)
    else:  # default
        ca_cb_arrows(color="lightblue", glycine_color="orange", radius=0.15, hlength=0.8, hradius=0.25, arrow_length=arrow_length)

def clear_edges():
    """Remove all CA-CA edges"""
    cmd.delete("edge_*")
    print("Cleared all CA-CA edges")

def clear_arrows():
    """Remove all CA-CB arrows"""
    cmd.delete("arrow_*")
    print("Cleared all CA-CB arrows")

def show_glycine_positions():
    """Visualize calculated pseudo CB positions for glycine residues"""
    # Create pseudo atoms for visualization
    model = cmd.get_model("name CA")
    for a in model.atom:
        resi = a.resi
        chain = a.chain
        
        # Check if this is glycine
        resn = cmd.get_model(f"chain {chain} and resi {resi} and name CA").atom[0].resn
        if resn == "GLY":
            pseudo_cb_coords = calculate_pseudo_cb(chain, resi)
            if pseudo_cb_coords is not None:
                # Create a pseudo atom for visualization
                cmd.pseudoatom(f"pseudo_cb_{chain}_{resi}", pos=pseudo_cb_coords)
    
    cmd.show("sphere", "pseudo_cb_*")
    cmd.color("orange", "pseudo_cb_*")
    cmd.set("sphere_scale", 0.3, "pseudo_cb_*")
    print("Created pseudo CB atoms for glycine residues")

# Execute the functions
ca_cb_arrows()  # Create CA->CB arrows including glycine pseudo CB
create_ca_edges(distance_cutoff=6.0)  # Create black edges between CA atoms within 6Å

# Optional: Apply red color scheme with specific arrow length
apply_color_scheme("raspberry", arrow_length=2.5)  # Red arrows, all 2.5 units long

# Optional: Show pseudo CB positions for glycine
#show_glycine_positions()

# Note: You can adjust parameters:
# - create_ca_edges(distance_cutoff=8.0, color="gray", radius=0.05)  # Thinner gray edges at 8Å
# - create_ca_edges(distance_cutoff=12.0, color="blue", radius=0.10)  # Blue edges at 12Å
# - Glycine arrows are colored differently (orange by default) to distinguish them
