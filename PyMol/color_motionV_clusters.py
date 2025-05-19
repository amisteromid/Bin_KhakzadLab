'''
dynamic_correlation_clustering /home/omokhtar/Desktop/data/DynaRepo_PDBbind/2G9H_A.pdb, force_k=5
'''
from pymol import cmd
import numpy as np
from scipy.cluster.vq import kmeans, vq
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import os

def dynamic_correlation_clustering(pdb_file, selection="name CA", max_clusters=10, outdir="./", force_k=0, min_clusters=2):
    """
    Calculate dynamic correlation vectors between CA atoms from multiple conformations in a PDB file,
    cluster them using the elbow method, and color each cluster differently in PyMOL.
    
    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file containing multiple models/conformations
    selection : str, optional
        PyMOL selection string to specify atoms for analysis (default: "name CA")
    max_clusters : int, optional
        Maximum number of clusters to test in the elbow method (default: 10)
    outdir : str, optional
        Directory to save output files (default: current directory)
    force_k : int, optional
        Force a specific number of clusters (overrides automatic detection) (default: 0 = auto)
    min_clusters : int, optional
        Minimum number of clusters to consider (default: 2)
    
    Returns:
    --------
    None
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Load the PDB file
    cmd.delete("all")
    cmd.load(pdb_file, "protein")
    
    # Get number of states/models
    n_states = cmd.count_states("protein")
    
    if n_states < 10:
        print(f"Warning: PDB file contains only {n_states} models, but at least 10 are recommended")
        if n_states < 2:
            print("Error: At least 2 models are required for dynamic analysis")
            return
    
    # Select atoms for analysis (typically CA atoms)
    cmd.select("atoms_to_analyze", selection)
    n_atoms = cmd.count_atoms("atoms_to_analyze")
    
    print(f"Analyzing {n_atoms} atoms across {n_states} conformations")
    
    # Extract coordinates for each state
    coords_by_state = []
    for state in range(1, n_states + 1):
        coords = cmd.get_coords("atoms_to_analyze", state)
        coords_by_state.append(coords)
    
    # Convert to numpy array for easier manipulation
    coords_array = np.array(coords_by_state)  # shape: (n_states, n_atoms, 3)
    
    # Calculate mean position for each atom across all states
    mean_positions = np.mean(coords_array, axis=0)  # shape: (n_atoms, 3)
    
    # Calculate movement vectors (displacement from mean)
    movement_vectors = coords_array - mean_positions  # shape: (n_states, n_atoms, 3)
    
    # Normalize movement vectors
    norms = np.linalg.norm(movement_vectors, axis=2, keepdims=True)
    # Avoid division by zero
    norms[norms == 0] = 1.0
    normalized_movement_vectors = movement_vectors / norms
    
    # Calculate Motion-V (pairwise cross products of movement vectors)
    n_pairs = n_atoms * (n_atoms - 1) // 2
    motion_v = np.zeros((n_pairs, 1))
    pair_index = 0
    
    # Initialize the full motion_v_matrix
    motion_v_matrix = np.zeros((n_atoms, n_atoms))
    
    # Compute pairwise cross-products
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            # Calculate average cross product across all states
            cross_products = []
            for state in range(n_states):
                v_i = normalized_movement_vectors[state, i]
                v_j = normalized_movement_vectors[state, j]
                cross_product = np.cross(v_i, v_j)
                cross_product_magnitude = np.linalg.norm(cross_product)
                cross_products.append(cross_product_magnitude)
            
            # Average cross product magnitude across all states
            avg_cross_product = np.mean(cross_products)
            motion_v[pair_index] = avg_cross_product
            
            # Store in the matrix (symmetric)
            motion_v_matrix[i, j] = avg_cross_product
            motion_v_matrix[j, i] = avg_cross_product
            
            pair_index += 1
    
    # Save the motion_v_matrix as a CSV file
    np.savetxt(f"{outdir}/motion_v_matrix.csv", motion_v_matrix, delimiter=",")
    
    # Find optimal number of clusters using elbow method
    distortions = []
    silhouette_scores = []
    K_range = range(min_clusters, min(max_clusters + 1, n_atoms))
    
    for k in K_range:
        # Run k-means
        centroids, _ = kmeans(motion_v, k)
        
        # Assign each data point to a cluster
        cluster_indices, _ = vq(motion_v, centroids)
        
        # Calculate distortion (within-cluster sum of squares)
        distortion = sum(np.linalg.norm(motion_v[i] - centroids[cluster_indices[i]]) ** 2 
                         for i in range(len(motion_v)))
        distortions.append(distortion)
        
        # Calculate silhouette score if there are enough samples
        if len(np.unique(cluster_indices)) > 1:
            score = silhouette_score(motion_v, cluster_indices)
            silhouette_scores.append(score)
        else:
            silhouette_scores.append(0)
    
    # Plot elbow curve
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(K_range, distortions, 'bo-')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Distortion')
    plt.title('Elbow Method for Optimal k')
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(K_range, silhouette_scores, 'ro-')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Method for Optimal k')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(f"{outdir}/elbow_method.png")
    
    # Find optimal k using both methods
    if int(force_k) > 0:
        # If user specified a forced number of clusters, use that
        optimal_k = int(force_k)
        print(f"Using user-specified k={optimal_k}")
    else:
        # Elbow method - find point of maximum curvature
        deltas = np.diff(distortions)
        if len(deltas) > 1:  # Need at least 2 points to calculate second derivative
            delta_deltas = np.diff(deltas)
            elbow_k = K_range[np.argmax(np.abs(delta_deltas)) + 1]
        else:
            elbow_k = K_range[0] if distortions else min_clusters
        
        # Silhouette method - find maximum score
        if silhouette_scores:
            silhouette_k = K_range[np.argmax(silhouette_scores)]
        else:
            silhouette_k = min_clusters
        
        # Default to higher number to avoid bias toward 2 clusters
        optimal_k = max(min_clusters + 1, round((elbow_k + silhouette_k) / 2))
        
        # If both methods suggest min_clusters and user didn't set min_clusters higher, respect that
        if elbow_k == min_clusters and silhouette_k == min_clusters and min_clusters == 2:
            optimal_k = min_clusters
            print(f"WARNING: Both methods strongly suggest only {min_clusters} clusters.")
            print("If you want more clusters, use force_k parameter to override")
        
        print(f"Elbow method suggests {elbow_k} clusters")
        print(f"Silhouette method suggests {silhouette_k} clusters")
    
    print(f"Using {optimal_k} clusters")
    
    # Perform k-means clustering with optimal k
    try:
        centroids, _ = kmeans(motion_v, optimal_k)
        cluster_indices, _ = vq(motion_v, centroids)
        
        # Debug information
        unique_clusters = np.unique(cluster_indices)
        print(f"Requested {optimal_k} clusters, got {len(unique_clusters)} clusters: {unique_clusters}")
        cluster_counts = [np.sum(cluster_indices == i) for i in range(optimal_k)]
        print(f"Cluster sizes: {cluster_counts}")
        
        # If kmeans didn't produce the requested number of clusters, try again with different initialization
        if len(unique_clusters) < optimal_k:
            print("Warning: K-means didn't produce the requested number of clusters. Trying again with different initialization...")
            for _ in range(5):  # Try a few times with different initializations
                centroids, _ = kmeans(motion_v, optimal_k, iter=50, thresh=1e-5)
                cluster_indices, _ = vq(motion_v, centroids)
                unique_clusters = np.unique(cluster_indices)
                if len(unique_clusters) >= optimal_k * 0.8:  # Accept if we get at least 80% of requested clusters
                    break
            
            print(f"After retrying, got {len(unique_clusters)} clusters")
            
            # If still unsuccessful, use a different approach
            if len(unique_clusters) < optimal_k * 0.8:
                print("Still not enough clusters. Using hierarchical clustering instead...")
                from scipy.cluster.hierarchy import linkage, fcluster
                Z = linkage(motion_v, method='ward')
                cluster_indices = fcluster(Z, optimal_k, criterion='maxclust') - 1  # Convert to 0-based indexing
    
    except Exception as e:
        print(f"Error during clustering: {e}")
        print("Falling back to simpler clustering method...")
        # Fallback to a simpler method if k-means has issues
        from sklearn.cluster import KMeans
        kmeans_model = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
        cluster_indices = kmeans_model.fit_predict(motion_v)
    
    # Create a helper function to directly extract cluster indices from motion vectors
    def extract_cluster_assignments(motion_v_matrix, n_clusters):
        """
        Directly extract cluster assignments using the full motion_v_matrix 
        rather than the pair-wise representation.
        
        Parameters:
        -----------
        motion_v_matrix : numpy.ndarray
            Square matrix of Motion-V values for all atom pairs
        n_clusters : int
            Number of clusters to extract
            
        Returns:
        --------
        numpy.ndarray
            Array of cluster assignments for each atom
        """
        from sklearn.cluster import KMeans

        # We'll use the rows of the matrix as feature vectors for clustering
        # This represents the correlation pattern of each atom with all other atoms
        kmeans_model = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        return kmeans_model.fit_predict(motion_v_matrix)
    
    # Map pair indices back to atom pairs
    if int(force_k) > 0:
        # If forcing a specific number of clusters, use the direct method
        print("Using direct atom clustering method with specified number of clusters...")
        atom_clusters = extract_cluster_assignments(motion_v_matrix, optimal_k)
    else:
        # Use the pair-based approach as before
        atom_clusters = np.zeros(n_atoms, dtype=int)
        
        # Use a different approach to assign clusters to atoms
        # We'll use the most common cluster assigned to pairs involving each atom
        cluster_counts = [dict() for _ in range(n_atoms)]
        pair_index = 0
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                cluster = cluster_indices[pair_index]
                
                # Count cluster occurrences for each atom
                if cluster not in cluster_counts[i]:
                    cluster_counts[i][cluster] = 0
                if cluster not in cluster_counts[j]:
                    cluster_counts[j][cluster] = 0
                    
                cluster_counts[i][cluster] += 1
                cluster_counts[j][cluster] += 1
                
                pair_index += 1
        
        # Assign the most frequent cluster to each atom
        for i in range(n_atoms):
            if cluster_counts[i]:
                atom_clusters[i] = max(cluster_counts[i].items(), key=lambda x: x[1])[0]
    
    # Debug: Print distribution of atoms in clusters
    unique_atom_clusters = np.unique(atom_clusters)
    print(f"\nAtom cluster distribution:")
    for c in range(optimal_k):
        count = np.sum(atom_clusters == c)
        print(f"  Cluster {c+1}: {count} atoms")
    
    # Create PyMOL selections for better visualization and debugging
    cmd.select("none")  # Clear any existing selection
    
    # First, create groups for all atoms and color them white
    cmd.create("structure", "protein", 1, 1)
    cmd.color("white", "structure")
    
    # Get atom indices from PyMOL selection
    atom_indices = []
    cmd.iterate("atoms_to_analyze", "atom_indices.append(ID)", space={"atom_indices": atom_indices})
    
    # Make sure we have colors for all clusters
    colors = [
        'red', 'green', 'blue', 'yellow', 'magenta', 'cyan', 'orange', 'pink', 
        'wheat', 'violet', 'slate', 'marine', 'olive', 'purple', 'teal', 
        'forest', 'firebrick', 'chocolate', 'chartreuse', 'splitpea',
        'lightblue', 'deepteal', 'hotpink', 'yelloworange', 'tv_green', 
        'ruby', 'tv_blue', 'lightpink', 'aquamarine', 'paleyellow'
    ]
    
    # Ensure we have enough colors by cycling if needed
    while len(colors) < optimal_k:
        colors.extend(colors[:optimal_k-len(colors)])
    
    # Show count of atoms in each cluster for debugging
    cluster_atom_counts = {}
    
    # Create selections and color them
    for cluster_id in range(optimal_k):
        # Get atom indices belonging to this cluster
        cluster_atom_indices = [str(atom_indices[i]) for i in range(n_atoms) if atom_clusters[i] == cluster_id]
        
        # Record count for debugging
        cluster_atom_counts[cluster_id] = len(cluster_atom_indices)
        
        if cluster_atom_indices:
            # Create selection string from atom indices
            selection_str = "id " + "+".join(cluster_atom_indices)
            group_name = f"cluster_{cluster_id+1}"
            
            # Create selection and color it
            cmd.select(group_name, selection_str)
            color_name = colors[cluster_id % len(colors)]
            cmd.color(color_name, group_name)
            
            # Print for debugging
            print(f"Cluster {cluster_id+1}: {len(cluster_atom_indices)} atoms, color: {color_name}")
    
    # Print cluster distribution
    print("\nCluster distribution (atoms per cluster):")
    for cluster_id, count in sorted(cluster_atom_counts.items(), key=lambda x: x[0]):
        print(f"  Cluster {cluster_id+1}: {count} atoms")
    
    # Display the structure with clusters
    cmd.hide("everything", "structure")
    cmd.show("cartoon", "structure")
    cmd.show("spheres", "atoms_to_analyze")
    
    # Create a special representation to highlight cluster boundaries
    cmd.set("sphere_scale", 0.5, "atoms_to_analyze")
    
    # Save cluster assignments with more useful information
    with open(f"{outdir}/cluster_assignments.txt", "w") as f:
        f.write("Atom_Index\tAtom_ID\tResidue\tCluster_ID\tCluster_Color\n")
        
        # Get more detailed atom information
        atom_info = []
        cmd.iterate("atoms_to_analyze", "atom_info.append((ID, resn, resi))", 
                    space={"atom_info": atom_info})
        
        for i in range(n_atoms):
            if i < len(atom_info):
                atom_id, resn, resi = atom_info[i]
                cluster_id = atom_clusters[i] + 1  # 1-based for output
                color = colors[atom_clusters[i] % len(colors)]
                f.write(f"{i+1}\t{atom_id}\t{resn}{resi}\t{cluster_id}\t{color}\n")
            else:
                f.write(f"{i+1}\tN/A\tN/A\t{atom_clusters[i]+1}\t{colors[atom_clusters[i] % len(colors)]}\n")
    
    # Save detailed motion vector information
    np.savetxt(f"{outdir}/motion_vectors.csv", motion_v, delimiter=",")
    
    # Save the actual clustering results
    cluster_info = np.column_stack((np.arange(n_atoms)+1, atom_clusters+1))
    np.savetxt(f"{outdir}/atom_clusters.csv", cluster_info, 
               delimiter=",", header="Atom_Index,Cluster_ID", 
               fmt="%d,%d", comments='')
    
    # Save PyMOL visualization state
    cmd.save(f"{outdir}/dynamic_clusters.pse")
    
    # Create a pymol script to reproduce the visualization
    with open(f"{outdir}/load_clusters.pml", "w") as f:
        f.write(f"# PyMOL script to load clusters\n")
        f.write(f"load {pdb_file}, protein\n")
        f.write(f"create structure, protein, 1, 1\n")
        f.write(f"select atoms_to_analyze, {selection}\n")
        f.write(f"hide everything\n")
        f.write(f"show cartoon, structure\n")
        f.write(f"show spheres, atoms_to_analyze\n")
        f.write(f"set sphere_scale, 0.5\n")
        
        # Add color commands for each cluster
        for cluster_id in range(optimal_k):
            color = colors[cluster_id % len(colors)]
            f.write(f"# Cluster {cluster_id+1}\n")
            
            # Get atom IDs for this cluster
            atom_ids = [str(atom_indices[i]) for i in range(n_atoms) if atom_clusters[i] == cluster_id]
            if atom_ids:
                selection_str = "id " + "+".join(atom_ids)
                f.write(f"select cluster_{cluster_id+1}, {selection_str}\n")
                f.write(f"color {color}, cluster_{cluster_id+1}\n")
    
    print(f"\nAnalysis complete. Results saved to {outdir}")
    print(f"The following files were created:")
    print(f"  - {outdir}/motion_v_matrix.csv: Raw motion correlation data")
    print(f"  - {outdir}/elbow_method.png: Plot for determining optimal cluster count")
    print(f"  - {outdir}/cluster_assignments.txt: Detailed cluster assignments for atoms")
    print(f"  - {outdir}/atom_clusters.csv: Simple atom-to-cluster mapping")
    print(f"  - {outdir}/dynamic_clusters.pse: PyMOL session with colored clusters")
    print(f"  - {outdir}/load_clusters.pml: PyMOL script to recreate visualization")
    print("\nClusters have been colored in PyMOL")
    
    return

# Add the function to PyMOL
cmd.extend("dynamic_correlation_clustering", dynamic_correlation_clustering)

# Help text
dynamic_correlation_clustering.__doc__ = """
DESCRIPTION
    Calculate dynamic correlation vectors between CA atoms from multiple conformations
    in a PDB file, cluster them using the elbow method, and color each cluster in PyMOL.

USAGE
    dynamic_correlation_clustering pdb_file [, selection [, max_clusters [, outdir [, force_k [, min_clusters]]]]]

ARGUMENTS
    pdb_file = string: path to PDB file with multiple models/conformations
    selection = string: atom selection {default: name CA}
    max_clusters = integer: maximum number of clusters to test {default: 10}
    outdir = string: directory to save output files {default: ./}
    force_k = integer: force specific number of clusters (0 = auto) {default: 0}
    min_clusters = integer: minimum number of clusters to consider {default: 2}

EXAMPLE
    dynamic_correlation_clustering /path/to/multimodel.pdb, name CA, 15, ./results
    dynamic_correlation_clustering /path/to/multimodel.pdb, name CA, 15, ./results, 5, 3

NOTES
    This function:
    1. Calculates mean positions for each selected atom across all conformations
    2. Computes movement vectors as deviations from the mean position
    3. Normalizes movement vectors
    4. Calculates Motion-V by averaging cross products of movement vectors
    5. Determines optimal cluster count using the elbow method
    6. Performs clustering and colors atoms by cluster in PyMOL
    
    Required packages: numpy, scipy, matplotlib, sklearn
"""

print("Dynamic correlation clustering function loaded successfully!")
print("Type 'help dynamic_correlation_clustering' for usage information")
print("To force a specific number of clusters, use the force_k parameter")
print("Example: dynamic_correlation_clustering pdb_file, selection, max_clusters, outdir, 5")

