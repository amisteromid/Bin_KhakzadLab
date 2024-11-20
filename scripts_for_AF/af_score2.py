from pymol import cmd
import os
import math
import csv

# Function to calculate AUC
def calculate_auc_roc(data):
    """
    Calculate AUC ROC manually.
    The data is a list of pairs: [true_label, predicted_probability]
    """
    # Sort the data by predicted probability in descending order
    data_sorted = sorted(data, key=lambda x: x[1], reverse=True)
    
    # Initialize variables for the calculations
    tp = 0  # True Positives
    fp = 0  # False Positives
    fn = sum([1 for label, _ in data if label == 1])  # False Negatives (count of actual positive labels)
    tn = len(data) - fn  # True Negatives (count of actual negative labels)
    
    # Initialize variables to calculate AUC
    auc = 0
    prev_fpr = 0
    prev_tpr = 0
    
    for i in range(len(data_sorted)):
        label, _ = data_sorted[i]
        
        if label == 1:
            tp += 1
            fn -= 1
        else:
            fp += 1
            tn -= 1
        
        # Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
        
        # Calculate the area of the trapezoid formed by the current point and the previous point
        auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2
        
        # Update previous FPR and TPR
        prev_fpr = fpr
        prev_tpr = tpr
    
    return auc
 # Sort the predicted probabilities and corresponding true labels
    sorted_indices = sorted(range(len(y_pred_prob)), key=lambda i: y_pred_prob[i], reverse=True)
    y_true_sorted = [y_true[i] for i in sorted_indices]
    y_pred_prob_sorted = [y_pred_prob[i] for i in sorted_indices]
    
    # Initialize variables for calculating the AUC
    tp = 0  # True Positives
    fp = 0  # False Positives
    fn = sum(y_true)  # False Negatives (all positives in y_true)
    tn = len(y_true) - fn  # True Negatives (all negatives in y_true)
    
    # Calculate the total area under the ROC curve
    auc = 0
    for i in range(len(y_true)):
        if y_true_sorted[i] == 1:
            tp += 1
            fn -= 1
        else:
            fp += 1
            tn -= 1
        
        # Calculate the change in the TPR and FPR
        if i < len(y_true) - 1 and y_pred_prob_sorted[i] != y_pred_prob_sorted[i + 1]:
            tpr = tp / (tp + fn) if tp + fn > 0 else 0
            fpr = fp / (fp + tn) if fp + tn > 0 else 0
            auc += (fpr - (fp - 1) / (fp + tn)) * (tpr + (tp - 1) / (tp + fn)) / 2
    
    return auc

    
def calculate_precision(y_true, y_pred):
    """
    Calculate Precision manually
    """
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    
    if tp + fp == 0:
        return 0
    
    precision = tp / (tp + fp)
    return precision

def calculate_recall(y_true, y_pred):
    """
    Calculate Recall (Sensitivity) manually
    """
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)
    
    if tp + fn == 0:
        return 0
    
    recall = tp / (tp + fn)
    return recall

def calculate_specificity(y_true, y_pred):
    """
    Calculate Specificity manually
    """
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    
    if tn + fp == 0:
        return 0
    
    specificity = tn / (tn + fp)
    return specificity

def calculate_accuracy(y_true, y_pred):
    """
    Calculate Accuracy manually
    """
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)
    
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    return accuracy

def calculate_mcc(y_true, y_pred):
    """
    Calculate Matthews Correlation Coefficient manually
    """
    tp = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 1)
    tn = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 0)
    fp = sum(1 for t, p in zip(y_true, y_pred) if t == 0 and p == 1)
    fn = sum(1 for t, p in zip(y_true, y_pred) if t == 1 and p == 0)
    
    if (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) == 0:
        return 0
    
    mcc = (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return mcc
    
def get_interface_residues(mobile, mobile_chain, target, color, cutoff=5.0):
    """Get interface residues between a chain and other chains"""
    cmd.select("temp_near", f"{mobile} and chain {mobile_chain} and byres all within {cutoff} of {target}")
    #cmd.color(color,'temp_near')
    myspace = {'residues': []}
    cmd.iterate("temp_near", "residues.append(resi)", space=myspace)
    cmd.delete("temp_near")
    return set(myspace['residues'])

def get_sequence_identity(sel1, sel2):
    # Get sequences for both selections
    seq1 = cmd.get_fastastr(sel1).split('\n')[1]
    seq2 = cmd.get_fastastr(sel2).split('\n')[1]
    
    # Check if sequences are of equal length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length. Consider aligning them first.")
    
    # Count matching positions
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    
    # Calculate identity percentage
    identity = (matches / len(seq1)) * 100
    
    return identity

def analyze_interfaces(benchmark_file, af_folder, output_file):
    # Create CSV file and write header
    with open(output_file, 'w') as f:
        f.write("pdb,ag_chain,af_chain,sequence_identity,rmsd,acc,precision,recall,specificity,mcc\n")
    
    with open(benchmark_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            try:
                pdb = row['pdb']
                cmd.remove("not polymer")
                if not os.path.exists(os.path.join(af_folder, f"folds_2024_11_03_12_38/{pdb}/fold_{pdb}_model_0.cif")): continue
                ag_chain = row['antigen_chain']
                h_chain = row['Hchain']
                l_chain = row['Lchain']
                scores=[]
                ps=[]
                for model_idx in range(5):
                    cmd.reinitialize()
                    
                    # Load structures
                    af_path = os.path.join(af_folder, f"folds_2024_11_03_12_38/{pdb}/fold_{pdb}_model_{model_idx}.cif")
                    cmd.load(af_path, f"af_{model_idx}")
                    cmd.fetch(pdb.upper())
                
                    # Initialize variables for best alignment
                    best_seq_identity = -1
                    best_chain = None
                    best_rmsd = float('inf')
                    
                    # Try matching with each chain in AF model based on sequence
                    for chain in cmd.get_chains(f"af_{model_idx}"):
                        # Select relevant portions
                        cmd.select("native_ag", f"{pdb.upper()} and chain {ag_chain} and polymer")
                        cmd.select("pred_chain", f"af_{model_idx} and chain {chain} and polymer")
                        
                        # Calculate sequence identity
                        seq_identity = get_sequence_identity("native_ag", "pred_chain")
                        
                        if seq_identity > best_seq_identity:
                            best_seq_identity = seq_identity
                            best_chain = chain
                            # Get RMSD for the best sequence match
                            best_rmsd = cmd.super("native_ag", "pred_chain", quiet=1)[0]
                    
                    # Final superposition with best chain
                    cmd.select("native_ag", f"{pdb.upper()} and chain {ag_chain}")
                    cmd.select("pred_chain", f"af_{model_idx} and chain {best_chain}")
                    cmd.super("native_ag", "pred_chain", quiet=1)
                    cmd.remove(f"{pdb.upper()} and chain {ag_chain} and polymer")
                
                    # Get interface residue
                    #cmd.color('gray60',pdb.upper())
                    #cmd.color('lightpink',f'af_{model_idx}')
                    ground_truth = get_interface_residues(f"af_{model_idx}", best_chain, f"{pdb.upper()} and (chain {h_chain} or chain {l_chain}) and polymer", color='green')
                    predicted = get_interface_residues(f"af_{model_idx}", best_chain, f"af_{model_idx} and not chain {best_chain} and polymer", color='ruby')
                
                    # Calculate MCC
                    all_residues = list(set([atom.resi for atom in cmd.get_model(f"af_{model_idx} and chain C").atom]))
                    y_true = [1 if str(res) in ground_truth else 0 for res in all_residues]
                    y_pred = [1 if str(res) in predicted else 0 for res in all_residues]
                    ps.append(y_pred)
                    mcc = calculate_mcc(y_true, y_pred)
                    specificity = calculate_specificity(y_true, y_pred)
                    precision = calculate_precision(y_true, y_pred)
                    recall = calculate_recall(y_true, y_pred)
                    accuracy = calculate_accuracy(y_true, y_pred)
                    scores.append([accuracy,precision,recall,specificity,mcc])
                
                # Write results to CSV
                ps = [sum(metric) / len(ps) for metric in zip(*ps)]
                roc_auc = calculate_auc_roc([[y_true[i], ps[i]] for i in range(len(y_true))])
                mean_metrics = [sum(metric) / len(scores) for metric in zip(*scores)]
                mean_accuracy, mean_precision, mean_recall, mean_specificity, mean_mcc = mean_metrics
                with open(output_file, 'a') as f:
                    f.write(f"{pdb.upper()},{ag_chain},{best_chain},{best_seq_identity:.4f},{best_rmsd:.2f},{mean_accuracy:.4f},{mean_precision:.4f},{mean_recall:.4f},{mean_specificity:.4f},{mean_mcc:.4f},{roc_auc:.4f}\n")
                
            except Exception as e:
                print(f"Error processing {pdb}_{ag_chain}: {str(e)}")
                continue
            finally:
                # Clean up
                cmd.delete("all")

# Usage
benchmark_file = "/home/omokhtar/Desktop/final_atom/data/benchmarks/final_ag/0_409.csv"
af_folder = "/home/omokhtar/Desktop/final_atom/data/benchmarks/final_ag/af"
output_file = "interface_analysis_results.csv"
analyze_interfaces(benchmark_file, af_folder, output_file)
