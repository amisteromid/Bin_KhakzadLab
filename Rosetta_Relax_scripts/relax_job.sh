#!/bin/bash

# Activate environment
conda activate dl

# Directories and paths
rosetta_relax_dir="/home/omokhtar/Downloads/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/relax.static.linuxgccrelease"
CIF_DOWNLOAD_URL="https://files.rcsb.org/download"
ALIGN_PY_SCRIPT_PATH="/home/omokhtar/Desktop/rosetta/output/align_conformers.py"
Model_extraction_PY="extract_one_model.py"
pdb_list="/home/omokhtar/Desktop/rosetta/output/dif_ids.txt"



while read id; do
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Strating with ${id}"
  # Process directory
  cd /home/omokhtar/Desktop/rosetta/output
  # Make directory with the pdb id
  mkdir -p "$id"
  # Change into the directory
  cd "$id"
  # Download the cif file in that directory
  if [ ! -f "${id}.cif" ]; then
    wget "${CIF_DOWNLOAD_URL}/${id}.cif"
  fi
  # Run the python script
  python3 "${Model_extraction_PY}" "${id}.cif"
  # Prepare the config file
  config_file="relax_configs"
  echo "-in:file:s ${id}.pdb" > "${config_file}"
  echo "-nstruct 5" >> "${config_file}"
  echo "-in:file:fullatom" >> "${config_file}"
  echo "-relax:thorough" >> "${config_file}"
  
  $rosetta_relax_dir @"${config_file}"
  
  # align them and make a single pdb file
  python3 "${ALIGN_PY_SCRIPT_PATH}" -i "${id}_00*.pdb" -o "${id}_relaxed.pdb"
   
  
done < "$pdb_list"
