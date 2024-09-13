import os

# Directory where your PDB files are located
directory = "/home/omokhtar/Desktop/Project/data/structures_relaxed"

# Function to process a single file
def process_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    new_lines = []
    for line in lines:
        # Remove lines that contain "END" followed by any amount of whitespace
        if line.strip() == "END":
            continue

        # Check for lines starting with "MODEL " and followed by a number greater than 9
        if line.startswith("MODEL ") and any(char.isdigit() for char in line):
            # Split the line at whitespace and check if the second part is a number > 9
            parts = line.split()
            if len(parts) > 1 and parts[1].isdigit() and int(parts[1]) > 9:
                line = line.replace("MODEL ", "MODEL", 1)
        
        new_lines.append(line)
    
    # Write the modified lines back to the file
    with open(filepath, 'w') as file:
        file.writelines(new_lines)

# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".pdb"):
        process_file(os.path.join(directory, filename))

print("Files have been processed.")

