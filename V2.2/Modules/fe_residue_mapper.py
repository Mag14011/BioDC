import os
import sys

def parse_pdb(filename):
    """Parse a PDB file to extract FE atom information."""
    fe_atoms = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                residue_number = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Check if it's an FE atom
                if atom_name == "FE":
                    fe_atoms[residue_number] = (x, y, z)
    return fe_atoms

def match_fe_atoms(pdb1_fe_atoms, pdb2_fe_atoms):
    """Match FE atoms between two PDBs based on coordinates."""
    residue_mapping = {}
    for res1, coord1 in pdb1_fe_atoms.items():
        for res2, coord2 in pdb2_fe_atoms.items():
            if coord1 == coord2:  # Compare coordinates to ensure they are the same atom
                residue_mapping[res1] = res2
                break
    return residue_mapping

def update_res_indexing(residue_mapping, input_file, output_file):
    """Update the residue numbering in the ResIndexing.txt file based on the mapping."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.split()
            if len(fields) < 6:
                outfile.write(line)
                continue
            
            original_res_num = int(fields[5])
            # Update the 6th field if a mapping exists
            new_res_num = residue_mapping.get(original_res_num, original_res_num)
            fields[5] = str(new_res_num)
            updated_line = " ".join(fields) + "\n"
            outfile.write(updated_line)

def process_files(original_pdb, reordered_pdb, indexing_dir, output_dir):
    res_indexing_file = os.path.join(indexing_dir, 'ResIndexing.txt')
    output_file = os.path.join(output_dir, "ResIndexing.txt")
    
    # Step 1: Read PDB files and extract FE atom data
    pdb1_fe_atoms = parse_pdb(original_pdb)
    pdb2_fe_atoms = parse_pdb(reordered_pdb)

    # Step 2: Match FE atoms between the two PDBs
    residue_mapping = match_fe_atoms(pdb1_fe_atoms, pdb2_fe_atoms)

    # Step 3: Update the ResIndexing.txt file
    update_res_indexing(residue_mapping, res_indexing_file, output_file)

    print(f"Updated {output_file} with the new residue numbers.")

if __name__ == '__main__':
    # Script can be called from the command line
    if len(sys.argv) != 4:
        print("Usage: python fe_residue_mapper.py <pdb1_file> <pdb2_file> <res_indexing_directory>")
        sys.exit(1)

    pdb1_file = sys.argv[1]
    pdb2_file = sys.argv[2]
    res_indexing_dir = sys.argv[3]

    process_files(pdb1_file, pdb2_file, res_indexing_dir)

