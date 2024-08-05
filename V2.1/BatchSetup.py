import os
import shutil
import subprocess

def create_folders_and_input_files(step, nproc=None, python_program_path=None, run_mode=None):
    # Get the current working directory
    directory = os.getcwd()
    
    # Get list of PDB files in the directory
    pdb_files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    
    if not pdb_files:
        # Check subdirectories if no PDB files are found in the current directory
        for subdir in os.listdir(directory):
            subdir_path = os.path.join(directory, subdir)
            if os.path.isdir(subdir_path):
                potential_pdb_file = f"{subdir}.pdb"
                potential_pdb_path = os.path.join(subdir_path, potential_pdb_file)
                if os.path.isfile(potential_pdb_path):
                    pdb_files.append(potential_pdb_path)
    
    if not pdb_files:
        print("No PDB files found in the current directory or subdirectories.")
        return
    
    # Create folders for each PDB file
    for pdb_file in pdb_files:
        if os.path.isfile(pdb_file):
            pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
            folder_path = os.path.join(directory, pdb_name)
            
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
            
            # Move the PDB file into the corresponding folder
            shutil.move(pdb_file, os.path.join(folder_path, os.path.basename(pdb_file)))
            
            # Write the appropriate input.txt file based on the step
            input_file_path = os.path.join(folder_path, 'input.txt')
            print(f"Writing to {input_file_path} for step {step}")
            
            with open(input_file_path, 'w') as file:
                if step == 1:
                    if nproc is None:
                        print("Number of processors (NProc) not provided for step 1.")
                        return
                    file.write(f"DivSel = 1\n")
                    file.write(f"ProgInPath = y\n")
                    file.write(f"OriginalPDB = {pdb_name}\n")
                    file.write(f"ConsecResID = y\n")
                    file.write(f"SelDisulfides = n\n")
                    file.write(f"ChooseMut = n\n")
                    file.write(f"SelCpH = y\n")
                    file.write(f"CreateResIndexingMethod = auto\n")
                    file.write(f"DistMin = 2.0\n")
                    file.write(f"DistMax = 4.0\n")
                    file.write(f"ProcessPDBChoice = yes\n")
                    file.write(f"RedoxState = R O R R\n")
                    file.write(f"OutPrefix = BioDC\n")
                    file.write(f"SolvEnv = e\n")
                    file.write(f"BoxShape = rec\n")
                    file.write(f"BufferSize = 5\n")
                    file.write(f"NaCount = 0\n")
                    file.write(f"ClCount = 0\n")
                    file.write(f"StructRelaxCompChoice = P\n")
                    file.write(f"NProc = {nproc}\n")
                elif step == 2:
                    file.write(f"DivSel = 2\n")
                    file.write(f"OutPrefix = BioDC\n")
                    file.write(f"SolvEnv = e\n")
                    file.write(f"PolySel = n\n")
                    file.write(f"NewSequence = 197 200 203 206\n")
                    file.write(f"CompLambda = yes\n")
                    file.write(f"CompDG = y\n")
                    file.write(f"DGmethod = pbsa\n")
                    file.write(f"ChooseRef = red\n")
                    file.write(f"DeltaGFromPBSACompChoice = p\n")
                    file.write(f"SelEpsin0 = y\n")
                    file.write(f"epsout0 = 78.2\n")
                    file.write(f"istrng0 = 100.0\n")
                    file.write(f"memb0 = n\n")
                    file.write(f"SelDelphi0 = n\n")
                    file.write(f"SelEpsin1 = y\n")
                    file.write(f"epsout1 = 78.2\n")
                    file.write(f"istrng1 = 100.0\n")
                    file.write(f"memb1 = n\n")
                    file.write(f"SelDelphi1 = n\n")
                    file.write(f"SelEpsin2 = y\n")
                    file.write(f"epsout2 = 78.2\n")
                    file.write(f"istrng2 = 100.0\n")
                    file.write(f"memb2 = n\n")
                    file.write(f"SelDelphi2 = n\n")
                    file.write(f"SelEpsin3 = y\n")
                    file.write(f"epsout3 = 78.2\n")
                    file.write(f"istrng3 = 100.0\n")
                    file.write(f"memb3 = n\n")
                    file.write(f"SelDelphi3 = n\n")
                    file.write(f"CompIntEng = y\n")
                    file.write(f"SelEpsin = y\n")
                    file.write(f"epsout = 78.2\n")
                    file.write(f"istrng = 100.0\n")
                    file.write(f"memb = n\n")
                    file.write(f"SelDelphi = n\n")
                    file.write(f"HemeHemeIntCompChoice = p\n")
                    file.write(f"CompHda = yes\n")
                    file.write(f"CompKet = yes\n")
                else:
                    print(f"Invalid step: {step}")
                    return
            
            # Launch the external Python program in the background from within the subdirectory
            if python_program_path:
                try:
                    if run_mode == 'parallel' or run_mode == 'p':
                        subprocess.Popen(["python", python_program_path, folder_path],
                        cwd=folder_path,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        print(f"Launched {python_program_path} for folder {folder_path} in parallel mode")
                    elif run_mode == 'serial' or run_mode == 's':
                        subprocess.run(["python", python_program_path, folder_path],
                        cwd=folder_path)
                        print(f"Executed {python_program_path} for folder {folder_path} in serial mode")
                except Exception as e:
                    print(f"Failed to launch the Python program: {e}")

if __name__ == "__main__":
    print("""
 This script assists with running BioDC in batch mode.

 1) To run in batch mode, place all the properly prepared PDBs into a folder. 
 2) From within that folder, run this scirpt. 

 The script
   * Moves each PDB into its own folder
   * Writes an input file based on whether you told the script you want to run step 1 or 2.
   * Launches BioDC in each folder, either in serial or parallel depending on your choice.
   * Within each sub-folder, BioDC will detect the input.txt file and appropriately run without any intervention the calculations.
""")

    try:
        step = int(input("\n Enter step number (1 or 2): "))
        nproc = None
        if step == 1:
            nproc = int(input("Enter the number of processors (NProc) for each minimization: "))
        python_program_path = input("Enter the path to BioDC: ")
        run_mode            = input("Run the calculations in parallel (p) or serial (s): ")
        create_folders_and_input_files(step, nproc, python_program_path, run_mode)
        print("Folders and input files created successfully.")
    except ValueError as e:
        print(f"Invalid input: {e}")

