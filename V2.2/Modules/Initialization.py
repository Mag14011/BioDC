################################################################################################################################################
# Generic Modules
import os
import sys
################################################################################################################################################

################################################################################################################################################

def Initialization(LaunchDir, InputDict):

    while True:
        if ("ProgInPath" in InputDict):
            ProgInPath = InputDict["ProgInPath"]
        else:
            ProgInPath = input("""
 This program requires VMD and the Amber Molecular Dynamics Suite 
 to be in your system's PATH variable. Are they (yes/no)? """)

        print(f"ProgInPath = {ProgInPath}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if (ProgInPath == 'YES') or (ProgInPath == 'Yes') or (ProgInPath == "yes") or (ProgInPath == "Y") or (ProgInPath == "y"):
            print(f""" 
 Good! Now, here are the PDBs in the launch directory 
 ({LaunchDir}) :\n""")

            for x in os.listdir(f"{LaunchDir}"):
                if x.endswith(".pdb"):
                    print(x)

            while True:
                if ("OriginalPDB" in InputDict):
                    OriginalPDB = InputDict["OriginalPDB"] 
                else:
                    OriginalPDB = input("""
 Which PDB would you like to setup 
 (omit the .pdb file extension)? """)

                print(f"OriginalPDB = {OriginalPDB}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                if (os.path.isfile(f"{LaunchDir}/{OriginalPDB}.pdb") == True):
                    print("\n That PDB was found! ")
                    OriginalPDB = f"{LaunchDir}/{OriginalPDB}"

                    while True:
                        if ("ConsecResID" in InputDict):
                            ConsecResID = InputDict["ConsecResID"] 
                        else:
                            ConsecResID = input("\n Does the PDB have consecutive residue numbering? ")

                        print(f"ConsecResID = {ConsecResID}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                        if (ConsecResID == "YES") or (ConsecResID == "Yes") or (ConsecResID == "yes") or (ConsecResID == "Y") or (ConsecResID == "y"):
                            print(" Perfect! That is exactly what's needed to use this module.")
                            break
                        elif (ConsecResID == "NO") or (ConsecResID == "No") or (ConsecResID == "no") or (ConsecResID == "N") or (ConsecResID == "n"):
                            sys.exit("""
 To use this module, please renumber the residues 
 in the PDB to have consecutive IDs. 

 This can be done, for example, with the pdb-tools 
 at the GitHub repository 
 https://github.com/haddocking/pdb-tools 
 using the following command:
 pdb_reres original_filename.pdb > new_filename.pdb \n""")
                        else:
                            print("Sorry, I didn't understand your response")

                    return OriginalPDB
                    break
                else: 
                    print("\n That PDB does not exist unfortunately")
            break
        elif (ProgInPath == "NO") or (ProgInPath == "No") or (ProgInPath == "no") or (ProgInPath == "N") or (ProgInPath == "n"):
            sys.exit("""
 Please make VMD and the Amber Molecular Dynamics Suite 
 findalbe in your system PATH variable. Then, please 
 re-run this program \n""")
        else: 
            print("\n Sorry, I didn't understand your response.")

################################################################################################################################################
