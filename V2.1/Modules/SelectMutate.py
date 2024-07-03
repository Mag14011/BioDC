################################################################################################################################################
# Generic Modules
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def SelectMutate(PDB, LaunchDir, InputDict):

    while True:
        try:
            if ("NumMut" in InputDict):
                NumMut = int(InputDict["NumMut"])
            else:
                NumMut = int(input(" How many residues do you want to mutate? "))
        except ValueError:
                print(" Your entry must be an integer.")
        else:
            break

    print(f"NumMut = {NumMut}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

    if (NumMut != 0):
        print(f"""
 mol new {PDB}.pdb""", file=open(f"Mutate.tcl", 'w'))

        print(""" 
 For each mutation, please enter the residue three-letter code 
 before the mutaiton, the residue ID, and the three later code
 after the mutaiton, each separated by a space.

 VMD will be used to change the indicated residue names and to 
 delete the origianl sidechains. Later, TLEaP will be used to 
 build-in the new sidechains based on the appropriate template
 in the selected force field library.""")

        idx = 0
        MutResList = " "
        MutResArray = []*NumMut
        print(f"MutResArray =", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        for n in range(NumMut): 
            if ("MutResArray" in InputDict):
                MutResArray = InputDict["MutResArray"][n]
                print(f"{MutResArray[0]}-{MutResArray[1]}-{MutResArray[2]}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                SelRes = str(input(f"  Mutation {idx+1}: "))
                MutResArray = list(map(str,SelRes.split(' ')))
                print(f"{MutResArray[0]}-{MutResArray[1]}-{MutResArray[2]}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            print(f"""
 set res [atomselect top "resname {MutResArray[0]} and resid {MutResArray[1]}"]
 $res set resname {MutResArray[2]}""", file=open(f"Mutate.tcl", 'a'))

            MutResList += MutResArray[1]+" "

            idx+=1

        print(f"\n", end='', file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        print(f"""
 set mut [atomselect top "(all and not resid {MutResList}) or (not sidechain and resid {MutResList})"]
 $mut writepdb mutated.pdb

 exit
 """, file=open(f"Mutate.tcl", 'a'))


    print(f"\n Generating PDB for the mutated structure  ...")
    subprocess.run(f"vmd -e Mutate.tcl > Mutate.log", shell=True)
    PDB = "mutated" 

    return PDB

################################################################################################################################################
