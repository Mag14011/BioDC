import os
import sys
import shutil
import subprocess
from subprocess import Popen
import StructRelax

def PreparedStructureSelection(LaunchDir, InputDict):
    if os.path.exists(f"{LaunchDir}/SPR"):
        StrucDir = f"{LaunchDir}/SPR"
    else:
        StrucDir = f"{LaunchDir}"

    print(f"\nThe following prmtop/rst7 files are in {StrucDir}:")
    for x in os.listdir(StrucDir):
        if x.endswith(".prmtop") or x.endswith(".rst7"):
            print(x)

    while True:
        if "OutPrefix" in InputDict:
            OutPrefix = InputDict["OutPrefix"]
            print(f"OutPrefix = {OutPrefix}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            OutPrefix = input("Prefix used for previously generated parm/rst7 ")
            print(f"OutPrefix = {OutPrefix}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop") and os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7"):
            print(f"Found {StrucDir}/{OutPrefix}_new.prmtop and {StrucDir}/{OutPrefix}_reord.rst7!")

            while True:
                if "SolvEnv" in InputDict:
                    SolvEnv = InputDict["SolvEnv"]
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SolvEnv = input("\nDoes your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                
                if SolvEnv.lower() in ["explicit", "implicit", "e", "i"]:
                    break
                else:
                    print("\nSorry, I didn't understand your selection for the type of solvent used.")

            if os.path.isfile(f"{StrucDir}/min.rst7"):
                print(f"Found the minimized structure ({StrucDir}/min.rst7). We are all set to proceed!")
            else:
                while True:
                    if "MinSel" in InputDict:
                        MinSel = InputDict["MinSel"]
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        MinSel = input(f"The minimized structure ({StrucDir}/min.rst7) is missing. Would you like to relax the structure (yes/no)? ")
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if MinSel.lower() in ["yes", "y"]:
                        StructRelax.StructRelax(LaunchDir, StrucDir, OutPrefix, SolvEnv, InputDict)
                        break
                    elif MinSel.lower() in ["no", "n"]:
                        sys.exit(f"Use of a structurally relaxed geometry is hard-coded into the program. If you wish to override this best-practice, please rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this module of BioDC.\n")
                    else:
                        print("Sorry, I didn't understand your response.")

            break
        elif os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop") and os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7"):
            print(f"Found {StrucDir}/{OutPrefix}_reord.prmtop and {StrucDir}/{OutPrefix}_reord.rst7!")

            while True:
                if "SolvEnv" in InputDict:
                    SolvEnv = InputDict["SolvEnv"]
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SolvEnv = input("\nDoes your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                if SolvEnv.lower() in ["explicit", "implicit", "e", "i"]:
                    break
                else:
                    print("\nSorry, I didn't understand your selection for the type of solvent used.")

            if os.path.isfile(f"{StrucDir}/min.rst7"):
                print(f"Found the minimized structure ({StrucDir}/min.rst7). We are all set to proceed!")
            else:
                while True:
                    if "MinSel" in InputDict:
                        MinSel = InputDict["MinSel"]
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        MinSel = input(f"The minimized structure ({StrucDir}/min.rst7) is missing. Would you like to relax the structure (yes/no)? ")
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if MinSel.lower() in ["yes", "y"]:
                        StructRelax.StructRelax(LaunchDir, StrucDir, OutPrefix, SolvEnv, InputDict)
                        break
                    elif MinSel.lower() in ["no", "n"]:
                        sys.exit(f"Use of a structurally relaxed geometry is hard-coded into the program. If you wish to override this best-practice, please rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this module of BioDC.\n")
                    else:
                        print("Sorry, I didn't understand your response.")

            break
        elif os.path.isfile(f"{StrucDir}/{OutPrefix}.prmtop") and os.path.isfile(f"{StrucDir}/{OutPrefix}.rst7"):
            print(f"Found {StrucDir}/{OutPrefix}.prmtop and {StrucDir}/{OutPrefix}.rst7!")

            while True:
                if "SolvEnv" in InputDict:
                    SolvEnv = InputDict["SolvEnv"]
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SolvEnv = input("\nDoes your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                    print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                if SolvEnv.lower() in ["explicit", "implicit", "e", "i"]:
                    break
                else:
                    print("\nSorry, I didn't understand your selection for the type of solvent used.")

            if os.path.isfile(f"{StrucDir}/min.rst7"):
                print(f"Found the minimized structure ({StrucDir}/min.rst7). We are all set to proceed!")
            else:
                while True:
                    if "MinSel" in InputDict:
                        MinSel = InputDict["MinSel"]
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        MinSel = input(f"The minimized structure ({StrucDir}/min.rst7) is missing. Would you like to relax the structure (yes/no)? ")
                        print(f"MinSel = {MinSel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if MinSel.lower() in ["yes", "y"]:
                        StructRelax.StructRelax(LaunchDir, StrucDir, OutPrefix, SolvEnv, InputDict)
                        break
                    elif MinSel.lower() in ["no", "n"]:
                        sys.exit(f"Use of a structurally relaxed geometry is hard-coded into the program. If you wish to override this best-practice, please rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this module of BioDC.\n")
                    else:
                        print("Sorry, I didn't understand your response.")
            break
        else:
            print("That pair of topology and coordinate files were NOT found! Please try again.")

    return OutPrefix, SolvEnv, StrucDir

def ReorderResByChain(LaunchDir, StrucDir, OutPrefix, InputDict):
    while True:
        if "PolySel" in InputDict:
            PolySel = InputDict["PolySel"]
            print(f"PolySel = {PolySel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            PolySel = input("\nIs your structure polymeric (yes/no)? ")
            print(f"PolySel = {PolySel}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if PolySel.lower() in ['yes', 'y']:
            print("""
The structure preparation stage required you to place all the heme
residues from all the chains at the end of the PDB with sequential
numbering. The programs in AmberTools that will be used to estimate
the charge transfer energetics want instead the residues to come in
the order of the connectivity; that is, the hemes of chain A should
come before any residue in chain B.

To oblige this different numbering convention, we'll use
CPPTRAJ of the AmberTools package to re-order the residues.
This process will write a new topology and coordinate file,
where the latter is of the structure you previously minimized.""")

            if os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop"):
                if os.path.isfile(f"{StrucDir}/min.rst7"):
                    print(f"\nFound the reordered topology ({StrucDir}/{OutPrefix}_new.prmtop) and coordinates ({StrucDir}/min.rst7)!")
                    OutPrefix = OutPrefix + "_new"
                    subprocess.run(f"ambpdb -p {StrucDir}/{OutPrefix}.prmtop -c {StrucDir}/min.rst7 > min.pdb", shell=True)
            elif os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop"):
                if os.path.isfile(f"{StrucDir}/min.rst7"):
                    print(f"\nFound the reordered topology ({StrucDir}/{OutPrefix}_reord.prmtop) and coordinates ({StrucDir}/min.rst7)!")
                    OutPrefix = OutPrefix + "_reord"
                    subprocess.run(f"ambpdb -p {StrucDir}/{OutPrefix}.prmtop -c {StrucDir}/{OutPrefix}.rst7 > min.pdb", shell=True)
            elif os.path.isfile(f"{StrucDir}/{OutPrefix}.prmtop"):
                print(f"""
parm {StrucDir}/{OutPrefix}.prmtop
trajin {StrucDir}/min.rst7
fixatomorder parmout {OutPrefix}_reord.prmtop
trajout {OutPrefix}_reord.rst7 topresnum
run
quit""", file=open("ReorderRes.in", "w"))

                subprocess.run("cpptraj -i ReorderRes.in > ReorderRes.log", shell=True)

                if os.path.isfile(f"{OutPrefix}_reord.prmtop"):
                    if os.path.isfile(f"{StrucDir}/min.rst7"):
                        print(f"\nCreated the reordered topology {OutPrefix}_reord.prmtop and found min.rst7!")
                        OutPrefix = OutPrefix + "_reord"
                        subprocess.run(f"ambpdb -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 > min.pdb", shell=True)
            break

        elif PolySel.lower() in ['no', 'n']:
            print("""
Great! Thanks for the clarification.
We will generate min.pdb without
re-ordering the residues.""")
            print(f"""
parm {StrucDir}/{OutPrefix}.prmtop
trajin {StrucDir}/min.rst7
trajout min.pdb
run
quit""", file=open("CreateMinPDB.in", "w"))

            subprocess.run("cpptraj -i CreateMinPDB.in > CreateMinPDB.log", shell=True)
            break
        else:
            print("Sorry, I didn't understand your response.")

    return PolySel, OutPrefix

def LinearizeHemeSequence(LaunchDir, InputDict):
    if os.path.exists(f"{LaunchDir}/SPR"):
        StrucDir = f"{LaunchDir}/SPR"
    else:
        StrucDir = f"{LaunchDir}"

    print("""
To compute the energetics for heme-to-heme electron transfer, we
need to know the linear sequence of hemes that will serve as
charge hopping sites. Typically, the linear sequence is NOT the
sequence of residue IDs in the PDB. We therefore need to specify
the linear sequence to compute the right electron transfer steps.\n""")

    if os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt"):
        shutil.copy(f"{StrucDir}/LinearizedHemeSequence.txt", f"{os.getcwd()}/LinearizedHemeSequence.txt")
        print(f"Found {StrucDir}/LinearizedHemeSequence.txt and copied it to the current (EE) directory")
    elif not os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt") and os.path.isfile(f"{StrucDir}/ResIndexing.txt"):
        shutil.copy(f"{StrucDir}/ResIndexing.txt", f"{os.getcwd()}/ResIndexing.txt")

        cidx = 0
        with open("ResIndexing.txt") as ri:
            Len_ri = len(ri.readlines())
            Entry = [0] * Len_ri
            HEM = [0] * Len_ri
            ri.seek(0)

            Lines_ri = ri.readlines()
            for line in Lines_ri:
                Entry[cidx] = line
                HEM[cidx] = int(line.strip().split(" ")[-3])
                cidx += 1

        print(f"The heme residue IDs in the present structure is/are: {HEM}\nThis sequence may not be the linear sequence in the structure, which is why we need you to enter the linear sequence.\n")

        while True:
            try:
                if "NewSequence" in InputDict:
                    New = list(map(int, InputDict["NewSequence"].strip().split()))
                    print(f"NewSequence = {' '.join(map(str, New))}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    New = list(map(int, input("Linear Sequence: ").strip().split()))
                    print(f"NewSequence = {' '.join(map(str, New))}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                x = len(New)
                RplNew = [0] * x
            except ValueError:
                print("Your entry must be a space-separated list of integers.")
            else:
                break

        for idx in range(0, x):
            if idx == 0:
                if New[idx] in HEM:
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'w'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'w'))
                else:
                    if "RplNew" in InputDict:
                        RplNew[idx] = InputDict["RplNew"]
                    else:
                        RplNew[idx] = input(f"You entered {New[idx]}, but a heme with that residue ID is not available. The available IDs are: {HEM}. You selected for analysis: {New}. What residue ID would you like in place of {New[idx]}: ")
                    New[idx] = RplNew[idx]
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'w'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'w'))
            else:
                if New[idx] in HEM:
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'a'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'a'))
                else:
                    if "RplNew" in InputDict:
                        RplNew[idx] = InputDict["RplNew"]
                    else:
                        RplNew[idx] = input(f"You entered {New[idx]}, but a heme with that residue ID is not available. The available IDs are: {HEM}. You selected for analysis: {New}. What residue ID would you like in place of {New[idx]}: ")
                    New[idx] = RplNew[idx]
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'a'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'a'))

    elif not os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt") and not os.path.isfile(f"{StrucDir}/ResIndexing.txt"):
        sys.exit(f"Both {StrucDir}/LinearizedHemeSequence.txt and {StrucDir}/ResIndexing.txt are missing. At least the latter of these files need to exist in order to proceed. Please re-run the Structure Preparation and Relaxation module to generate ResIndexing.txt and then re-run the current module.\n")
