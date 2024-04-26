################################################################################################################################################
# Generic Modules
import os
import sys
import shutil
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import StructRelax
################################################################################################################################################

################################################################################################################################################
def PreparedStructureSelection(LaunchDir):

    if (os.path.exists(f"{LaunchDir}/SPR") == True):
        StrucDir = f"{LaunchDir}/SPR"
    else:
        StrucDir = f"{LaunchDir}"

    print(f"\n The following prmtop/rst7 files are in {StrucDir}:")
    for x in os.listdir(StrucDir):
        if x.endswith(".prmtop") or x.endswith(".rst7"):
            print(x)

    while True:
        OutPrefix = input(" Prefix used for previously generated parm/rst7 ")

        if (os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop") == True) and (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7") == True):
            print(f""" 
 Found {StrucDir}/{OutPrefix}_new.prmtop and {StrucDir}/{OutPrefix}_reord.rst7! """)

            while True:
                SolvEnv = input("\n Does your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                if (SolvEnv == "Explicit") or (SolvEnv == "Implicit") or (SolvEnv == "explicit") or (SolvEnv == "implicit") or (SolvEnv == "E") or (SolvEnv == "I") or (SolvEnv == "e") or (SolvEnv == "i"):
                    pass
                    break
                else:
                    print("""\n Sorry, I didn't understand your selection for the type of solvent used.""")

            if (os.path.isfile(f"{StrucDir}/min.rst7") == True):
                    print(f""" 
 Found the minimized structure ({StrucDir}/min.rst7).
 We are all set to proceed!""")
            elif (os.path.isfile(f"{StrucDir}/min.rst7") == False):
                while True:
                    MinSel = input(f"""
 The minimized structure ({StrucDir}/min.rst7) is missing.
 Would you like to relax the structure (yes/no)? """)
                    if (MinSel == "YES") or (MinSel == "Yes") or (MinSel == "yes") or (MinSel == "Y") or (MinSel == "y"):
                        StructRelax.StructRelax(StrucDir, OutPrefix, SolvEnv)
                        break
                    elif (MinSel == "NO") or (MinSel == "No") or (MinSel == "no") or (MinSel == "N") or (MinSel == "n"):
                        sys.exit(f"""
 Use of a structually relaxed geometry is hard-corded into the 
 program. If you wish to override this best-practice, please
 rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this
 module of BioDC. \n""")
                    else:
                        print(" Sorry, I didn't understand your response.")

            break
        elif (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop") == True) and (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7") == True): 
            print(f""" 
 Found {StrucDir}/{OutPrefix}_reord.prmtop and {StrucDir}/{OutPrefix}_reord.rst7! """)

            while True:
                SolvEnv = input("\n Does your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                if (SolvEnv == "Explicit") or (SolvEnv == "Implicit") or (SolvEnv == "explicit") or (SolvEnv == "implicit") or (SolvEnv == "E") or (SolvEnv == "I") or (SolvEnv == "e") or (SolvEnv == "i"):
                    pass
                    break
                else:
                    print("""\n Sorry, I didn't understand your selection for the type of solvent used.""")

            if (os.path.isfile(f"{StrucDir}/min.rst7") == True):
                print(f"""
 Found the minimized structure ({StrucDir}/min.rst7).
 We are all set to proceed!""")
            elif (os.path.isfile(f"{StrucDir}/min.rst7") == False):
                while True:
                    MinSel = input(f"""
 The minimized structure ({StrucDir}/min.rst7) is missing.
 Would you like to relax the structure (yes/no)? """)
                    if (MinSel == "YES") or (MinSel == "Yes") or (MinSel == "yes") or (MinSel == "Y") or (MinSel == "y"):
                        StructRelax.StructRelax(StrucDir, OutPrefix, SolvEnv)
                        break
                    elif (MinSel == "NO") or (MinSel == "No") or (MinSel == "no") or (MinSel == "N") or (MinSel == "n"):
                        sys.exit(f"""
 Use of a structually relaxed geometry is hard-corded into the 
 program. If you wish to override this best-practice, please
 rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this
 module of BioDC. \n""")
                    else:
                        print(" Sorry, I didn't understand your response.")


            break
        elif (os.path.isfile(f"{StrucDir}/{OutPrefix}.prmtop") == True) and (os.path.isfile(f"{StrucDir}/{OutPrefix}.rst7") == True): 
            print(f"""
 Found {StrucDir}/{OutPrefix}.prmtop and {StrucDir}/{OutPrefix}.rst7! """)

            while True:
                SolvEnv = input("\n Does your structure have an explicit or implicit solvent present (explicit/implicit)? ")
                if (SolvEnv == "Explicit") or (SolvEnv == "Implicit") or (SolvEnv == "explicit") or (SolvEnv == "implicit") or (SolvEnv == "E") or (SolvEnv == "I") or (SolvEnv == "e") or (SolvEnv == "i"):
                    pass
                    break
                else:
                    print("""\n Sorry, I didn't understand your selection for the type of solvent used.""")

            if (os.path.isfile(f"{StrucDir}/min.rst7") == True):
                print(f"""
 Found the minimized structure ({StrucDir}/min.rst7).
 We are all set to proceed!""")
            elif (os.path.isfile(f"{StrucDir}/min.rst7") == False):
                while True:
                    MinSel = input(f"""
 The minimized structure ({StrucDir}/min.rst7) is missing.
 Would you like to relax the structure (yes/no)? """)
                    if (MinSel == "YES") or (MinSel == "Yes") or (MinSel == "yes") or (MinSel == "Y") or (MinSel == "y"):
                        StructRelax.StructRelax(StrucDir, OutPrefix, SolvEnv)
                        break
                    elif (MinSel == "NO") or (MinSel == "No") or (MinSel == "no") or (MinSel == "N") or (MinSel == "n"):
                        sys.exit(f"""
 Use of a structually relaxed geometry is hard-corded into the 
 program. If you wish to override this best-practice, please
 rename your prepared {OutPrefix}.rst7 to min.rst7 and re-run this
 module of BioDC. \n""")
                    else:
                        print(" Sorry, I didn't understand your response.")
            break
        else:
            print(""" 
 That pair of topology and coordinate files were NOT found! 
 Please try again.""")

#           sys.exit(""" 
#A pairt of topology (prmtop) and coordinate (rst7) files with that
#output-prefix is missing. 
#
#Please try running the Structure Preparation and Relaxation module 
#before proceeding""")

    return OutPrefix, SolvEnv, StrucDir
################################################################################################################################################

################################################################################################################################################

def ReorderResByChain(StrucDir, OutPrefix):

    while True:
        PolySel = input("""
 Is your structure polymeric (yes/no)? """)

        if (PolySel == 'YES') or (PolySel == 'Yes') or (PolySel == "yes") or (PolySel == "Y") or (PolySel == "y"):
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

            if (os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop") == True):
                if (os.path.isfile(f"{StrucDir}/min.rst7") == True):
                    print(f"\n Found the reordered topology ({StrucDir}/{OutPrefix}_new.prmtop) and coordinates ({StrucDir}/min.rst7)!")
                    OutPrefix = OutPrefix+"_new"
                    subprocess.run(f"ambpdb -p {StrucDir}/{OutPrefix}.prmtop -c {StrucDir}/min.rst7 > min.pdb", shell=True)
            elif (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop") == True):
                if (os.path.isfile(f"{StrucDir}/min.rst7") == True):
                    print(f"\n Found the reordered topology ({StrucDir}/{OutPrefix}_reord.prmtop) and coordinates ({StrucDir}/min.rst7)!")
                    OutPrefix = OutPrefix+"_reord"
                    subprocess.run(f"ambpdb -p {StrucDir}/{OutPrefix}.prmtop -c {StrucDir}/{OutPrefix}.rst7 > min.pdb", shell=True)
            elif (os.path.isfile(f"{StrucDir}/{OutPrefix}.prmtop") == True):
                print(f"""
parm {StrucDir}/{OutPrefix}.prmtop
trajin {StrucDir}/min.rst7
fixatomorder parmout {OutPrefix}_reord.prmtop
trajout {OutPrefix}_reord.rst7 topresnum
run
quit""", file=open("ReorderRes.in", "w"))

                subprocess.run("cpptraj -i ReorderRes.in > ReorderRes.log", shell=True)

                if (os.path.isfile(f"{OutPrefix}_reord.prmtop") == True):
                    if (os.path.isfile("{StrucDir}/min.rst7") == True):
                        print(f"\n Created the reordered topology {OutPrefix}_reord.prmtop and found min.rst7!")
                        OutPrefix = OutPrefix+"_reord"
                        subprocess.run(f"ambpdb -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 > min.pdb", shell=True)
            break

        elif (PolySel == 'NO') or (PolySel == 'No') or (PolySel == "no") or (PolySel == "N") or (PolySel == "n"):
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
            print(" Sorry, I didn't understand your response.")

    return PolySel, OutPrefix
################################################################################################################################################

################################################################################################################################################

def LinearizeHemeSequence(LaunchDir):

    if (os.path.exists(f"{LaunchDir}/SPR") == True):
        StrucDir = f"{LaunchDir}/SPR"
    else:
        StrucDir = f"{LaunchDir}"

    print("""
 To compute the energetics for heme-to-heme electron transfer, we 
 need to know the linear sequence of hemes that will serve as 
 charge hopping sites. Typically, the linear sequence is NOT the 
 sequence of residue IDs in the PDB. We therefore need to specify
 the linear sequence to compute the right electron transfer steps. \n""")

    if (os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt") == True):
        shutil.copy(f"{StrucDir}/LinearizedHemeSequence.txt", f"{os.getcwd()}/LinearizedHemeSequence.txt")
        print(f" Found {StrucDir}/LinearizedHemeSequence.txt and copied it to the current (EE) directory")
    elif (os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt") == False) and (os.path.isfile(f"{StrucDir}/ResIndexing.txt") == True):
        shutil.copy(f"{StrucDir}/ResIndexing.txt", f"{os.getcwd()}/ResIndexing.txt")

        cidx=0
        with open("ResIndexing.txt") as ri:
            Len_ri = len(ri.readlines())
            Entry = [0]*Len_ri
            HEM = [0]*Len_ri
            ri.seek(0)

            Lines_ri = ri.readlines()
            for line in Lines_ri:
                Entry[cidx] = line
                HEM[cidx] = int(line.strip().split(" ")[-3])
                cidx+=1

        print(f""" The heme residue IDs in the present structure is/are: 
 {HEM}

 This sequence may not be the linear sequence in the structure, whcih is
 why we need you to enter the linear sequence.\n""")

        while True:
            try:
                New=list(map(int, input(" Linear Sequence: ").strip().split()))
                x = len(New)
                RplNew = [0]*x
            except ValueError:
                print(" Your enter must be a space-separate list of integers.")
            else:
                break
        
        for idx in range(0, x):
            if (idx == 0):
                if ( New[idx] in HEM ): 
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'w'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'w'))
                elif ( New[idx] not in HEM ):
                    RplNew[idx] = input(f""" 
 You entered {New[idx]}, but a heme with that residue ID is not available. 

     The available IDs are: {HEM}
 You selected for analysis: {New}

 What residue ID would you like inplace of {New[idx]}: """)
                    New[idx] = RplNew[idx]
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'w'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'w'))
            else:
                if ( New[idx] in HEM ): 
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'a'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'a'))
                elif ( New[idx] not in HEM ):
                    RplNew[idx] = input(f""" 
 You entered {New[idx]}, but a heme with that residue ID is not available. 

     The available IDs are: {HEM}
 You selected for analysis: {New}

 What residue ID would you like inplace of {New[idx]}: """)
                    New[idx] = RplNew[idx]
                    print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'a'))
                    print(Entry[HEM.index(int(New[idx]))], end='', file=open('SelResIndexing.txt', 'a'))

    elif (os.path.isfile(f"{StrucDir}/LinearizedHemeSequence.txt") == False) and (os.path.isfile(f"{StrucDir}/ResIndexing.txt") == False):
        sys.exit(f"""
 Both {StrucDir}/LinearizedHemeSequence.txt and {StrucDir}/ResIndexing.txt are missing. 
 At least the latter of these files need to exist in order to proceed.
 Please re-run teh Structure Preparation and Relaxation module to
 generate ResIndexing.txt and then re-run the current module.\n""")

################################################################################################################################################
