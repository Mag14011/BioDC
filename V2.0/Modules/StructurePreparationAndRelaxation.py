################################################################################################################################################
# Generic Modules
import os
################################################################################################################################################
# Custom Modules 
import Initialization
import SelectDisulfides
import SelectMutate
import SelectpHActiveSites
import CreateResIndexing
import ProcessResIndexing
import ReBuildStructure
import StructRelax
################################################################################################################################################

def StructurePreparationAndRelaxation(LaunchDir, ForceFieldDir, InputDict):

    print("""
 We will ask a series of questions about your structure
 in order to properly prepare it. 
    """, end=" ")

    PDB = Initialization.Initialization(LaunchDir, InputDict)

    while True:
        if ("SelDisulfides" in InputDict):
            SelDisulfides = InputDict["SelDisulfides"]
        else:
            SelDisulfides = input("\n Are there disulfide linkages in your structure (yes/no)? ") 

        print(f"SelDisulfides = {SelDisulfides}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if (SelDisulfides == "YES") or (SelDisulfides == "Yes") or (SelDisulfides == "yes") or (SelDisulfides == "Y") or (SelDisulfides == "y"):
            DisulfList, DisulfArray = SelectDisulfides.SelectDisulfides(InputDict, LaunchDir)
            break
        elif (SelDisulfides == "NO") or (SelDisulfides == "No") or (SelDisulfides == "no") or (SelDisulfides == "N") or (SelDisulfides == "n"):
            DisulfList = "" 
            DisulfArray = []   
            print(" No disulfide linkages will be present in the prepared structure.")
            break
        else:
            print(" Sorry, I didn't understand your response.")

    while True:
        if ("ChooseMut" in InputDict):
            ChooseMut = InputDict["ChooseMut"]
        else:
            ChooseMut = input("\n Would you like to mutate a residue? ")

        print(f"ChooseMut = {ChooseMut}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if (ChooseMut == "YES") or (ChooseMut == "Yes") or (ChooseMut == "yes") or (ChooseMut == "Y") or (ChooseMut == "y"):
            PDB = SelectMutate.SelectMutate(PDB, LaunchDir, InputDict)
            break
        elif (ChooseMut == "NO") or (ChooseMut == "No") or (ChooseMut == "no") or (ChooseMut == "N") or (ChooseMut == "n"):
            print(" Structure will not be mutated.")
            break
        else:
            print(" Sorry, I didn't understand your response.")

    while True:
        SelASPIDs = ""
        SelGLUIDs = ""
        SelHISIDs = ""
        SelLYSIDs = ""
        SelTYRIDs = ""
        SelPRNIDs = ""

        if ("SelCpH" in InputDict):
            SelCpH = InputDict["SelCpH"]
        else:
            SelCpH = input(f""" 
 Do you intend to run molecular dynamics where
 pH active residues are titrated (yes/no)? """)

        print(f"SelCpH = {SelCpH}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if (SelCpH == "YES") or (SelCpH == "Yes") or (SelCpH == "yes") or (SelCpH == "Y") or (SelCpH == "y"):
            if ("SelASPIDs_1" in InputDict) and ("SelGLUIDs_1" in InputDict) and ("SelHISIDs_1" in InputDict):
                SelASPIDs = InputDict["SelASPIDs_1"]
                print(f"SelASPIDs_1 = {SelASPIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                SelGLUIDs = InputDict["SelGLUIDs_1"]
                print(f"SelGLUIDs_1 = {SelGLUIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                SelHISIDs = InputDict["SelHISIDs_1"]
                print(f"SelHISIDs_1 = {SelHISIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                print(PDB, "FirstPass", LaunchDir)
                SelASPIDs, SelGLUIDs, SelHISIDs, SelLYSIDs, SelTYRIDs, SelPRNIDs = SelectpHActiveSites.SelectpHActiveSites(PDB, "FirstPass", InputDict, LaunchDir)
            break

        elif (SelCpH == "NO") or (SelCpH == "No") or (SelCpH == "no") or (SelCpH == "N") or (SelCpH == "n"):
            print(""" No residues in the prepared structure will be titratable 
 if you decide to run molecular dynamics.""")
            break
        else:
            print(" Sorry, I didn't understand your response.")

    PDB = CreateResIndexing.CreateResIndexing(PDB, InputDict, LaunchDir)
    ProcessResIndexing.ProcessResIndexing(PDB, DisulfList, SelASPIDs, SelGLUIDs, SelHISIDs, SelLYSIDs, SelTYRIDs, SelPRNIDs, InputDict, LaunchDir)
    OutPrefix, SolvEnv = ReBuildStructure.ReBuildStructure(ForceFieldDir, PDB, SelCpH, InputDict, LaunchDir)
    StructRelax.StructRelax(os.getcwd(), OutPrefix, SolvEnv)

    return OutPrefix, SolvEnv
################################################################################################################################################
