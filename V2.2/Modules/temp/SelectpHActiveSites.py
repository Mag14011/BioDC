################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def SelectpHActiveSites(PDB, switch, InputDict, LaunchDir):
    print(""" 
 Residues ASP, GLU, LYS, HIS, and TYR can be titrated.""")

    if (switch == "FirstPass"):
        print(f""" 
 The residue names of ASP, GLU, and HIS residues
 need to be changed if they are to be titrated in
 CpHMD: ASP -> AS4, GLU -> GL4, and HIS -> HIP.

 We'll present a list of residue IDs for ASP, 
 GLU, and HIS residues please give a space-
 separated list of the residue IDs you would 
 like titrated  during molecular dynamics.

 Please note the following points: 
 1) N- or C-terminal residues cannot be titrated
    because the forcefield parameters are not 
    available. C-terminal residue IDs of the 
    above residue types will ONLY be excluded
    from the presented lists if an OXT atom type 
    is present in the PDB file. N-terminal 
    residues will not be automatically excluded, 
    so be careful not to select them. 
 
 2) The list of His residue IDs comprises ALL 
    His residues name HIS (not HIE, HID, or HIP), 
    but the His residues that are coordinated to 
    a heme group cannot be titrated.
        
    If you select the ID for a coordinated His, the
    program will override your choice. This feature
    allows you to blindly select all His residues, 
    and to end up with a structure with all but the 
    coordinated His residues titratable.

 Reading {PDB}.pdb ...""")

        RESLIST = ["ASP", "GLU", "HIS"]

    if (switch == "SecondPass"):
        print(f""" 
 We'll present a list of residue IDs of each type.
 For each type, please give a space-separated 
 list of the residue IDs you would like titrated 
 during molecular dynamics.

 Please note the following points: 
 1) N- or C-terminal residues cannot be titrated
    because the force field parmaeters are not 
    available. C-terminal residue IDs of the 
    above residue types will ONLY be excluded
    from the presented lists if an OXT atom type 
    is present in the PDB file. N-terminal 
    residues will not be automatically excluded, 
    so be careful not to select them. 
 
 2) The list of His residue IDs does not include
    the His residues coordinated to heme groups 
    because those residues are not titratable.

 Reading {PDB}.pdb ...""")

        RESLIST = ["AS4", "GL4", "HIP", "LYS", "TYR", "PRN"]

    SelASPIDs = ''
    SelGLUIDs = ''
    SelHISIDs = '' 
    SelLYSIDs = '' 
    SelTYRIDs = '' 
    SelPRNIDs = ''
    RESIDLIST = ''
    for RES in RESLIST:

        RESID = " "
        TERM = " "
        word1 = f"CA  {RES}"
        word2 = f"OXT {RES}"
        word3 = f"TER   "
        if (os.path.isfile(f"{PDB}.pdb") == True):

            with open(f'{PDB}.pdb', 'r') as fp:
                lines = fp.readlines()
                for line in lines:
                    if (line.find(word2) != -1) or ((line.find(word3) != -1) and (switch == "FirstPAss")):
                        idx2 = 0
                        for x in line:
                            if (idx2 >= 22) and (idx2 <= 25):
                                TERM+=x
                            idx2+=1
                        TERM+=" "
                    if (line.find(word1) != -1):
                        idx1 = 0
                        for x in line:
                            if (idx1 >= 22) and (idx1 <= 25):
                                RESID+=x
                            idx1+=1
                        RESID+=" "

            TERMLIST=list(TERM.split())
            RESIDLIST=list(RESID.split())

            for i in RESIDLIST:
                if i in TERMLIST:
                    TERMLIST.remove(i)
                    RESIDLIST.remove(i)
        if (RES != "PRN") and (len(RESIDLIST) != 0):
            print(f"\n There are {len(RESIDLIST)} {RES} residues that can be titrated with IDs: ", end=" ")
            print(*RESIDLIST)
        elif (RES == "PRN") and (switch == "SecondPass"):
            print(f"\n There are {len(RESIDLIST)} {RES} residues that can be titrated with IDs: ", end=" ")
            print(*RESIDLIST)

        if (RES == "ASP" and len(RESIDLIST) != 0) or (RES == "AS4" and len(RESIDLIST) !=0):
            if ("SelCpH" in InputDict) and ("SelASPIDs_1" in InputDict) and (switch == "FirstPass"):
                SelASPIDs = InputDict["SelASPIDs_1"]
                print(f"SelASPIDs_1 = {SelASPIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelASPIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelASPIDs = ' '.join(RESIDLIST) 
                print(f"SelASPIDs_1 = {SelASPIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelASPIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelASPIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelASPIDs_1 = {SelASPIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if ("SelCpH" in InputDict) and ("SelASPIDs_2" in InputDict) and (switch == "SecondPass"):
                SelASPIDs = InputDict["SelASPIDs_2"]
                print(f"SelASPIDs_2 = {SelASPIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelASPIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelASPIDs = ' '.join(RESIDLIST) 
                print(f"SelASPIDs_2 = {SelASPIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelASPIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelASPIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelASPIDs_2 = {SelASPIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        elif (RES == "GLU" and len(RESIDLIST) !=0) or (RES == "GL4" and len(RESIDLIST) != 0):
            if ("SelCpH" in InputDict) and ("SelGLUIDs_1" in InputDict) and (switch == "FirstPass"):
                SelGLUIDs = InputDict["SelGLUIDs_1"]
                print(f"SelGLUIDs_1 = {SelGLUIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelGLUIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelGLUIDs = ' '.join(RESIDLIST) 
                print(f"SelGLUIDs_1 = {SelGLUIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelGLUIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelGLUIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelGLUIDs_1 = {SelGLUIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if ("SelCpH" in InputDict) and ("SelGLUIDs_2" in InputDict) and (switch == "SecondPass"):
                SelGLUIDs = InputDict["SelGLUIDs_2"]
                print(f"SelGLUIDs_2 = {SelGLUIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelGLUIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelGLUIDs = ' '.join(RESIDLIST) 
                print(f"SelGLUIDs_2 = {SelGLUIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelGLUIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelGLUIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelGLUIDs_2 = {SelGLUIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        elif (RES == "HIS" and len(RESIDLIST) != 0) or (RES == "HIP" and len(RESIDLIST) !=0):
            if ("SelCpH" in InputDict) and ("SelHISIDs_1" in InputDict) and (switch == "FirstPass"):
                SelHISIDs = InputDict["SelHISIDs_1"]
                print(f"SelHISIDs_1 = {SelHISIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelHISIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelHISIDs = ' '.join(RESIDLIST) 
                print(f"SelHISIDs_1 = {SelHISIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelHISIDs_1" not in InputDict) and (switch == "FirstPass"):
                SelHISIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelHISIDs_1 = {SelHISIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if ("SelCpH" in InputDict) and ("SelHISIDs_2" in InputDict) and (switch == "SecondPass"):
                SelHISIDs = InputDict["SelHISIDs_2"]
                print(f"SelHISIDs_2 = {SelHISIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelHISIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelHISIDs = ' '.join(RESIDLIST) 
                print(f"SelHISIDs_2 = {SelHISIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelHISIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelHISIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelHISIDs_2 = {SelHISIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        elif (RES == "LYS" and len(RESIDLIST) != 0):
            if ("SelCpH" in InputDict) and ("SelLYSIDs_2" in InputDict) and (switch == "SecondPass"):
                SelLYSIDs = InputDict["SelLYSIDs_2"]
                print(f"SelLYSIDs_2 = {SelLYSIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelLYSIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelLYSIDs = ' '.join(RESIDLIST) 
                print(f"SelLYSIDs_2 = {SelLYSIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelLYSIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelLYSIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelLYSIDs_2 = {SelLYSIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        elif (RES == "TYR" and len(RESIDLIST) != 0):
            if ("SelCpH" in InputDict) and ("SelTYRIDs_2" in InputDict) and (switch == "SecondPass"):
                SelTYRIDs = InputDict["SelTYRIDs_2"]
                print(f"SelTYRIDs_2 = {SelTYRIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelTYRIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelTYRIDs = ' '.join(RESIDLIST) 
                print(f"SelTYRIDs_2 = {SelTYRIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelTRYIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelTYRIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelTYRIDs_2 = {SelTYRIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        elif (RES == "PRN") and (switch == "SecondPass"):
            if ("SelCpH" in InputDict) and ("SelPRNIDs_2" in InputDict) and (switch == "SecondPass"):
                SelPRNIDs = InputDict["SelPRNIDs_2"]
                print(f"SelPRNIDs_2 = {SelPRNIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" in InputDict) and ("SelPRNIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelPRNIDs = ' '.join(RESIDLIST) 
                print(f"SelPRNIDs_2 = {SelPRNIDs}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            elif ("SelCpH" not in InputDict) and ("SelTRYIDs_2" not in InputDict) and (switch == "SecondPass"):
                SelPRNIDs = input(f"  pH active {RES} residue IDs: ")
                print(f"SelPRNIDs_2 = {SelPRNIDs.lstrip()}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

    return SelASPIDs, SelGLUIDs, SelHISIDs, SelLYSIDs, SelTYRIDs, SelPRNIDs 
################################################################################################################################################
