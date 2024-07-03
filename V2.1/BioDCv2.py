################################################################################################################################################
# Generic Modules
import os
import sys
import shutil
from datetime import datetime
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
#ProgDir = sys.path[0].split('BioDCv2')[1]
ProgDir = sys.path[0]
if (os.path.exists(f"{ProgDir}/Modules") == True) and (os.path.exists(f"{ProgDir}/ForceFieldLib") == True):
    sys.path.append(f"{ProgDir}/Modules")
elif (os.path.exists(f"{ProgDir}/Modules") == False) or (os.path.exists(f"{ProgDir}/ForceFieldLib") == False):
    sys.exit("""
 Either the Modules or the ForceFieldLib sub-directories
 are missing from the BioDCv2 directory. Please try 
 re-downloading BioDC from the GitHub repository to 
 fix this problem.\n""")

import StructurePreparationAndRelaxation as SPR
import AntechamberToEnergeticEvalulation as A2EE
import EnergeticEvalulation as EE
import RedoxCurrentPrediction as RCP
import ReadInput
################################################################################################################################################

ForceFieldDir = f"{ProgDir}/ForceFieldLib"
LaunchDir = os.getcwd()
TimeStamp = datetime.now()

#print("Enter your selections, one per line. Press Ctrl+D (or Ctrl+Z on Windows) when done:")
#selections = sys.stdin.read().splitlines()
#if not selections:
#    print("No selections were provided.")
#else:
#   for selection in selections:
#       print(f"You selected: {selection}")
#   pass
################################################################################################################################################

print(f"""
 ================================================================== 
                        Weclome to BioDC
             A program that automates and accelerates
               the computation of redox currents in
                (polymeric) multi-heme cytochormes 

              Written by Matthew J. Guberman-Pfeffer
                  Last Updated: 12/20/2023

 Start time: {TimeStamp}
 Directory paths:
   Program:                    {ProgDir}
   Heme forcefield parameters: {ForceFieldDir}
   Current working directory:  {LaunchDir}
 ================================================================== 
""", end=" ")

print(f"""
 BioDC presents a highly modular workflow that has three 
 major divisions: 
    (1) Structure Preparaiton & Relaxation
    (2) Energetic Estimation
    (3) Redox Current Prediction
""", end=" ")

if (os.path.isfile(f"input.txt") == True):
    InputDict = ReadInput.ReadInput()

    for key, value in InputDict.items():
        print(f"{key}: {value}")
else:
    InputDict = ''

while True:
    try:
        if ("DivSel" in InputDict):
            DivSel = int(InputDict["DivSel"])
        else:
            DivSel = int(input("""
 Which of these divisions would you like to perform?
 (Enter zero "0" to be guided through the entire
 workflow.) (0/1/2/3) """))
        break
    except ValueError:
        print("""
 Sorry, you must enter 0, 1, 2, or 3.""")
    except NameError:
        print("""
 Sorry, you must enter 0, 1, 2, or 3.""")


if (os.path.isfile(f"{LaunchDir}/InteractiveInput.txt") == True):
    shutil.move(f"{LaunchDir}/InteractiveInput.txt", f"{LaunchDir}/PriorInteractiveInput.txt") 

print(f"DivSel = {DivSel}", file=open("InteractiveInput.txt", 'w'))

if (DivSel == 0) or (DivSel == 1):
    print("""
 ===================================================================  
 First: Structure Preparation & Relaxation
 ===================================================================\n""") 

    try:
        print(""" Creating a directory called SPR for this module.""")
        os.makedirs("SPR")
        os.chdir('SPR')
        print(f""" Now, the current working directory is: {os.getcwd()}""")
    except FileExistsError:
        while True:
            if ("DirFate" in InputDict):
                DirFate = InputDict["DirFate"]
            else:
                DirFate = input("""
   A directory named SPR already exists. 
   It may have been created from a prior 
   session with BioDC. Would you like to 
   delete it and start over? """)

            print(f"DirFate = {DirFate}", file=open("InteractiveInput.txt", 'a'))

            if ( DirFate == "YES" ) or ( DirFate == "Yes" ) or ( DirFate == "yes" ) or ( DirFate == "Y" ) or ( DirFate == "y" ):
                shutil.rmtree("SPR")
                print("\n   Old directory deleted.")
                os.makedirs("SPR")
                print("   New directory Created.")
                os.chdir('SPR')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            elif ( DirFate == "NO" ) or ( DirFate == "No" ) or ( DirFate == "no" ) or ( DirFate == "N" ) or ( DirFate == "n" ):
                os.chdir('SPR')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            else:
                print("Sorry, I didn't understand your response.")

    OutPrefix, SolvEnv  = SPR.StructurePreparationAndRelaxation(LaunchDir, ForceFieldDir, InputDict)

if (DivSel == 0) or (DivSel == 2):
    if (DivSel == 0):
        os.chdir('../')

    print("""
 ===================================================================  
 Antechamber to Energetic Evalulation 
 ===================================================================\n""")

    try:
        print(""" Creating a directory called EE for this module.""")
        os.makedirs("EE")
        os.chdir('EE')
        print(f""" Now, the current working directory is: {os.getcwd()}""")
    except FileExistsError:
        while True:
            DirFate = input("""
   A directory named EE already exists. 
   It may have been created from a prior 
   session with BioDC. Would you like to 
   delete it and start over? """)
            if ( DirFate == "YES" ) or ( DirFate == "Yes" ) or ( DirFate == "yes" ) or ( DirFate == "Y" ) or ( DirFate == "y" ):
                shutil.rmtree("EE")
                print("\n   Old directory deleted.")
                os.makedirs("EE")
                print("   New directory Created.")
                os.chdir('EE')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            elif ( DirFate == "NO" ) or ( DirFate == "No" ) or ( DirFate == "no" ) or ( DirFate == "N" ) or ( DirFate == "n" ):
                os.chdir('EE')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            else:
                print("Sorry, I didn't understand your response.")

    if (DivSel == 0):
        StrucDir = f"{LaunchDir}/SPR"
        PolySel, OutPrefix = A2EE.ReorderResByChain(StrucDir, OutPrefix)
        A2EE.LinearizeHemeSequence(LaunchDir)
        print("""
 ===================================================================  
 Second: Energetic Evalulation
 ===================================================================\n""", end=" ")
        EE.EnergeticEvalulation(LaunchDir, ForceFieldDir, OutPrefix, SolvEnv, PolySel)

    elif (DivSel == 2):
        OutPrefix, SolvEnv, StrucDir = A2EE.PreparedStructureSelection(LaunchDir)
        PolySel, OutPrefix = A2EE.ReorderResByChain(StrucDir, OutPrefix)
        A2EE.LinearizeHemeSequence(LaunchDir)
        print("""
 ===================================================================  
 Second: Energetic Evalulation
 ===================================================================\n""", end=" ")
        EE.EnergeticEvalulation(LaunchDir, ForceFieldDir, OutPrefix, SolvEnv, PolySel)

if (DivSel == 0) or (DivSel == 3):
    if (DivSel == 0):
        os.chdir('../')

    print("""
 ===================================================================  
 Third: Redox Current Prediction
 ===================================================================\n""")

    try:
        print(""" Creating a directory called RCP for this module.""")
        os.makedirs("RCP")
        os.chdir('RCP')
        print(f""" Now, the current working directory is: {os.getcwd()}""")
    except FileExistsError:
        while True:
            DirFate = input("""
   A directory named RCP already exists. 
   It may have been created from a prior 
   session with BioDC. Would you like to 
   delete it and start over? """)
            if ( DirFate == "YES" ) or ( DirFate == "Yes" ) or ( DirFate == "yes" ) or ( DirFate == "Y" ) or ( DirFate == "y" ):
                shutil.rmtree("RCP")
                print("\n   Old directory deleted.")
                os.makedirs("RCP")
                print("   New directory Created.")
                os.chdir('RCP')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            elif ( DirFate == "NO" ) or ( DirFate == "No" ) or ( DirFate == "no" ) or ( DirFate == "N" ) or ( DirFate == "n" ):
                os.chdir('RCP')
                print(f"   Now, the current working directory is: {os.getcwd()}")
                break
            else:
                print("Sorry, I didn't understand your response.")

    RCP.RedoxCurrentPrediction(LaunchDir)

if (DivSel != 0) and (DivSel != 1) and (DivSel != 2) and (DivSel != 3):
    sys.exit("""
 Please re-launch BioDC and select one of the available modules. \n""")
################################################################################################################################################
