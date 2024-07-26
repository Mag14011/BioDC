import os
import sys
import subprocess
from subprocess import Popen
import LambdaFromSASA
import DeltaGFromPBSA
import HemeHemeInt
import AssignCouplingFromGeom
import ComputeMarcusRates

def EnergeticEvaluation(LaunchDir, ForceFieldDir, OutPrefix, SolvEnv, PolySel, InputDict):

    while True:
        if "CompLambda" in InputDict:
            CompLambda = InputDict["CompLambda"]
            print(f"CompLambda = {CompLambda}", file=open("InteractiveInput.txt", 'a'))
        else:
            CompLambda = input("\nShould we compute the reorganization energy (yes/no)? ")
            print(f"CompLambda = {CompLambda}", file=open("InteractiveInput.txt", 'a'))

        if CompLambda.lower() in ['yes', 'y']:
            Lambda = LambdaFromSASA.LambdaFromSASA(OutPrefix)
            break
        elif CompLambda.lower() in ['no', 'n']:
            print("\nAn array where each element is the lambda for a charge transfer is needed.")
            if "NumSteps" in InputDict:
                NumSteps = int(InputDict["NumSteps"])
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))
            else:
                NumSteps = int(input("Enter the total number of charge transfer steps? "))
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))

            print("\nEnter lambda for each charge transfer step followed by return ")
            Lambda = [0] * NumSteps
            for idx in range(NumSteps):
                if f"Lambda_{idx}" in InputDict:
                    Lambda[idx] = float(InputDict[f"Lambda_{idx}"])
                    print(f"Lambda_{idx} = {Lambda[idx]}", file=open("InteractiveInput.txt", 'a'))
                else:
                    Lambda[idx] = float(input(" "))
                    print(f"Lambda_{idx} = {Lambda[idx]}", file=open("InteractiveInput.txt", 'a'))
            break
        else:
            print("Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        if "CompDG" in InputDict:
            CompDG = InputDict["CompDG"]
            print(f"CompDG = {CompDG}", file=open("InteractiveInput.txt", 'a'))
        else:
            CompDG = input("\nShould we compute the reaction free energy (yes/no)? ")
            print(f"CompDG = {CompDG}", file=open("InteractiveInput.txt", 'a'))

        if CompDG.lower() in ['yes', 'y']:
            if SolvEnv.lower() in ["explicit", "e"]:
                print("""
 Two different methods are implemented to estimate heme redox 
 potentials and thereby reaction free energies when an explicit
 solvent is present:

    (1) Compute the change in electrostatic interaction energy upon 
        heme oxidation using the Linear Interaction Energy method 
        in CPPTRAJ.

    (2) Compute the change in electrostatic interaction energy upon 
        heme oxidation using the Poisson–Boltzmann Surface Area
        method implemented in AmberTools. In this case, the explicit 
        solvent prepared with the structure is discarded.

    Recommendation:
        The LIE approach has the advantages that it is considerably 
        faster than the PBSA approach, and the overall change in 
        electrostatic energy is decomposed into contributions from 
        different groups of residues. 
	
	However, the LIE approach is only recommended for the analysis
	of multiple snapshots. Electrostatic interactions need to be 
	thermally averaged to reduce the strong dependence of fixed
	charge (non-polarizable) interactions on a given configuration.
	This averaging is effectively done in the PBSA approach, at 
	least for the solvent.

	PBSA should be used to analyze a single or multiple configurations.
	LIE should be used to analyze only multiple configurations.

	Both methods should give better results as the number of examined 
	configurations increases.\n""")

                while True:
                    if "DGmethod" in InputDict:
                        DGmethod = InputDict["DGmethod"]
                        print(f"DGmethod = {DGmethod}", file=open("InteractiveInput.txt", 'a'))
                    else:
                        DGmethod = input("Should we use the LIE or PBSA method (lie/pbsa)? ")
                        print(f"DGmethod = {DGmethod}", file=open("InteractiveInput.txt", 'a'))

                    if DGmethod.lower() in ["lie", "1"]:
                        # DG = DeltaGFromLIE(ForceFieldDir, OutPrefix, HemeType)
                        # break
                        print("Sorry, this module is only implemented in BioDC version 1.0 at the moment.")
                    elif DGmethod.lower() in ["pbsa", "2"]:
                        print(""" 
 CAUTION: Many parameters must be set to run PBSA calculations. 

 BioDC queries you for some of them, but generally adopts default 
 settings recommended in the Amber manual. It is **strongly**
 recommended to read sections 6.1 through 6.3 of the Amber manual
 and to modify the hard-coded parameter entries in the code.
 """)
                        DG, SelRefRedoxState, FFchoice = DeltaGFromPBSA.DeltaGFromPBSA(ForceFieldDir, OutPrefix, SolvEnv, PolySel)
                        break
                    else:
                        print("Sorry, I didn't understand your response.\n")
                break

            if SolvEnv.lower() in ["implicit", "i"]:
                print("""
 The method implemented to estimate heme redox potentials and thereby
 reaction free energies when an explicit solvent is NOT present is 
 to compute the change in electrostatic interaction energy upon 
 heme oxidation using the Poisson–Boltzmann Surface Area
 module in the AmberTools package (essentially a Delphi-type
 calculation). \n""")
                DG, SelRefRedoxState, FFchoice = DeltaGFromPBSA.DeltaGFromPBSA(ForceFieldDir, OutPrefix, SolvEnv, PolySel)
                break
        elif CompDG.lower() in ['no', 'n']:
            print("\nAn array where each element is the DG for a charge transfer is needed.")

            if "NumSteps" in InputDict:
                NumSteps = int(InputDict["NumSteps"])
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))
            else:
                NumSteps = int(input("Enter the total number of charge transfer steps? "))
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))

            print("\nEnter DG for each charge transfer step followed by return ")
            DG = [0] * NumSteps
            for idx in range(NumSteps):
                if f"DG_{idx}" in InputDict:
                    DG[idx] = float(InputDict[f"DG_{idx}"])
                    print(f"DG_{idx} = {DG[idx]}", file=open("InteractiveInput.txt", 'a'))
                else:
                    DG[idx] = float(input(" "))
                    print(f"DG_{idx} = {DG[idx]}", file=open("InteractiveInput.txt", 'a'))
            break
        else:
            print("Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        if "CompIntEng" in InputDict:
            CompIntEng = InputDict["CompIntEng"]
            print(f"CompIntEng = {CompIntEng}", file=open("InteractiveInput.txt", 'a'))
        else:
            CompIntEng = input("\nShould we compute heme-heme interaction energies (yes/no)? ")
            print(f"CompIntEng = {CompIntEng}", file=open("InteractiveInput.txt", 'a'))

        if CompIntEng.lower() in ['yes', 'y']:
            IntEng = HemeHemeInt.HemeHemeInt(ForceFieldDir, FFchoice, OutPrefix, SelRefRedoxState)
            break 
        elif CompIntEng.lower() in ['no', 'n']:
            print("Skipping the computation of heme-heme interaction energies.")
            break 
        else: 
            print("Sorry, I didn't understand your response.")
            
    while True:
        if "CompHda" in InputDict:
            CompHda = InputDict["CompHda"]
            print(f"CompHda = {CompHda}", file=open("InteractiveInput.txt", 'a'))
        else:
            CompHda = input("\nShould we estimate the electronic coupling from the geometry (yes/no)? ")
            print(f"CompHda = {CompHda}", file=open("InteractiveInput.txt", 'a'))

        if CompHda.lower() in ['yes', 'y']:
            Hda = AssignCouplingFromGeom.AssignCouplingFromGeom(LaunchDir, OutPrefix)
            break
        elif CompHda.lower() in ['no', 'n']:
            print("\nAn array where each element is the Hda for a charge transfer is needed.")

            if "NumSteps" in InputDict:
                NumSteps = int(InputDict["NumSteps"])
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))
            else:
                NumSteps = int(input("Enter the total number of charge transfer steps? "))
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))

            print("\nEnter Hda for each charge transfer step followed by return ")
            Hda = [0] * NumSteps
            for idx in range(NumSteps):
                if f"Hda_{idx}" in InputDict:
                    Hda[idx] = float(InputDict[f"Hda_{idx}"])
                    print(f"Hda_{idx} = {Hda[idx]}", file=open("InteractiveInput.txt", 'a'))
                else:
                    Hda[idx] = float(input(" "))
                    print(f"Hda_{idx} = {Hda[idx]}", file=open("InteractiveInput.txt", 'a'))

                if idx == 0:
                    print("Hda(Site-%0d <-> Site-%0d) Hda = %6.3f meV" % (idx, idx + 1, Hda[idx]), file=open('Hda.txt', 'w'))
                else:
                    print("Hda(Site-%0d <-> Site-%0d) Hda = %6.3f meV" % (idx, idx + 1, Hda[idx]), file=open('Hda.txt', 'a'))
            break
        else:
            print("Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        if "CompKet" in InputDict:
            CompKet = InputDict["CompKet"]
            print(f"CompKet = {CompKet}", file=open("InteractiveInput.txt", 'a'))
        else:
            CompKet = input("\nShould we compute the non-adiabatic Marcus-theory rates (yes/no)? ")
            print(f"CompKet = {CompKet}", file=open("InteractiveInput.txt", 'a'))

        if CompKet.lower() in ['yes', 'y']:
            ComputeMarcusRates.ComputeMarcusRates(Hda, Lambda, DG)
            break
        elif CompKet.lower() in ['no', 'n']:
            print("\nAn array where each element is the Ket for a charge transfer is needed.")

            if "NumSteps" in InputDict:
                NumSteps = int(InputDict["NumSteps"])
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))
            else:
                NumSteps = int(input("Enter the total number of charge transfer steps? "))
                print(f"NumSteps = {NumSteps}", file=open("InteractiveInput.txt", 'a'))

            print("\nEnter the forward Ket, return, the reverse Ket, return, and so on for each reaction ")
            ketf = [0] * NumSteps
            ketb = [0] * NumSteps
            idxf = 0
            idxb = 0
            for idx in range(2 * NumSteps):
                if idx == 0 or idx % 2 == 0:
                    if f"ketf_{idxf}" in InputDict:
                        ketf[idxf] = float(InputDict[f"ketf_{idxf}"])
                        print(f"ketf_{idxf} = {ketf[idxf]}", file=open("InteractiveInput.txt", 'a'))
                    else:
                        ketf[idxf] = float(input(" "))
                        print(f"ketf_{idxf} = {ketf[idxf]}", file=open("InteractiveInput.txt", 'a'))
                    idxf += 1
                else:
                    if f"ketb_{idxb}" in InputDict:
                        ketb[idxb] = float(InputDict[f"ketb_{idxb}"])
                        print(f"ketb_{idxb} = {ketb[idxb]}", file=open("InteractiveInput.txt", 'a'))
                    else:
                        ketb[idxb] = float(input(" "))
                        print(f"ketb_{idxb} = {ketb[idxb]}", file=open("InteractiveInput.txt", 'a'))
                    idxb += 1

            for idx in range(NumSteps):
                if idx == 0:
                    print("%.3E,%.3E" % (ketf[idx], ketb[idx]), file=open("rates.txt", "w"))
                else:
                    print("%.3E,%.3E" % (ketf[idx], ketb[idx]), file=open("rates.txt", 'a'))
            break
        else:
            print("Sorry, I didn't understand your response.")