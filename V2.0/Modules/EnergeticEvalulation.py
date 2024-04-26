################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import LambdaFromSASA
import DeltaGFromPBSA
import HemeHemeInt
import AssignCouplingFromGeom
import ComputeMarcusRates
################################################################################################################################################

def EnergeticEvalulation(LaunchDir, ForceFieldDir, OutPrefix, SolvEnv, PolySel):

    while True:
        CompLambda = input("""
 Should we compute the reorganization energy (yes/no)? """)

        if (CompLambda == 'YES') or (CompLambda == "Yes") or (CompLambda == "yes") or (CompLambda == "Y") or (CompLambda == "y"):
            Lambda = LambdaFromSASA.LambdaFromSASA(OutPrefix)
            break

        elif (CompLambda == 'NO') or (CompLambda == "No") or (CompLambda == "no") or (CompLambda == "N") or (CompLambda == "n"):
            print("""
 An array where each eleemnt is the lambda for a charge transfer is needed.""")

            NumSteps = int(input(" Enter the total number of charge transfer steps? "))
            print("""
 Enter lambda for each charge transfer step followed by return """)

            Lambda = [0]*(NumSteps)
            for idx in range(0, NumSteps):
                Lambda[idx] = float(input(""" """))
            break
        else:
            print(" Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        CompDG = input("""
 Should we compute the reaction free energy (yes/no)? """)

        if (CompDG == 'YES') or (CompDG == "Yes") or (CompDG == "yes") or (CompDG == "Y") or (CompDG == "y"):

            if (SolvEnv == "explicit") or (SolvEnv == "Explicit") or (SolvEnv == "e") or (SolvEnv == "E"):
                print("""
 Two different methods are implemented to estimate heme redox 
 potentials and thereby reaction free energies when an explicit
 solvent is present:

    (1) Compute the change in electrostatic interaction energy upon 
        heme oxidation using the Linear Interaction Energy method 
        in CPPTRAJ.

    (2) Compute the change in electrostatic interaciton energy upon 
        heme oxidation using the Poisson–Boltzmann Surface Area
        method implemented in AmberTools. In this case, the explicit 
        solvent prepared with the sturcutre is discarded.

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
                    DGmethod = input(""" Should we use the LIE or PBSA method (lie/pbsa)? """)

                    if (DGmethod == "LIE") or (DGmethod == "lie") or (DGmethod == str(1)):
                        #DG = DeltaGFromLIE(ForceFieldDir, OutPrefix, HemeType)
                        #break
                        print(" Sorry, this module is only implemented in BioDC version 1.0 at the moment.")
                    elif (DGmethod == "PBSA") or (DGmethod == "pbsa") or (DGmethod == str(2)):
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
                        print(" Sorry, I didn't understand your response.\n")
                break

            if (SolvEnv == "implicit") or (SolvEnv == "Implicit") or (SolvEnv == "i") or (SolvEnv == "I"):
                print("""
 The method implemented to estimate heme redox potentials and thereby
 reaction free energies when an explicit solvent is NOT presnet is 
 to compute the change in electrostatic interaciton energy upon 
 heme oxidation using the Poisson–Boltzmann Surface Area
 module in the AmberTools package (essentially a Delphi-type
 calculation). \n""")
                DG, SelRefRedoxState, FFchoice = DeltaGFromPBSA.DeltaGFromPBSA(ForceFieldDir, OutPrefix, SolvEnv, PolySel)
                break
        elif (CompDG == 'NO') or (CompDG == "No") or (CompDG == "no") or (CompDG == "N") or (CompDG == "n"):
            print("""
 An array where each eleemnt is the DG for a charge transfer is needed.""")

            NumSteps = int(input(" Enter the total number of charge transfer step? "))
            print("""
 Enter DG for each charge transfer step followed by return """)

            DG = [0]*(NumSteps)
            for idx in range(0, NumSteps):
                DG[idx] = float(input(""" """))
            break
        else:
            print(" Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        CompIntEng = input("""
 Should we compute heme-heme interaction energies (yes/no)? """)

        if (CompIntEng == 'YES') or (CompIntEng == "Yes") or (CompIntEng == "yes") or (CompIntEng == "Y") or (CompIntEng == "y"):
            IntEng = HemeHemeInt.HemeHemeInt(ForceFieldDir, FFchoice, OutPrefix, SelRefRedoxState)
            break 
        elif (CompIntEng == 'NO') or (CompIntEng == "No") or (CompIntEng == "no") or (CompIntEng == "N") or (CompIntEng == "n"):
            print(" Skipping the computation of heme-heme interaction energies.")
            break 
        else: 
            print(" Sorry, I didn't understand your response.")
            
    while True:
        CompHda = input("""
 Should we estimate the electronic coupling from the geometry (yes/no)? """)

        if (CompHda == 'YES') or (CompHda == "Yes") or (CompHda == "yes") or (CompHda == "Y") or (CompHda == "y"):
            Hda = AssignCouplingFromGeom.AssignCouplingFromGeom(LaunchDir, OutPrefix)
            break
        elif (CompHda == 'NO') or (CompHda == "No") or (CompHda == "no") or (CompHda == "N") or (CompHda == "n"):
            print("""
 An array where each eleemnt is the Hda for a charge transfer is needed.""")

            idx = 0
            NumSteps = int(input(" Enter the total number of charge transfer step? "))
            print("""
 Enter Hda for each charge transfer step followed by return """)

            Hda = [0]*(NumSteps)
            for idx in range(0, NumSteps):
                Hda[idx] = float(input(""" """))

                if (idx == 0):
                    print("Hda(Site-%0d <-> Site-%0d) Hda = %6.3f meV" %(idx,  idx+1, Hda[idx]), file=open('Hda.txt', 'w'))
                else:
                    print("Hda(Site-%0d <-> Site-%0d) Hda = %6.3f meV" %(idx,  idx+1, Hda[idx]), file=open('Hda.txt', 'a'))
            break
        else:
            print(" Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        CompKet = input("""
 Should we compute the non-adiabatic Marcus-theory rates (yes/no)? """)

        if (CompKet == 'YES') or (CompKet == "Yes") or (CompKet == "yes") or (CompKet == "Y") or (CompKet == "y"):
            ComputeMarcusRates.ComputeMarcusRates(Hda, Lambda, DG)
            break
        elif (CompKet == 'NO') or (CompKet == "No") or (CompKet == "no") or (CompKet == "N") or (CompKet == "n"):
            print("""
 An array where each eleemnt is the Ket for a charge transfer is needed.""")

            NumSteps = int(input(" Enter the total number of charge transfer steps? "))
            print("""
 Enter the forward Ket, return, the reverse Ket, return, and so on 
 for each reaction """)

            ketf = [0]*(NumSteps)
            ketb = [0]*(NumSteps)
            idxf = 0; idxb= 0; 
            for idx in range(0, (2*NumSteps)):
                if (idx == 0) or (idx % 2 == 0):
                    ketf[idxf] = float(input(""" """))
                    idxf += 1
                elif (idx != 0) or (idx % 2 != 0):
                    ketb[idxb] = float(input(""" """))
                    idxb += 1

            for idx in range(0, NumSteps):
                if (idx == 0):
                    print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "w"))
                else:
                    print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "a"))
            break
        else:
            print(" Sorry, I didn't understand your response.")

################################################################################################################################################
