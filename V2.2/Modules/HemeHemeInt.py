################################################################################################################################################
# Generic Modules
import os
import sys
import itertools
import subprocess
import argparse
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import DefineRefState
import PairedChargeAssignment
import ReadInput
import AnalyzeHemeCooperativity as AHC
################################################################################################################################################

def HemeHemeInt(LaunchDir, ForceFieldDir, FFchoice, OutPrefix, SelRefRedoxState, InputDict):

    if not os.path.isfile("RefState.prmtop"):
        print("""
 Generating the reference state topology where all hemes
 are in the reduced state.
 """)
        SelRefRedoxState = DefineRefState(LaunchDir, OutPrefix, InputDict)
    else:
        print("""
 Found RefState.prmtop
 """)

    Es = []
    if os.path.isfile("Lambda.txt"):
        with open('Lambda.txt', 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                if 'Es        =' in line:
                    Es.append(float(line.strip().split()[2]))

        if Es:
            min_epsin = round(min(Es), 3)
            max_epsin = round(max(Es), 3)
            avg_epsin = round(sum(Es) / len(Es), 3)

            while True:
                if "SelEpsin" in InputDict:
                    SelEpsin = InputDict["SelEpsin"]
                    print(f"SelEpsin = {SelEpsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SelEpsin = str(input(f""" 
 The internal dielectric constants range from {min_epsin} to {max_epsin}
 The average internal dielectric constant is {avg_epsin}
 
 Only one value for the internal dielectric constant can be used for all  
 hemes in these calculations.

 Would you like to use this average value for the PBSA calculations to
 compute heme-heme interactions? """))
                    print(f"SelEpsin = {SelEpsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                if SelEpsin.lower() in ["yes", "y"]:
                    epsin = avg_epsin
                    break
                elif SelEpsin.lower() in ["no", "n"]:
                    while True:
                        if "epsin" in InputDict:
                            epsin = InputDict["epsin"]
                            print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                            break
                        else:
                            try:
                                epsin = round(float(input(" What should the average internal static dielectric constant be? ")), 3)
                                print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                                break
                            except ValueError:
                                print(" Your entry needs to be a floating-poiint number.")
                    break
                else:
                    print(" Sorry, I didn't understand your response.")

        if not Es:
            while True:
                if "epsin" in InputDict:
                    epsin = InputDict["epsin"]
                    print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    try:
                        epsin = round(float(input(""" 
 An average interior static dielectric constant
 cannot be assigned for the PBSA calculations from 
 Lambda.txt because you choose to enter, instead of
 compute reorganization energies.

 What value would you like to use? """)), 3)
                        print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        break
                    except ValueError:
                        print(" Your entry needs to be a floating-poiint number.")
                break

    elif not os.path.isfile("Lambda.txt"):
        while True:
            if "epsin" in InputDict:
                epsin = InputDict["epsin"]
                print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                try:
                    epsin = round(float(input(""" 
 An average interior static dielectric constant
 cannot be assigned for the PBSA calculations from 
 Lambda.txt because the file does not exist. 

 what value would you like to use? """)), 3)
                    print(f"epsin = {epsin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    break
                except ValueError:
                    print(" Your entry needs to be a floating-poiint number.")
            break

    while True:
        if "epsout" in InputDict:
            epsout = InputDict["epsout"]
            print(f"epsout = {epsout}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            try:
                epsout = round(float(input(" What should the external static dielectric be? ")), 3)
                print(f"epsout = {epsout}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                break
            except ValueError:
                print(" Your entry needs to be a floating-poiint number.")
        break

    while True:
        if "istrng" in InputDict:
            istrng = InputDict["istrng"]
            print(f"istrng = {istrng}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            try:
                istrng = float(input(f" What ionic strength should be used in mM? "))
                print(f"istrng = {istrng}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                break
            except ValueError:
                print(" Your entry needs to be a floating-poiint number.")
        break

    while True:
        if "memb" in InputDict:
            memb = InputDict["memb"]
            print(f"memb = {memb}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            memb = input(f" Should there be an implicit slab membrane? ")
            print(f"memb = {memb}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if memb.lower() in ["yes", "y"]:
            membraneopt = 1

            while True:
                if "epsmem" in InputDict:
                    epsmem = InputDict["epsmem"]
                    print(f"epsmem = {epsmem}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    try:
                        epsmem = float(input(f" What should be the value of the membrane dielectric constant? "))
                        print(f"epsmem = {epsmem}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        break
                    except ValueError:
                        print(" Your entry needs to be a floating-poiint number.")
                break

            IPB = 1
            INP = 0
            ivalence = 0
            bcopt = 10
            eneopt = 1
            sasopt = 0
            solvopt = 1
            smoothopt = 1
            maxitn = 200
            nfocus = 1

            while True:
                if "mthick" in InputDict:
                    mthick = InputDict["mthick"]
                    print(f"mthick = {mthick}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    try:
                        mthick = float(input(f" What is the thickness of the desired membrane (Å)? "))
                        print(f"mthick = {mthick}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        break
                    except ValueError:
                        print(" Your entry needs to be a floating-poiint number.")
                break

            while True:
                if "SelPoretype" in InputDict:
                    SelPoretype = InputDict["SelPoretype"]
                    print(f"SelPoretype = {SelPoretype}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SelPoretype = input(f" Does the protein have a solvent-filled channel region that should be automatically detected? ")
                    print(f"SelPoretype = {SelPoretype}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                if SelPoretype.lower() in ["yes", "y"]:
                    poretype = 1
                    break
                elif SelPoretype.lower() in ["no", "n"]:
                    poretype = 0
                    break
                else:
                    print(" Sorry, I didn't understand your respond.")
            break
        elif memb.lower() in ["no", "n"]:
            while True:
                if "SelDelphi" in InputDict:
                    SelDelphi = InputDict["SelDelphi"]
                    print(f"SelDelphi = {SelDelphi}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    SelDelphi = input(f" Should a solution-phease Delphi-like calculation be performed? ")
                    print(f"SelDelphi = {SelDelphi}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                if SelDelphi.lower() in ["yes", "y"]:
                    membraneopt = 0
                    poretype = 0
                    mthick = 40.0
                    epsmem = epsout
                    IPB = 1
                    INP = 0
                    ivalence = 1
                    bcopt = 5
                    eneopt = 2
                    sasopt = 0
                    solvopt = 1
                    smoothopt = 2
                    maxitn = 100
                    nfocus = 1
                    break

                elif SelDelphi.lower() in ["no", "n"]:
                    membraneopt = 0
                    poretype = 0
                    mthick = 40.0
                    epsmem = epsout
                    IPB = 2
                    INP = 2
                    ivalence = 0
                    bcopt = 5
                    eneopt = 2
                    sasopt = 0
                    solvopt = 1
                    smoothopt = 1
                    maxitn = 100
                    nfocus = 2
                    break

                else:
                    print(" Sorry, I didn't understand your response.")
            break
        else:
            print(" Sorry, I didn't understand your response.")

    if os.path.isfile(f"pbsa_epsin{epsin}_epsout{epsout}.key"):
        print(f""" 
 Found pbsa_epsin{epsin}_epsout{epsout}.key to be used in PBSA calculation""")
    else:
        print(f"""
 Single point PB calculation
 &cntrl
  IPB={IPB},             ! Dielectric interface model with the level-set funciton
  INP={INP},             ! Non-polar solvation free energy method   
  ntx=1,             ! Read formatted coordinates from inpcrd 
  imin=1,            ! Single-point energy evalulation 
 /

 &pb
  pbtemp=300,        ! Temperature for salt effects in PB equation 
  ivalence={ivalence},        ! 
  istrng={istrng},      ! Ionic strength in mM for PB equation 
  epsin={epsin},       ! Solute region dielectric constant 
  epsout={epsout},       ! Solvent region dielectric constant 
  epsmem={epsmem},       ! Membrane dielectric constant
  membraneopt={membraneopt},     ! Turn off/on implicit slab membrane
  mthick={mthick},    ! Membrane thickness in Å
  mctrdz=0,          ! Membrane center in Z direction Å; 0 = centered at center of protein 
  poretype={poretype},        ! Turn off(0)/on(1) pore-searching algorithm
  radiopt=0,         ! Atomic radii from topology used; optimized radius (choice 1) for FE is missing 
  dprob=1.4,         ! Solvent probe radius for molecular surface definition  
  iprob=2.0,         ! Mobile ion probe radius used to define the Stern layer. 
  mprob=2.7,         ! Membrane lipid probe radius 
  sasopt=0,          ! Use solvent-excluded surface type for solute 
  triopt=1,          ! Use trimer arc dots to map analytical solvent excluded surface  
  arcres=0.25,       ! Resolution of dots (in Å) used to represent solvent-accessible arcs 
  maxarcdot=15000,    ! 
  smoothopt={smoothopt},       ! Use weighted harmonic average of epsin and epsout for boundary grid edges across solute/solvent dielectric boundary 
  saopt=1,           ! Compute solute surface area 
  decompopt=2,       ! sigma decomposiiton scheme for non-polar solvation 
  use_rmin=1,        ! Use rmin for van der waals radi, improves agreement with TIP3P
  sprob=0.557,       ! Compute dispersion term using solvent probe radius (in Å) for solvent accessible surface area 
  vprob=1.300,       ! Compute non-polar cavity solvation free energy using olvent probe radius (in Å) for molecular volume 
  rhow_effect=1.129, ! Effective water density for non-polar dispersion term
  use_sav=1,         ! Use molecular volume for cavity term
  maxsph=400,        ! Approximate number of dots to represent the maximum atomic solvent accessible surface
  npbopt=0,          ! Linear PB equation is solved  
  solvopt=1,         ! ICCG/PICCG iterative solver 
  accept=0.001,      ! Iteration convergence criterion   
  maxitn={maxitn},        ! Maximum number of iterations for finite difference solver 
  fillratio=1.5,     ! ratio between longest dimension of rectangular finite-difference grid and that of the solute    
  space=0.5,         ! Grid spacing for finite-difference solver 
  nfocus={nfocus},          ! Number of successive FD calculations for electrostatic focusing  
  fscale=8,          ! Ratio between coarse and fine grid spacings in electrostatic focussing 
  npbgrid=1,         ! Frequency for regenerating finite-difference grid 
  bcopt={bcopt},           ! Boundary grid potentials computed using all grid charges
  eneopt={eneopt},          ! Reaction field energy computed using dielectric boundary surface charges 
  frcopt=2,          ! reaction field forces and dielectric boundary forces computed with dielectric boundary surface polarized charges 
  scalec=0,          ! Dielectric boundary surface charges are not scaled before computing reaction field energy and forces
  cutfd=5,           ! Atom-based cutoff distance to remove short-range finite-difference interactions and to add pairwise charge-based interactions
  cutnb=0,           ! Atom-based cutoff distance for van der Waals interactions and pairwise Coulombic interactions when eneopt=2   
  !phiout=1,         ! Output spatial distribution of electrostatic potential f
  !phiform=2,        ! DX format of the electrostatic potential file for VMD
  !outlvlset=true,   ! Output total level set, used in locating interfaces between regions of differing dielectric constant
  !outmlvlset=true,  ! Output membrane level set, used in locating interfaces between regions of differing dielectric constant
  !npbverb=1,        ! Output verbosity; 1 is verbose 
  !isurfchg=1,       ! Save surface changes to a file
 /
        """, file=open(f"pbsa_epsin{epsin}_epsout{epsout}.key", 'w'))

    while True:
        if "HemeHemeIntCompChoice" in InputDict:
            HemeHemeIntCompChoice = InputDict["HemeHemeIntCompChoice"]
            print(f"HemeHemeIntCompChoice = {HemeHemeIntCompChoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            HemeHemeIntCompChoice = input(" Do you wish to run any needed computations in serial or parallel (s/p)? ")
            print(f"HemeHemeIntCompChoice = {HemeHemeIntCompChoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if HemeHemeIntCompChoice.lower() in ["serial", "s", "parallel", "p"]:
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    try:
        print("""
 ================================================
 Generating topologies for Redox Microstates
 ================================================""")
        PairedChargeAssignment.PairedChargeAssignment(ForceFieldDir, FFchoice, SelRefRedoxState, InputDict)
    except Exception as e:
        print(f"""
Detailed error in PairedChargeAssignment module:
Error type: {type(e).__name__}
Error message: {str(e)}

Please check the error messages above and logs for more details.""")
        raise  # This will show the full traceback

    idx = 0
    if os.path.isfile("SelResIndexing.txt"):
        with open("SelResIndexing.txt") as fp:
            NumHEC = len(fp.readlines())
            HEM = [0] * NumHEC
            fp.seek(0)

            Lines = fp.readlines()
            for line in Lines:
                HEM[idx] = int(line.strip().split(" ")[-3])
                idx += 1
        PairCount = list(itertools.combinations(HEM, r=2))
    else:
        sys.exit("""
 SelResIndexing.txt is missing.
 We cannot proceed without this file.""")

    M = [[0 for column in range(len(HEM))] for row in range(len(HEM))]
    N = [[0 for column in range(4)] for row in range(4)]

    print("State Energies", file=open('StateEnergies.txt', 'w'))
    print("""
 Each line in this file indicates from left-to-right:

 Fields 1 & 2: 
    The redox state (o = oxidized/r = reduced) and 
    zero-based index of the first heme in the 
    considered pair
 Fields 3 & 4: 
    The redox state (o = oxidized/r = reduced) and 
    zero-based index of the second heme in the 
    considered pair
 Field 5: 
    The total system energy (in eV) computed with 
    PBSA for the specified redox microstates. 
    If there are more than two hemes, the hemes 
    not specified on a given line are in the 
    reduced state.

 The end of the file presents the matrix of state energies
 where diagonal elements are oxidation energies and 
 off-diagonal elements are interaction energies.

 All energies in the matrix are in meV and relative to 
 the fully reduced system. Positive interaction energies 
 indicate how much the oxidation of one heme is disfavored 
 by the oxidation of the adjacent heme. \n
 """, file=open('StateEnergies.txt', 'a'))

    for i in range(len(PairCount)):
        Hi = PairCount[i][0] 
        Hj = PairCount[i][1]
        i = HEM.index(Hi) 
        j = HEM.index(Hj) 

        idxc = 0
        command = [''] * 4
        print(f"""
 ================================================
 Analyzing pair Heme-{Hi} - Heme-{Hj}...
 ================================================""")

        for k in ("o", "r"):
            for l in ("o", "r"):
                if not os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") or not os.path.isfile(f"{k}{Hi}-{l}{Hj}.rst7"):
                    if os.path.isfile(f"GeneratePairIntTopologiesForHems{Hi}-{Hj}.in"):
                        print(f"""
 Unable to find the prmtop and/or rst7 file for 
 the {k}{Hi}-{l}{Hj} micro-redox state, but found
 GeneratePairIntTopologiesForHems{Hi}-{Hj}.in.
 We will try to re-run TLEaP to generate the 
 needed files.""")
                        subprocess.run(f"tleap -s -f GeneratePairIntTopologiesForHems{Hi}-{Hj}.in > GeneratePairIntTopologiesForHems{Hi}-{Hj}.log", shell=True)
                        if not os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") or not os.path.isfile(f"{k}{Hi}-{l}{Hj}.rst7"):
                            sys.exit(f"""
 TLEaP failed. Please check 
 GeneratePairIntTopologiesForHems{Hi}-{Hj}.log""")
                    else:
                        sys.exit("""
 Something went wrong in the 
 PairedChargeAssignment module""")

        for k in ("o", "r"):
            for l in ("o", "r"):
                if os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop"):
                    if os.path.isfile(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out"):
                        print(f""" 
 Found pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out from a prior execution.
 This prior output will be used for the analysis.""")
                    else:
                        print(f""" 
 Did not find pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out from a prior execution.
 The calculation will be submitted.""")
                        if HemeHemeIntCompChoice.lower() in ["serial", "s"]:
                            print(f" Running PBSA calculation for {k}{Hi}-{l}{Hj} with epsin {epsin} and epsout {epsout}...")
                            print(f"  pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7")
                            subprocess.run(f"pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7", shell=True)
                        elif HemeHemeIntCompChoice.lower() in ["parallel", "p"]:
                            command[idxc] = f"pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7"
                            idxc += 1

        if idxc != 0:
            commandrev = list(filter(None, command))
            print("\n Submitting " + str(len(commandrev)) + " PBSA calculations in parallel")
            print(*commandrev, sep='\n')
            procs = [subprocess.Popen(i, shell=True) for i in commandrev]

            for p in procs:
                p.wait()
                print("  Finished: " + str(p))

        chk = 4
        if not os.path.isfile(f"pbsa_o{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out"):
            chk -= 1
            print(f" Something went wrong! pbsa_o{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out is missing.")
        if not os.path.isfile(f"pbsa_o{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out"):
            chk -= 1
            print(f" Something went wrong! pbsa_o{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out is missing.")
        if not os.path.isfile(f"pbsa_r{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out"):
            chk -= 1
            print(f" Something went wrong! pbsa_r{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out is missing.")
        if not os.path.isfile(f"pbsa_r{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out"):
            chk -= 1
            print(f" Something went wrong! pbsa_r{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out is missing.")

        if chk == 4:
            # print(" All four files found") 

            for k in ("o", "r"):
                for l in ("o", "r"):
                    if k == "o" and l == "o":
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                if 'Etot' in line and idx == 0:
                                    N[0][0] = float(line.strip().split()[2]) * 0.043
                                    print(k, i, l, j, round(N[0][0], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

                    if k == "o" and l == "r":
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                if 'Etot' in line and idx == 0:
                                    N[0][1] = float(line.strip().split()[2]) * 0.043
                                    print(k, i, l, j, round(N[0][1], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

                    if k == "r" and l == "o":
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                if 'Etot' in line and idx == 0:
                                    N[1][0] = float(line.strip().split()[2]) * 0.043
                                    print(k, i, l, j, round(N[1][0], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

                    if k == "r" and l == "r":
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                if 'Etot' in line and idx == 0:
                                    N[1][1] = float(line.strip().split()[2]) * 0.043
                                    print(k, i, l, j, round(N[1][1], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

        M[i][i] = round(((N[0][1] - N[1][1])) * 1000, 3)
        M[i][j] = round(((N[0][0] - N[1][0]) - (N[0][1] - N[1][1])) * 1000, 3)
        M[j][i] = round(((N[0][0] - N[0][1]) - (N[1][0] - N[1][1])) * 1000, 3)
        M[j][j] = round(((N[1][0] - N[1][1])) * 1000, 3)

    print("""
 Matrix of State Energies:
    Diagonal terms = oxidation energies
    Off-diagonal terms = interaction energies

    Energies in the matrix are in meV and relative to the 
    fully reduced system. Positive interaction energies 
    indicate how much the oxidation of one heme is 
    disfavored by the oxidation of the adjacent heme.\n """)

    print(*M, sep='\n')

    print("""
 This data and the individual state energies on which 
 it is based is saved to StateEnergies.txt.""")

    print("\n", *M, sep='\n', file=open('StateEnergies.txt', 'a'))

    print("""
 This data and the individual state energies on which
 it is based is saved to StateEnergies.txt.""")

    # Add analysis of heme cooperativity
    print("""
 ================================================
 Analyzing Heme Cooperativity
 ================================================""")

    # Get make_periodic parameter
    while True:
        if "MakePeriodic" in InputDict:
            make_periodic = InputDict["MakePeriodic"].lower() in ['true', 'yes', 'y', '1']
            print(f"MakePeriodic = {make_periodic}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            break
        else:
            choice = input(" Include step from last heme back to first? (yes/no): ").lower()
            print(f"MakePeriodic = {choice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            if choice in ['yes', 'y', 'no', 'n']:
                make_periodic = choice in ['yes', 'y']
                break
            else:
                print(" Please answer yes or no.")

    # Determine appropriate model based on user input
    while True:
        if "CooperativityModel" in InputDict:
            model = InputDict["CooperativityModel"]
            print(f"CooperativityModel = {model}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            model = input("""
 Which model should be used to analyze cooperativity?
 Options:
   'geo' (geometric pathway)
   'seq' (sequential pathway)
   'both' (analyze both pathways)
 Your choice: """).lower()
            print(f"CooperativityModel = {model}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if model in ['geo', 'seq', 'both']:
            break
        else:
            print(" Invalid choice. Please select 'geo', 'seq', or 'both'.")

    # Get energy shift parameter
    while True:
        if "EnergyShift" in InputDict:
            energy_shift = InputDict["EnergyShift"]
            print(f"EnergyShift = {energy_shift}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            break
        else:
            try:
                energy_shift = float(input(" Enter the energy shift value (in eV): "))
                print(f"EnergyShift = {energy_shift}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                break
            except ValueError:
                print(" Please enter a valid number.")

    # Get interaction scale parameter
    while True:
        if "InteractionScale" in InputDict:
            interaction_scale = InputDict["InteractionScale"]
            print(f"InteractionScale = {interaction_scale}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            break
        else:
            try:
                interaction_scale = float(input(" Enter the interaction scale factor: "))
                print(f"InteractionScale = {interaction_scale}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                break
            except ValueError:
                print(" Please enter a valid number.")

    # Get adjacent_only parameter
    while True:
        if "AdjacentOnly" in InputDict:
            adjacent_only = InputDict["AdjacentOnly"].lower() in ['true', 'yes', 'y', '1']
            print(f"AdjacentOnly = {adjacent_only}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            break
        else:
            choice = input(" Analyze adjacent pairs only? (yes/no): ").lower()
            print(f"AdjacentOnly = {choice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            if choice in ['yes', 'y', 'no', 'n']:
                adjacent_only = choice in ['yes', 'y']
                break
            else:
                print(" Please answer yes or no.")

    # Get plot option parameter
    while True:
        if "PlotOption" in InputDict:
            plot_option = InputDict["PlotOption"]
            print(f"PlotOption = {plot_option}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            break
        else:
            plot_option = input("""
 Which curves should be plotted?
 Options:
   'ox' (oxidized only)
   'red' (reduced only)
   'both' (both curves)
 Your choice: """).lower()
            print(f"PlotOption = {plot_option}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            if plot_option in ['ox', 'red', 'both']:
                break
            else:
                print(" Invalid choice. Please select 'ox', 'red', or 'both'.")

    # Process the analysis
    try:
        print("\nProcessing cooperativity analysis...")
        
        # Create path to StateEnergies.txt in EE directory
        state_energies_file = os.path.join(LaunchDir, "EE", "StateEnergies.txt")
        
        # Create output directory if it doesn't exist
        output_dir = os.path.join(LaunchDir, "EE")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        if model == 'both':
            # Process both models
            seq_result = AHC.process_matrix(
                name="BioDC", 
                filepath=state_energies_file,  # Full path to input file
                model='seq',
                energy_shift=energy_shift,
                interaction_scale=interaction_scale,
                adjacent_only=adjacent_only,
                make_periodic=make_periodic,
                output_dir=output_dir  # Full path to output directory
            )
            
            geo_result = AHC.process_matrix(
                name="BioDC",
                filepath=state_energies_file,  # Full path to input file
                model='geo',
                energy_shift=energy_shift,
                interaction_scale=interaction_scale,
                adjacent_only=adjacent_only,
                make_periodic=make_periodic,
                output_dir=output_dir  # Full path to output directory
            )
            
            if seq_result and geo_result:
                # Create comparison plot in same directory
                AHC.create_comparison_plot(seq_result, geo_result, output_dir)
                
        else:
            # Process single model
            result = AHC.process_matrix(
                name="BioDC",
                filepath=state_energies_file,  # Full path to input file
                model=model,
                energy_shift=energy_shift,
                interaction_scale=interaction_scale,
                adjacent_only=adjacent_only,
                make_periodic=make_periodic,
                output_dir=output_dir  # Full path to output directory
            )

        if model == 'both':
            print(f"""
 Analysis complete. Output files have been generated in the EE directory:
 - Redox plots: redox_plot_seq.png and redox_plot_geo.png
 - Data files: redox_data_BioDC_seq.txt and redox_data_BioDC_geo.txt
 - ΔG analysis files: DG_ind_BioDC.txt, DG_seq_BioDC.txt, and DG_geo_BioDC.txt
 - Potential progression plots: potential_progression_seq_BioDC.png and potential_progression_geo_BioDC.png""")
        else:
            print(f"""
 Analysis complete. Output files have been generated in the EE directory:
 - Redox plot: redox_plot_{model}.png
 - Data file: redox_data_BioDC_{model}.txt
 - ΔG analysis files: DG_ind_BioDC.txt and DG_{model}_BioDC.txt
 - Potential progression plot: potential_progression_{model}_BioDC.png""")
            
        if adjacent_only:
            print(f" - ΔG landscape plot: DG_landscape_BioDC_{model}.png")
            
    except Exception as e:
        print(f"\nError during cooperativity analysis: {str(e)}")
        print("Please check the input parameters and try again.")

    return None
################################################################################################################################################
