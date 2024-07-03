################################################################################################################################################
# Generic Modules
import os
import sys
import itertools
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import PairedChargeAssignment
################################################################################################################################################

def HemeHemeInt(ForceFieldDir, FFchoice, OutPrefix, SelRefRedoxState):

    if (os.path.isfile(f"RefState.prmtop") == False):
        print("""
 Generating the reference state topology where all hemes
 are in the reduced state.
 """)
        SelRefRedoxState = DefineRefState(OutPrefix)
    elif (os.path.isfile(f"RefState.prmtop") == True):
        print("""
 Found RefState.prmtop
 """)

    Es = []
    if (os.path.isfile("Lambda.txt") == True):
        with open('Lambda.txt', 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                word1 = 'Es        ='

                if (line.find(word1) != -1):
                    Es.append(float(line.strip().split()[2]))

        if(len(Es) != 0):
            min_epsin = round(min(Es), 3)
            max_epsin = round(max(Es), 3)
            avg_epsin = round((sum(Es)/len(Es)), 3) 

            while True:
                SelEpsin = str(input(f""" 
 The internal dielectric constants range from {min_epsin} to {max_epsin}
 The average internal dielectric constant is {avg_epsin}
 
 Only one value for the internal dielectric constant can be used for all  
 hemes in these calculations.

 Would you like to use this average value for the PBSA calculations to
 compute heme-heme interactions? """))
                if (SelEpsin == "YES") or (SelEpsin == "Yes") or (SelEpsin == "yes") or (SelEpsin == "Y") or (SelEpsin == "y"):
                    epsin = avg_epsin
                    break

                elif (SelEpsin == "NO") or (SelEpsin == "No") or (SelEpsin == "no") or (SelEpsin == "N") or (SelEpsin == "n"):
                    while True:
                        try:
                            epsin = round(float(input(" What should the average internal static dielectric constant be? ")), 3)
                        except ValueError:
                            print(" Your entry needs to be a floating-poiint number.")
                        else:
                            break
                        break
                    break
                else:
                    print(" Sorry, I didn't understand your response.")

        if(len(Es) == 0):
            while True:
                try:
                    epsin = round(float(input(""" 
 An average interior static dielectric constant
 cannot be assigned for the PBSA calculations from 
 Lambda.txt because you choose to enter, instead of
 compute reorganization energies.

 What value would you like to use? """)), 3)
                except ValueError:
                    print(" Your entry needs to be a floating-poiint number.")
                else:
                    break
                break

    elif (os.path.isfile("Lambda.txt") == False):
        while True:
            try:
                epsin = round(float(input(""" 
 An average interior static dielectric constant
 cannot be assigned for the PBSA calculations from 
 Lambda.txt because the file does not exist. 

 what value would you like to use? """)), 3)
            except ValueError:
                print(" Your entry needs to be a floating-poiint number.")
            else:
                break
            break

    while True:
        try:
            epsout = round(float(input(" What should the external static dielectric be? ")), 3)
        except ValueError:
            print(" Your entry needs to be a floating-poiint number.")
        else:
            break
        break

    while True:
        try:
            istrng = float(input(f" What ionic strength should be used in mM? "))
        except ValueError:
            print(" Your entry needs to be a floating-poiint number.")
        else:
            break
   
    while True:
        memb = input(f" Should there be an implicit slab membrane? ")

        if (memb == "YES") or (memb == "Yes") or (memb == "yes") or (memb == "Y") or (memb == "y"):
            membraneopt = 1

            while True:
                try:
                    epsmem = float(input(f" What should be the value of the membrane dielectric constant? "))
                except ValueError:
                    print(" Your entry needs to be a floating-poiint number.")
                else:
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
                try:
                    mthick = float(input(f" What is the thickness of the desired membrane (Å)? "))
                except ValueError:
                    print(" Your entry needs to be a floating-poiint number.")
                else:
                    break

            while True:
                SelPoretype = input(f" Does the protein have a solvent-filled channel region that should be automatically detected? ")
                if (SelPoretype == "YES") or (SelPoretype == "Yes") or (SelPoretype == "yes") or (SelPoretype == "Y") or (SelPoretype == "y"):
                    poretype = 1
                    break
                elif (SelPoretype == "NO") or (SelPoretype == "No") or (SelPoretype == "no") or (SelPoretype == "N") or (SelPoretype == "n"):
                    poretype = 0
                    break
                else:
                    print(" Sorry, I didn't understand your respond.")
            break
        elif (memb == "NO") or (memb == "No") or (memb == "no") or (memb == "N") or (memb == "n"):
            while True:
                SelDelphi = input(f" Should a solution-phease Delphi-like calculation be performed? ")

                if (SelDelphi == "YES") or (SelDelphi == "Yes") or (SelDelphi == "yes") or (SelDelphi == "Y") or (SelDelphi == "y"):
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

                elif (SelDelphi == "NO") or (SelDelphi == "No") or (SelDelphi == "n") or (SelDelphi == "N") or (SelDelphi == "n"):
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

    if (os.path.isfile(f"pbsa_epsin{epsin}_epsout{epsout}.key") == True):
        print(f""" 
 Found pbsa_epsin{epsin}_epsout{epsout}.key to be used in PBSA calculation""")
    elif (os.path.isfile(f"pbsa_epsin{epsin}_epsout{epsout}.key") == False):
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
  mthick={mthick}    ! Membrane thickness in Å
  mctrdz=0,          ! Membrane center in Z direction Å; 0 = centered at center of protein 
  poretype={poretype},        ! Turn off(0)/on(1) pore-searching algorithm
  radiopt=0,         ! Atomic radii from topology used; optimized radius (choice 1) for FE is missing 
  dprob=1.4,         ! Solvent probe radius for molecular surface definition  
  iprob=2.0,         ! Mobile ion probe radius used to define the Stern layer. 
  mprob=2.7,         ! Membrane lipid probe radius 
  sasopt=0,          ! Use solvent-excluded surface type for solute 
  triopt=1,          ! Use trimer arc dots to map analytical solvent excluded surface  
  arcres=0.25,       ! Resolution of dots (in Å) used to represent solvent-accessible arcs 
  maxarcdot=15000    ! 
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
        CompChoice = input(" Do you wish to run any needed computations in serial or parallel (s/p)? ")

        if (CompChoice == "SERIAL") or (CompChoice == "Serial") or (CompChoice == "serial") or (CompChoice == "S") or (CompChoice == "s") or (CompChoice == "PARALLEL") or (CompChoice == "Parallel") or (CompChoice == "parallel") or (CompChoice == "P") or (CompChoice == "p"):
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    print(f"""
 ================================================
 Generating topologies for Redox Microstates
 ================================================""")
    PairedChargeAssignment.PairedChargeAssignment(ForceFieldDir, FFchoice, SelRefRedoxState)

    idx = 0
    if (os.path.isfile("SelResIndexing.txt") == True):
        with open("SelResIndexing.txt") as fp:
            NumHEC = len(fp.readlines())
            HEM = [0]*NumHEC
            fp.seek(0)

            Lines = fp.readlines()
            for line in Lines:
                HEM[idx] = int(line.strip().split(" ")[-3])
                idx+=1
        PairCount = list(itertools.combinations(HEM, r=2))
    elif (os.path.isfile("SelResIndexing.txt") == False):
        sys.exit("""
 SelResIndexing.txt is missing.
 We cannot proveed witihout this file.""")

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
        command = ['']*(4)
        print(f"""
 ================================================
 Analyzing pair Heme-{Hi} - Heme-{Hj}...
 ================================================""")

        for k in ("o", "r"):
            for l in ("o", "r"):
                if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == False) or (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == False):
                    if (os.path.isfile(f"GeneratePairIntTopologiesForHems{Hi}-{Hj}.in") == True):
                        print(f"""
 Unable to find the prmtop and/or rst7 file for 
 tje {k}{Hi}-{l}{Hj} micro-redox state, but found
 GeneratePairIntTopologiesForHems{Hi}-{Hj}.in.
 We weill try to re-run TLEaP to generate the 
 needed files.""")
                        subprocess.run(f"tleap -s -f GeneratePairIntTopologiesForHems{Hi}-{Hj}.in > GeneratePairIntTopologiesForHems{Hi}-{Hj}.log", shell=True)
                        if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == False) or (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == False):
                            sys.exit(f"""
 TLEaP failed. please check 
 GeneratePairIntTopologiesForHems{Hi}-{Hj}.log""")
                    if (os.path.isfile(f"GeneratePairIntTopologiesForHems{Hi}-{Hj}.in") == False):
                        sys.exit("""
 Something went wrong in the 
 PairedChargeAssignment module""")

        for k in ("o", "r"):
            for l in ("o", "r"):
                if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == True):
                    if (os.path.isfile(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out") == True):
                        print(f""" 
 Found pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out from a prior execution.
 This prior output will be used for the analysis.""")
                    elif (os.path.isfile(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out") == False):
                        print(f""" 
 Did not find pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out from a prior execution.
 The calculation will be submitted.""")
                        if (CompChoice == "SERIAL") or (CompChoice == "Serial") or (CompChoice == "serial") or (CompChoice == "S") or (CompChoice == "s"):
                            print(f" Running PBSA calculation for {k}{Hi}-{l}{Hj} with epsin {epsin} and epsout {epsout}...")
                            print(f"  pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7")
                            subprocess.run(f"pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7", shell=True)
                        elif (CompChoice == "PARALLEL") or (CompChoice == "Parallel") or (CompChoice == "parallel") or (CompChoice == "P") or (CompChoice == "p"):
                            command[idxc] = f"pbsa -O -i pbsa_epsin{epsin}_epsout{epsout}.key -o pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out -p {k}{Hi}-{l}{Hj}.prmtop -c {k}{Hi}-{l}{Hj}.rst7"
                            idxc += 1
            
        if (idxc != 0):
            commandrev = list(filter(None, command))
            print("\n Submitting "+str(len(commandrev))+" PBSA calculations in parallel")
            print(*commandrev,sep='\n')
            procs = [ subprocess.Popen(i, shell=True) for i in commandrev ]

            for p in procs:
                p.wait()
                print("  Finished: "+str(p))

        chk = 4
        if (os.path.isfile(f"pbsa_o{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out") == False):
            chk -= 1
            print(f" Something went wrong! pbsa_o{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out is missing.""")
        if (os.path.isfile(f"pbsa_o{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out") == False):
            chk -= 1
            print(f" Something went wrong! pbsa_o{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out is missing.""")
        if (os.path.isfile(f"pbsa_r{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out") == False):
            chk -= 1
            print(f" Something went wrong! pbsa_r{Hi}-o{Hj}_epsin{epsin}_epsout{epsout}.out is missing.""")
        if (os.path.isfile(f"pbsa_r{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out") == False):
            chk -= 1
            print(f" Something went wrong! pbsa_r{Hi}-r{Hj}_epsin{epsin}_epsout{epsout}.out is missing.""")

        if (chk == 4):
            #print(" All four files found") 
                     
            for k in ("o", "r"):
                for l in ("o", "r"):

                    if (k == "o") and (l == "o"):
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                word1 = 'Etot'

                                if (line.find(word1) != -1) and (idx == 0):
                                    N[0][0] = float(line.strip().split()[2]) * 0.043 
                                    print(k, i, l, j, round(N[0][0], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

                    if (k == "o") and (l == "r"):
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                word1 = 'Etot'

                                if (line.find(word1) != -1) and (idx == 0):
                                    N[0][1] = float(line.strip().split()[2]) * 0.043 
                                    print(k, i, l, j, round(N[0][1], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

                    if (k == "r") and (l == "o"):
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                word1 = 'Etot'

                                if (line.find(word1) != -1) and (idx == 0):
                                    N[1][0] = float(line.strip().split()[2]) * 0.043 
                                    print(k, i, l, j, round(N[1][0], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1
                                    
                    if (k == "r") and (l == "r"):
                        idx = 0
                        with open(f"pbsa_{k}{Hi}-{l}{Hj}_epsin{epsin}_epsout{epsout}.out", 'r', encoding="latin-1") as fp:
                            lines = fp.readlines()
                            for line in lines:
                                word1 = 'Etot'

                                if (line.find(word1) != -1) and (idx == 0):
                                    N[1][1] = float(line.strip().split()[2]) * 0.043 
                                    print(k, i, l, j, round(N[1][1], 3), file=open('StateEnergies.txt', 'a'))
                                    idx += 1

        M[i][i] = round(((N[0][1] - N[1][1]))*1000, 3)
        M[i][j] = round(((N[0][0] - N[1][0]) - (N[0][1] - N[1][1]))*1000, 3)
        M[j][i] = round(((N[0][0] - N[0][1]) - (N[1][0] - N[1][1]))*1000, 3)
        M[j][j] = round(((N[1][0] - N[1][1]))*1000, 3)
    print("""
 Matrix of State Energies:
    Diagonal terms = oxidation energies
    Off-diagonal terms = interaction energies

    Energies in the matirx are in meV and relative to the 
    fully reduced system. Positive interaction energies 
    indicate how much the oxidation of one heme is 
    disfavored by the oxidation of the adjacent heme.\n """)

    print(*M, sep='\n')

    print("""
 This data and the individual state energies on whcih 
 it is based is saved to StateEnergies.txt.""")

    print("\n", *M, sep='\n', file=open('StateEnergies.txt', 'a'))

################################################################################################################################################
