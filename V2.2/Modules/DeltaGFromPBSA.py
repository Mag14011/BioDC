# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
# Custom Modules 
import DefineRefState
import GenerateRedoxStateTopologies

def DeltaGFromPBSA(LaunchDir, ForceFieldDir, OutPrefix, SolvEnv, PolySel, InputDict):
    SelRefRedoxState, FFchoice = DefineRefState.DefineRefState(LaunchDir, ForceFieldDir, OutPrefix, InputDict)

    print(f"""
 Starting from the reference state, oxidized- and 
 reduced-state topologies for each heme will now 
 be generated""")
    GenerateRedoxStateTopologies.GenerateRedoxStateTopologies(ForceFieldDir, FFchoice, SelRefRedoxState)

    while True:
        if "DeltaGFromPBSACompChoice" in InputDict:
            DeltaGFromPBSACompChoice = InputDict["DeltaGFromPBSACompChoice"]
            print(f"DeltaGFromPBSACompChoice = {DeltaGFromPBSACompChoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        else:
            DeltaGFromPBSACompChoice = input("\n Do you wish to run any needed computations in serial or parallel (s/p)? ")
            print(f"DeltaGFromPBSACompChoice = {DeltaGFromPBSACompChoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        if DeltaGFromPBSACompChoice.lower() in ["serial", "s", "parallel", "p"]:
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    Es = []
    if os.path.isfile(f"Lambda.txt"):
        with open('Lambda.txt', 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                if 'Es        =' in line:
                    Es.append(float(line.strip().split()[2]))
    else:
        print(""" 
 Interior static dielectric constants cannot be assigned
 from Lambda.txt for the PBSA calculations because the file
 was not found. You will be asked to specify this parameter.
 """)

    idx = 0; idxc = 0; 
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEM         = [0]*x
        SelEpsin    = [0]*x
        epsin       = [0]*x
        epsout      = [0]*x
        istrng      = [0]*x
        memb        = [0]*x
        epsmem      = [0]*x
        mthick      = [0]*x
        SelPoretype = [0]*x
        poretype    = [0]*x
        SelDelphi   = [0]*x
        RedoxState  = [0]*2
        command     = [' ']*(2*x)

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEM[idx] = int(line.strip().split(" ")[1])

            if len(Es) != 0:
                if idx == 0:
                    epsin[idx] = round(Es[idx], 3) 
                    if f"SelEpsin{idx}" in InputDict:
                        SelEpsin[idx] = InputDict[f"SelEpsin{idx}"]
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        SelEpsin[idx] = input(f""" 
 The internal static dielectric constant for heme-{HEM[idx]} is {epsin[idx]}.
 Should this value be used? """)
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if SelEpsin[idx].lower() in ["no", "n"]:
                        while True:
                            try:
                                epsin[idx] = round(float(input(f" What should the internal dielectric constant be for the PBSA calculation on heme-{HEM[idx]}? ")), 3)
                                print(f"epsin{idx} = {epsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                            except ValueError:
                                print(" Your entry needs to be a floating-poiint number.")
                            else:
                                break
                elif idx != 0 and idx != len(HEM) - 1:
                    epsin[idx] = round(((Es[idx-1] + Es[idx]) / 2), 3)
                    if f"SelEpsin{idx}" in InputDict:
                        SelEpsin[idx] = InputDict[f"SelEpsin{idx}"]
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        SelEpsin[idx] = input(f""" 
 The internal static dielectric constant for heme-{HEM[idx]} is {epsin[idx]}.
   This value is an average of the static dielectric constants estimated for 
   the (i-1, i) and (i, i+1) heme pairs, where i = heme-{HEM[idx]}.
 Should this value be used? """)
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if SelEpsin[idx].lower() in ["no", "n"]:
                        while True:
                            try:
                                epsin[idx] = round(float(input(f" What should the internal dielectric constant be for the PBSA calculation on heme-{HEM[idx]}? ")), 3)
                                print(f"epsin{idx} = {epsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                            except ValueError:
                                print(" Your entry needs to be a floating-poiint number.")
                            else:
                                break
                elif idx == len(HEM) - 1:
                    if PolySel.lower() in ["no", "n"]:
                        epsin[idx] = round(Es[len(Es)-1], 3)
                    if PolySel.lower() in ["yes", "y"]:
                        epsin[idx] = round(((Es[len(Es)-1] + Es[0]) / 2), 3)
                    if f"SelEpsin{idx}" in InputDict:
                        SelEpsin[idx] = InputDict[f"SelEpsin{idx}"]
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        SelEpsin[idx] = input(f""" 
 The internal static dielectric constant for heme-{HEM[idx]} is {epsin[idx]}.
   This value is an average of the static dielectric constants estimated for 
   the (i-1, i) and (i, i+1) heme pairs, where i = heme-{HEM[idx]}.

   If you inidcated the structure is polymeric and you included the first
   heme in the next sub-unit in the linear sequence to be considered,
   i-1 is the penultimate heme in one sub-unit and i+1 is the second heme 
   in the next sub-unit.

   If you inidcated the structure is polymeric but did NOT include the 
   first heme of the next sub-unit in the linear sequence to be considered,
   the static dielectric constant is the average of the values for the last
   pair in one sub-unit and the first pair in the next sub-unit.

   If your structure is not polymeric, the static dielectric constant is the
   value estimated for the last heme pair in the specified linear sequence.
   
 Should this value be used? """)
                        print(f"SelEpsin{idx} = {SelEpsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if SelEpsin[idx].lower() in ["no", "n"]:
                        while True:
                            try:
                                epsin[idx] = round(float(input(f" What should the internal dielectric constant be for the PBSA calculation on heme-{HEM[idx]}? ")), 3)
                                print(f"epsin{idx} = {epsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                            except ValueError:
                                print(" Your entry needs to be a floating-poiint number.")
                            else:
                                break
            else:
                if f"epsin{idx}" in InputDict:
                    epsin[idx] = round(float(InputDict[f"epsin{idx}"]), 3)
                    print(f"epsin{idx} = {epsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    while True:
                        try:
                            epsin[idx] = round(float(input(f" What should the internal dielectric constant be for the PBSA calculation on heme-{HEM[idx]}? ")), 3)
                            print(f"epsin{idx} = {epsin[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        except ValueError:
                            print(" Your entry needs to be a floating-poiint number.")
                        else:
                            break

            if f"epsout{idx}" in InputDict:
                epsout[idx] = round(float(InputDict[f"epsout{idx}"]), 3)
                print(f"epsout{idx} = {epsout[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                while True:
                    try:
                        epsout[idx] = round(float(input(f" What should the external dielectric constant be for heme-{HEM[idx]}? ")), 3)
                        print(f"epsout{idx} = {epsout[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    except ValueError:
                        print(" Your entry needs to be a floating-poiint number.")
                    else:
                        break

            if f"istrng{idx}" in InputDict:
                istrng[idx] = float(InputDict[f"istrng{idx}"])
                print(f"istrng{idx} = {istrng[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                while True:
                    try:
                        istrng[idx] = float(input(f" What ionic strength should be used in mM? "))
                        print(f"istrng{idx} = {istrng[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    except ValueError:
                        print(" Your entry needs to be a floating-poiint number.")
                    else:
                        break

            if f"memb{idx}" in InputDict:
                memb[idx] = InputDict[f"memb{idx}"]
                print(f"memb{idx} = {memb[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
            else:
                memb[idx] = input(f" Should there be an implicit slab membrane? ")
                print(f"memb{idx} = {memb[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if memb[idx].lower() in ["yes", "y"]:
                membraneopt = 1
                if f"epsmem{idx}" in InputDict:
                    epsmem[idx] = float(InputDict[f"epsmem{idx}"])
                    print(f"epsmem{idx} = {epsmem[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    while True:
                        try:
                            epsmem[idx] = float(input(f" What should be the value of the membrane dielectric constant? "))
                            print(f"epsmem{idx} = {epsmem[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        except ValueError:
                            print(" Your entry needs to be a floating-point number.")
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

                if f"mthick{idx}" in InputDict:
                    mthick[idx] = float(InputDict[f"mthick{idx}"])
                    print(f"mthick{idx} = {mthick[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    while True:
                        try:
                            mthick[idx] = float(input(f" What is the thickness of the desired membrane (Å)? "))
                            print(f"mthick{idx} = {mthick[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        except ValueError:
                            print(" Your entry needs to be a floating-point number.")
                        else:
                            break

                if f"SelPoretype{idx}" in InputDict:
                    SelPoretype[idx] = InputDict[f"SelPoretype{idx}"]
                    print(f"SelPoretype{idx} = {SelPoretype[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                else:
                    while True:
                        SelPoretype[idx] = input(f" Does the protein have a solvent-filled channel region that should be automatically detected (yes/no)? ")
                        print(f"SelPoretype{idx} = {SelPoretype[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        if SelPoretype[idx].lower() in ["yes", "y"]:
                            poretype[idx] = 1
                            break
                        elif SelPoretype[idx].lower() in ["no", "n"]:
                            poretype[idx] = 0
                            break
                        else:
                            print(" Sorry, I didn't understand your response.")
            else:
                while True:
                    if f"SelDelphi{idx}" in InputDict:
                        SelDelphi[idx] = InputDict[f"SelDelphi{idx}"]
                        print(f"SelDelphi{idx} = {SelDelphi[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    else:
                        SelDelphi[idx] = input(f" Should a solution-phase Delphi-like calculation be performed? ")
                        print(f"SelDelphi{idx} = {SelDelphi[idx]}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    if SelDelphi[idx].lower() in ["yes", "y"]:
                        membraneopt = 0
                        poretype[idx] = 0
                        mthick[idx] = 40.0
                        epsmem[idx] = epsout[idx]
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

                    elif SelDelphi[idx].lower() in ["no", "n"]:
                        membraneopt = 0
                        poretype[idx] = 0
                        mthick[idx] = 40.0
                        epsmem[idx] = epsout[idx]
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

            if os.path.isfile(f"pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key"):
                print(f""" 
 Found pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key to be used in PBSA calculation for heme-{HEM[idx]}""")
            else:
                print(f"""
 Single point PB calculation
 &cntrl
  IPB={IPB},             ! Dielectric interface model with the level-set function
  INP={INP},             ! Non-polar solvation free energy method   
  ntx=1,             ! Read formatted coordinates from inpcrd 
  imin=1,            ! Single-point energy evaluation 
 /

 &pb
  pbtemp=300,        ! Temperature for salt effects in PB equation 
  ivalence={ivalence},        ! 
  istrng={istrng[idx]},      ! Ionic strength in mM for PB equation 
  epsin={epsin[idx]},       ! Solute region dielectric constant 
  epsout={epsout[idx]},       ! Solvent region dielectric constant 
  epsmem={epsmem[idx]},       ! Membrane dielectric constant
  membraneopt={membraneopt},     ! Turn off/on implicit slab membrane
  mthick={mthick[idx]}    ! Membrane thickness in Å
  mctrdz=0,          ! Membrane center in Z direction Å; 0 = centered at center of protein 
  poretype={poretype[idx]},        ! Turn off(0)/on(1) pore-searching algorithm
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
  decompopt=2,       ! sigma decomposition scheme for non-polar solvation 
  use_rmin=1,        ! Use rmin for van der waals radii, improves agreement with TIP3P
  sprob=0.557,       ! Compute dispersion term using solvent probe radius (in Å) for solvent accessible surface area 
  vprob=1.300,       ! Compute non-polar cavity solvation free energy using solvent probe radius (in Å) for molecular volume 
  rhow_effect=1.129, ! Effective water density for non-polar dispersion term
  use_sav=1,         ! Use molecular volume for cavity term
  maxsph=400,        ! Approximate number of dots to represent the maximum atomic solvent accessible surface
  npbopt=0,          ! Linear PB equation is solved  
  solvopt=1,         ! ICCG/PICCG iterative solver 
  accept=0.001,      ! Iteration convergence criterion   
  maxitn={maxitn},        ! Maximum number of iterations for finite difference solver 
  fillratio=1.5,     ! Ratio between longest dimension of rectangular finite-difference grid and that of the solute    
  space=0.5,         ! Grid spacing for finite-difference solver 
  nfocus={nfocus},          ! Number of successive FD calculations for electrostatic focusing  
  fscale=8,          ! Ratio between coarse and fine grid spacings in electrostatic focusing 
  npbgrid=1,         ! Frequency for regenerating finite-difference grid 
  bcopt={bcopt},           ! Boundary grid potentials computed using all grid charges
  eneopt={eneopt},          ! Reaction field energy computed using dielectric boundary surface charges 
  frcopt=2,          ! Reaction field forces and dielectric boundary forces computed with dielectric boundary surface polarized charges 
  scalec=0,          ! Dielectric boundary surface charges are not scaled before computing reaction field energy and forces
  cutfd=5,           ! Atom-based cutoff distance to remove short-range finite-difference interactions and to add pairwise charge-based interactions
  cutnb=0,           ! Atom-based cutoff distance for van der Waals interactions and pairwise Coulombic interactions when eneopt=2   
  !phiout=1,         ! Output spatial distribution of electrostatic potential 
  !phiform=2,        ! DX format of the electrostatic potential file for VMD
  !outlvlset=true,   ! Output total level set, used in locating interfaces between regions of differing dielectric constant
  !outmlvlset=true,  ! Output membrane level set, used in locating interfaces between regions of differing dielectric constant
  !npbverb=1,        ! Output verbosity; 1 is verbose 
  !isurfchg=1,       ! Save surface changes to a file
 /
                """, file=open(f"pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key", 'w'))

            for state in ["ox", "red"]:
                if state == "ox":
                    RedoxState[0] = "o"
                if state == "red":
                    RedoxState[0] = "r"

                if os.path.isfile(f"pbsa_{RedoxState[0]}{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out"):
                    print(f""" 
 Found pbsa_{RedoxState[0]}{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out from a prior execution.
 This prior output will be used for the analysis.""")
                else: 
                    if DeltaGFromPBSACompChoice.lower() in ["serial", "s"]:
                        print(f" Running PBSA calculation for {state}. Heme-{HEM[idx]} ...")
                        print(f"  pbsa -O -i pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key -o pbsa_{RedoxState[0]}{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out -p {RedoxState[0]}{HEM[idx]}.prmtop -c {RedoxState[0]}{HEM[idx]}.rst7")
                        subprocess.run(f"pbsa -O -i pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key -o pbsa_{RedoxState[0]}{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out -p {RedoxState[0]}{HEM[idx]}.prmtop -c {RedoxState[0]}{HEM[idx]}.rst7", shell=True)
                    if DeltaGFromPBSACompChoice.lower() in ["parallel", "p"]:
                        command[idxc] = f"pbsa -O -i pbsa-{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.key -o pbsa_{RedoxState[0]}{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out -p {RedoxState[0]}{HEM[idx]}.prmtop -c {RedoxState[0]}{HEM[idx]}.rst7"
                        idxc += 1

            idx += 1

        if idxc != 0:
            commandrev = list(filter(None, command))
            print("\n Submitting "+str(len(commandrev))+" PBSA calculations in parallel")
            print(*commandrev, sep='\n')
            procs = [subprocess.Popen(i, shell=True) for i in commandrev]

            for p in procs:
                p.wait()
                print("  Finished: "+str(p))

    DEtot = [0]*(len(HEM))
    DEelec = [0]*(len(HEM))
    DG = [0]*(len(HEM)-1)
    for idx in range(len(HEM)):

        idx1 = 0; idx2 = 0;
        word1 = 'Etot'; word2 = 'EELEC'
        try:
            with open(f"pbsa_o{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out", 'r', encoding="utf-8") as fp:
                lines = fp.readlines()

                for line in lines:
                    if (line.find(word1) != -1) and (idx1 == 0):
                        EtotOx = float(line.strip().split()[2]) * 0.043 
                        idx1 += 1

                    if (line.find(word2) != -1) and (idx2 == 0):
                        EelecOx = float(line.strip().split()[2]) * 0.043
                        idx2 += 1
        except UnicodeDecodeError:
            with open(f"pbsa_o{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out", 'r', encoding="latin-1") as fp:
                lines = fp.readlines()

                for line in lines:
                    if (line.find(word1) != -1) and (idx1 == 0):
                        EtotOx = float(line.strip().split()[2]) * 0.043 
                        idx1 += 1

                    if (line.find(word2) != -1) and (idx2 == 0):
                        EelecOx = float(line.strip().split()[2]) * 0.043
                        idx2 += 1

        idx3 = 0; idx4 = 0;
        try:
            with open(f"pbsa_r{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out", 'r', encoding="utf-8") as fp:
                lines = fp.readlines()

                for line in lines:
                    if (line.find(word1) != -1) and (idx3 == 0):
                        EtotRed = float(line.strip().split()[2]) * 0.043
                        idx3 += 1

                    if (line.find(word2) != -1) and (idx4 == 0):
                        EelecRed = float(line.strip().split()[2]) * 0.043
                        idx4 += 1

        except UnicodeDecodeError:
            with open(f"pbsa_r{HEM[idx]}_epsin{epsin[idx]}_epsout{epsout[idx]}.out", 'r', encoding="latin-1") as fp:
                lines = fp.readlines()

                for line in lines:
                    if (line.find(word1) != -1) and (idx3 == 0):
                        EtotRed = float(line.strip().split()[2]) * 0.043
                        idx3 += 1

                    if (line.find(word2) != -1) and (idx4 == 0):
                        EelecRed = float(line.strip().split()[2]) * 0.043
                        idx4 += 1

        DEtot[idx] = (EtotOx - EtotRed) 
        DEelec[idx] = (EelecOx - EelecRed) 

        if idx == 0:
            print("\n Result:")

        print(f"""  step={idx} HEM-{HEM[idx]} EtotOx={EtotOx:.3f} eV EtotRed={EtotRed:.3f} eV DEtot={DEtot[idx]:.3f} eV EelecOx={EelecOx:.3f} eV EelecRed={EelecRed:.3f} eV DEelec={DEelec[idx]:.3f} eV""")

    for idx in range(len(HEM)-1):
        DG[idx] = -1 * ((-1 * DEtot[idx]) + (DEtot[idx+1])) 

        if idx == 0:
            print(f"(HEM-{HEM[idx]} = {DEtot[idx]:.3f} eV) -> (HEM-{HEM[idx+1]} = {DEtot[idx+1]:.3f} eV); DG = {DG[idx]:.3f} eV", file=open('DG.txt', 'w'))
        else:
            print(f"(HEM-{HEM[idx]} = {DEtot[idx]:.3f} eV) -> (HEM-{HEM[idx+1]} = {DEtot[idx+1]:.3f} eV); DG = {DG[idx]:.3f} eV", file=open('DG.txt', 'a'))

    return DG, SelRefRedoxState, FFchoice
