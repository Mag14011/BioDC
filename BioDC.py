import os
import sys
import subprocess
from subprocess import Popen
import math 
import numpy as np
import derrida
import blumberger     
from pandas import read_csv

################################################################################################################################################
# Written by Matthew J. Guberman-Pfeffer on 05/23 - 06/02/2023; latest revision: 08/30/2023 
################################################################################################################################################

################################################################################################################################################

def CheckProgInPath():

    while True:
        ProgInPath = input("""
 This program requires VMD and the AmberTools package to be 
 in your system's PATH variable. Are they? (yes/no)? """)

        if (ProgInPath == 'Yes') or (ProgInPath == "yes") or (ProgInPath == "Y") or (ProgInPath == "y"):
            print(""" 
 Good! Now, here are the PDBs in the present direcotry:\n""")
            return ProgInPath
            break
        elif (ProgInPath == "No") or (ProgInPath == "no") or (ProgInPath == "N") or (ProgInPath == "n"):
            sys.exit("""
 Please make VMD findalbe in your system PATH variable 
 and then re-run this program \n""")
        else: 
            print(" Sorry, I didn't understand your response.")

################################################################################################################################################

def InitialStructureSelection():
    for x in os.listdir():
        if x.endswith(".pdb"):
            print(x)

    while True:
        OriginalPDB = input("""
 Which PDB would you like to setup 
 (omit the .pdb file extension)? """)

        if (os.path.isfile(OriginalPDB + ".pdb") == True):
            print(" That PDB was found! ")
            return OriginalPDB
            break
        else: 
            print(" That PDB does not exist unfortunately")

################################################################################################################################################

def CreateResIndexing(PDB):
    if (os.path.isfile("CorrectedResIndexing.txt") == True):
        print(""" 
 Found CorrectedResIndexing.txt! 
 Copying CorrectedResIndexing.txt to ResIndexing.txt.""")
        subprocess.run("cp CorrectedResIndexing.txt ResIndexing.txt", shell=True)
    else:
        while True:
            CreateResIndexingMethod = input("""
 We need to create a file (ResIndexing.txt) that identifies 
 the IDs of the Cys and His residues bound to the heme macrocycle.
 Would you like to create it automatically or manually (auto/man)? """)

            if (CreateResIndexingMethod == "auto"):
                DistThresh = input("""
 The automated creation of ResIndexing.txt has two parts:
   1) Writting CreateResIndexing.tcl 
   2) Submitting the TCL script to Visual Molecular Dynamics (VMD)

 The TCL script identifies the Cys and His residues within a 
 distance cutoff of each heme group and assumes that the residues 
 found within that cutoff are bonded to that heme.

 What distance threshold would you liek to use? 
 (recommended = 2.5 angstroms)                  """)

                print("""
 Please make sure that the correct residues are identified by, 
 for example, creating representations in VMD with the residue IDs 
 givenon each line of the ResIndexing.txt file. 

 If the wrong residues are identified, the setup later with TLEaP 
 will fail because the bond definitions will be wrong. In this case, 
 please correct the residue IDs and save the changes to 
 CorrectedResIndexing.txt. When you re-run this python scirpt, the
 CorrectedResIndexing.txt file will be detected and used to replace
 ResIndexing.txt. 
                """)

                print("""
 set DistThresh {}
 mol new {}.pdb""".format(DistThresh, PDB), file=open('CreateResIndexing.tcl', 'w'))

                print("""
 set HEC [[atomselect top "resname HEC HEM and name FE"] get resid]

 set out [open "ResIndexing.txt" w]

 set i {0}
 foreach HResID $HEC {
    set ShiftHResID [expr {$HResID + 1000}]

    set heme [atomselect top "resname HEC HEM and resid $HResID"]
    $heme set resid $ShiftHResID

    set NewHResID [expr {$HResID + ($i * 2)}]
    set HIS [lsort -integer [[atomselect top "resname HIS and name NE2 and within $DistThresh of resname HEC HEM and resid $ShiftHResID"] get resid]]
    set CYS [lsort -integer [[atomselect top "resname CYS and name SG  and within $DistThresh of resname HEC HEM and resid $ShiftHResID"] get resid]]
    set CYSb [lindex $CYS 0]
    set CYSc [lindex $CYS 1]
    set HISp [expr {$CYSc + 1}]

    if {[lindex $HIS 0] != $HISp} {
        set HISd [lindex $HIS 0]
    } else {
        set HISd [lindex $HIS 1]    }

    puts      "$CYSb $CYSc $HISp $HISd $ShiftHResID $NewHResID"
    puts $out "$CYSb $CYSc $HISp $HISd $ShiftHResID $NewHResID"
    set i [expr {$i+1}]
 }
 close $out

 set all [atomselect top all] """, file=open('CreateResIndexing.tcl', 'a'))

                print("""
 $all writepdb {}_renumd.pdb""".format(PDB), file=open('CreateResIndexing.tcl', 'a')
                )

                print("""
 exit
                """, file=open('CreateResIndexing.tcl', 'a'))
        
                subprocess.run("vmd -e CreateResIndexing.tcl > CreateResIndexing.log", shell=True)
                break
            elif (CreateResIndexingMethod == "man") or (CreateResIndexingMethod == "manual"): 
                print("""
 To create ResIndexing.txt by hand:
    Create a txt file with an editor of your choosing (e.g. 
    vi ResIndexing.txt). In this file,there must be one line for 
    each heme cofactor in your structure. Each line should contain 
    of six numbers separated by a single space. The first four 
    numbers from left-to-right should be the residue IDs for the
        (1) Cys attached to the B-ring of the heme macrocycle,
        (2) Cys attached to the C-ring of the heme macrocycle,
        (3) Proximal His ligated to the Fe center of the heme,
        (4) Distal His ligated to the Fe center of the heme,

    The fifth and sixth numbers on a line in ResIndexing.txt pertain 
    to the heme group itself.

    The fifth number is the current residue ID for the heme, shifted 
    by +1000. A shift in the residue ID is needed to (potentially) 
    avoid assigning the propionates, whcih are made into separate 
    residues to whatever came immediately after the heme in the 
    original PDB. (See the below example.)

    The sixth number is the residue ID the heme will be assigned in 
    the final and properly formatted PDB for use with TLEap. This 
    number should be the original residue ID for the heme + (i * 2), 
    where i is a zero-based index that counts the number of hemes 
    in your system.

    To given an example:
        Let's say you have a di-heme system with the residue IDs of 
        the hemes being 114 and 115. The entries in ResIndexing.txt 
        should be:

            CysB1 CysC1 HisP1 HisD1 1114 114
            CysB2 CysC2 HisP2 HisD2 1115 117

        where CysB, CysC, HisP, and HisD are respectively defined by 
        points 1-4 above, and the 1 and 2 are used to distinguish the 
        Cys/His for the first and second heme, respectively.

    Note that the renumbering allows the propionates of the first 
    heme to be assigned to residues 115 and 116 without overlapping 
    with the second heme, which originally was residue 115.
                    """)
                if (os.path.isfile("ResIndexing.txt") == True):
                    print(" Found ResIndexing.txt! Moving to the next step...")
                    break
                else: 
                    print(" ResIndexing.txt not found!")
                    sys.exit(""" 
 Please re-run this scirpt once you have created ResIndexing.txt and 
 re-select the manual option. \n""")
            else:
                print("""
 Sorry, I didn't understand your response.
 Let's try again.""")

################################################################################################################################################

def ProcessPDB(PDB):
    print("""
 Your ResIndexing file will now be used with VMD
 to create the HEH, HIO, and CYO residues.
 
 We will write and submit a script called SetupStructure.tcl to 
 perform this magic.

 CAUTION: The magic of the script is only as good as the information 
 in ResIndexing.txt. If the wrong residue IDs are specified, 
 everything from hereon out will be, put politely, junk!
    """, end=" ")

    while True:
        ProcessPDBChoice = input(""" 
 Shall we venture forward with SetupStructure.tcl (yes/no)? """)

        if (ProcessPDBChoice == "yes") or (ProcessPDBChoice == "Yes") or (ProcessPDBChoice == "y") or (ProcessPDBChoice == "Y"):
            print("""
 #-------------------------------------------------------------
 #Input
 mol load pdb {}_renumd.pdb
 #-------------------------------------------------------------""".format(PDB), file=open('SetupStructure.tcl', 'w'))

            print("""
 set INPUT   [open  "ResIndexing.txt" r]
 set CONTENT [read -nonewline $INPUT]
 close $INPUT

 set RECORDS [split $CONTENT "\\n"]
 foreach REC $RECORDS {
    set FIELDS [split $REC " "]
    lassign $FIELDS CYS1b CYS1c HIS1p HIS1d HEH1 idx
    puts "$CYS1b $CYS1c $HIS1p $HIS1d $HEH1 $idx "
    set chk1 [[atomselect top "alpha and resid $HIS1p $HIS1d"] get {resname resid}]
    set chk2 [[atomselect top "alpha and resid $CYS1b $CYS1c"] get {resname resid}]
    puts "Selected ligands: $chk1 $chk2"

    #Rename heme to HEH
    set  HEH [atomselect top "resid $HEH1 and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D"]
    $HEH set resname HEH

    #Reindex HEH ResID
    set H1  [atomselect top "resname HEH and resid $HEH1  "]
    $H1 set resid $idx
 #-------------------------------------------------------------
    #Relabel Propinates
    set PRN1A  [atomselect top "resid $HEH1  and name CAA CBA CGA O1A O2A"]
    set PRN1D  [atomselect top "resid $HEH1  and name CAD CBD CGD O1D O2D"]
    $PRN1A set resname PRN; $PRN1D set resname PRN
    $PRN1A set resid [expr {$idx + 1}];  $PRN1D set resid [expr {$idx + 2}]
    set CA    [atomselect top "(resid [expr {$idx + 1}] and name CAA) or (resid [expr {$idx + 2}] and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resid [expr {$idx + 1}] and name CBA) or (resid [expr {$idx + 2}] and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resid [expr {$idx + 1}] and name CGA) or (resid [expr {$idx + 2}] and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resid [expr {$idx + 1}] and name O1A) or (resid [expr {$idx + 2}] and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resid [expr {$idx + 1}] and name O2A) or (resid [expr {$idx + 2}] and name O2D)"]
    $O2 set name O2

 #-------------------------------------------------------------
    #Rename His Residues

    set HISs [atomselect top "resid $HIS1p $HIS1d and name CB CG CD2 ND1 CE1 NE2"]
    $HISs set resname HEH
    set HISb [atomselect top "resid $HIS1p $HIS1d and name N CA C O"]
    $HISb set resname HIO
 #-------------------------------------------------------------
    #Rename Proximal Histidine

    set CB1  [atomselect top "name CB  and resid $HIS1p"]
    $CB1 set name CB1
    set CG1  [atomselect top "name CG  and resid $HIS1p"]
    $CG1 set name CG1
    set ND11 [atomselect top "name ND1 and resid $HIS1p"]
    $ND11 set name ND11
    set CD21 [atomselect top "name CD2 and resid $HIS1p"]
    $CD21 set name CD21
    set CE11 [atomselect top "name CE1 and resid $HIS1p"]
    $CE11 set name CE11
    set NE21 [atomselect top "name NE2 and resid $HIS1p"]
    $NE21 set name NE21
 #-------------------------------------------------------------
    #Rename Distal Histidine

    set CB2  [atomselect top "name CB and resid $HIS1d"]
    $CB2 set name CB2
    set CG2  [atomselect top "name CG and resid $HIS1d"]
    $CG2 set name CG2
    set ND12 [atomselect top "name ND1 and resid $HIS1d"]
    $ND12 set name ND12
    set CD22 [atomselect top "name CD2 and resid $HIS1d"]
    $CD22 set name CD22
    set CE12 [atomselect top "name CE1 and resid $HIS1d"]
    $CE12 set name CE12
    set NE22 [atomselect top "name NE2 and resid $HIS1d"]
    $NE22 set name NE22
 #-------------------------------------------------------------

    set  HH1   [atomselect top "resid $HIS1p $HIS1d and resname HEH"]
    $HH1 set   resid $idx
 #-------------------------------------------------------------
    #Rename Cys Residues

    set CYSS [atomselect top "resid $CYS1b $CYS1c and name CB SG"]
    $CYSS set resname HEH
    set CYSB [atomselect top "resid $CYS1b $CYS1c and name N CA C O"]
    $CYSB set resname CYO
 #-------------------------------------------------------------
    #Rename Cysteine Attached to Ring B

    set CBB2 [atomselect top "name CB and resid $CYS1b"]
    $CBB2 set name CBB2
    set SGB2 [atomselect top "name SG and resid $CYS1b"]
    $SGB2 set name SGB2
#-------------------------------------------------------------
    #Rename Cysteine Attached to Ring C

    set CBC1 [atomselect top "name CB and resid $CYS1c"]
    $CBC1 set name CBC1
    set SGC1 [atomselect top "name SG and resid $CYS1c"]
    $SGC1 set name SGC1
 #-------------------------------------------------------------

    set CH1 [atomselect top "resname HEH and resid $CYS1b $CYS1c"]
    $CH1 set resid $idx
 #-------------------------------------------------------------

    set HEM1 [atomselect top "resname HEH and resid $idx"]
    $HEM1 num
    set PRN1A  [atomselect top "resname PRN and resid [expr {$idx + 1}]"]
    set PRN1D  [atomselect top "resname PRN and resid [expr {$idx + 2}]"]
    $HEM1  writepdb HEM$idx.pdb; $PRN1A  writepdb PRNA$idx.pdb; $PRN1D  writepdb PRND$idx.pdb}
    set prot [atomselect top "noh and all and not resname HEH PRN"]
    $prot writepdb prot.pdb
 #-------------------------------------------------------------

    exit
            """, file=open('SetupStructure.tcl', 'a'))
            subprocess.run("vmd -e SetupStructure.tcl > SetupStructure.log", shell=True)
            print("""
 VMD finished. Please check SetupStructure.log for any erros. You may 
 also want to inspect the generated PDBs for the protein,each heme, 
 and each heme propionic acid group.
            """, end=" ")
            break

        elif (ProcessPDBChoice == "no") or (ProcessPDBChoice == "No") or (ProcessPDBChoice == "n") or (ProcessPDBChoice == "N"):
            sys.exit("""
 I'm sorry but I don't know how to proceed without running 
 SetupStructure.tcl using VMD. Hopefully this program was helpful 
 up to this point! If you have another way to format the PDB and 
 generate the TLEaP input file, that's great! If not, please don't 
 hestitate to re-run this program. I promise it'll save you many 
 headaches and much time! 
            """)
            break
        else:
            print("""
 I'm sorry but I didn't understand your selection.
 It's just a bit of artificial stupidity!
 Let's try agin.
            """, end=" ")
################################################################################################################################################

def ReBuildStructure(PDB):
    print("""
 Now, we need to stitch the edited PDBs of the protein, hemes, 
 and heme propionic acid groups into a single PDB. Then, this 
 re-constructed PDB of the multi-heme protein will be processed
 with TLEaP of the AmberTools package to generate topology and 
 coordinate files.
    """, end=" ")
    
    OutPrefix = input(" \n Prefix for output parm/rst7 ")

    idx = 0
    with open("ResIndexing.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x
        HISp = [0]*x
        HISd = [0]*x
        CYSb = [0]*x
        CYSc = [0]*x

        fp.seek(0)
        Lines = fp.readlines()

        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[5])
            HISp[idx] = int(line.strip().split(" ")[2])
            HISd[idx] = int(line.strip().split(" ")[3])
            CYSb[idx] = int(line.strip().split(" ")[0])
            CYSc[idx] = int(line.strip().split(" ")[1])
            idx += 1
    
        subprocess.run("cat prot.pdb > temp.pdb", shell=True)

        for res in HEH:
            subprocess.run("cat HEM"+str(res)+".pdb >> temp.pdb", shell=True)
            subprocess.run("cat PRNA"+str(res)+".pdb >> temp.pdb", shell=True)
            subprocess.run("cat PRND"+str(res)+".pdb >> temp.pdb", shell=True)

        subprocess.run("grep -v CRYST1 temp.pdb | grep -v END >"+OutPrefix+"-"+PDB+".pdb", shell=True)

    print("""
 source leaprc.constph
 source leaprc.conste
 source leaprc.water.tip3p

 # Load parameters for ions
 loadAmberParams frcmod.ionsjc_tip3p

 # Load PDB
 %s = loadpdb %s-%s.pdb

 # Connecting HEH to the protein""" %(OutPrefix, OutPrefix, PDB), file=open('tleap.in', 'w'))

    for idx in range(x):
        print("""
 bond %s.%0d.CB1 %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, HISp[idx]), end=" ", file=open('tleap.in', 'a'))
        print("""
 bond %s.%0d.CB2 %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, HISd[idx]), end=" ", file=open('tleap.in', 'a'))
        print("""
 bond %s.%0d.CBB2 %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, CYSb[idx]), end=" ", file=open('tleap.in', 'a'))
        print("""
 bond %s.%0d.CBC1 %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, CYSc[idx]), end=" ", file=open('tleap.in', 'a'))
        print("""
 bond %s.%0d.C2A %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, HEH[idx]+1), end=" ", file=open('tleap.in', 'a'))
        print("""
 bond %s.%0d.C3D %s.%0d.CA""" %(OutPrefix, HEH[idx], OutPrefix, HEH[idx]+2), end=" ", file=open('tleap.in', 'a'))

    while True:
        SolvEnv = input(""" Should the structure be prepared with 
 an explicit or implicit solvent (explicit/implicit)? """)

        if (SolvEnv == "explicit") or (SolvEnv == "Explicit") or (SolvEnv == "e") or (SolvEnv == "E"):
            while True:
                BoxShape = input("  Using a rectangular or an octahedral box (rec/octahed)? ")
                BufferSize = int(input("  With how much of a solvent buffer (in angstroms)? "))
            
                if (BoxShape == "rectangular") or (BoxShape == "rec"):
                    print("""

 #Solvate
 solvateBox %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                    break
                elif (BoxShape == "octahedral") or (BoxShape == "octahed"):
                    print("""

 #Solvate
 solvateOct %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                    break
                else:
                    print("  Sorry, I didn't understand your response.")

            NaCount = int(input("  And how many Na+ ions; 0 = enough for charge neutrality? "))
            ClCount = int(input("  And how many Cl- ions; 0 = enough for charge neutrality? "))

            print("""

 #Add ions
 addions %s Na+ %0d""" %(OutPrefix, NaCount), end=" ", file=open('tleap.in', 'a'))
            print("""
 addions %s Cl- %0d""" %(OutPrefix, ClCount), end=" ", file=open('tleap.in', 'a'))

            break
        elif (SolvEnv == "implicit") or (SolvEnv == "Implicit") or (SolvEnv == "i") or (SolvEnv == "I"):
            break
        else:
            print("  Sorry, I didn't understand your response.")

    print("""

 # Save topology and coordinate files
 saveamberparm %s %s.prmtop %s.rst7""" %(OutPrefix, OutPrefix, OutPrefix), end=" ", file=open('tleap.in', 'a'))

    print("""

 quit""", end=" ", file=open('tleap.in', 'a'))

    print("""
 The re-compiled structure will now be processed with TLEaP.
    """, end=" ")
    subprocess.run("tleap -s -f tleap.in > tleap.log", shell=True)
    print("""
 TLEaP finished! 
 Please inspect the structure to make sure it is correct.
    """)
    return OutPrefix, SolvEnv 

################################################################################################################################################

def StructRelax(Output, Solv):
    if (os.path.isfile(Output+".prmtop") == True and os.path.isfile(Output+".rst7")):
        print(" Found %s.prmtop and %s.rst7" %(Output, Output))
        print(" Preparing to relax the geometry")

        if (Solv == "explicit") or (Solv == "Explicit") or (Solv == "e") or (Solv == "E"):
            print("""
 Energy Minimization Stage in Implicit Solvent
 &cntrl
  imin=1,            ! Perform an energy minimization
  ntb=1,             ! Constant volumne
  cut=10.0,          ! Non-bonded cutoff in angstroms
  ntmin=1,           ! Steepest descent + conjugate gradient method
  ncyc=1000,         ! Number of steepest descent cycles
  maxcyc=5000,       ! Maximum number of minimization cycles
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N&!:WAT|(:HEH,PRN)@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
 /
            """, file=open('min.in', 'w'))

            print(" Running minimization ...")
            subprocess.run("mpirun -np 64 pmemd.MPI -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
            print(" Minimization finished!")
        else:
            pass

        if (Solv == "implicit") or (Solv == "Implicit") or (Solv == "i") or (Solv == "I"):
            print("""
Energy Minimization in Implicit Solvent
&cntrl
  imin=1,            ! Perform an energy minimization
  ntb=0,             ! Non-periodic
  cut=9999           ! Non-bonded cutoff in Å
  ntmin=1,           ! steepest descent + conjugate gradient method 
  ncyc=200,          ! Number of steepest descent cycles
  maxcyc=500,        ! Maximum number of minimization cycles
  igb=2,             ! Generalized Born implicit solvent model
  saltcon=0.1,       ! salt concentration in M
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N|(:HEH,PRN)@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
/
            """, file=open('min.in', 'w'))

            print(" Running minimization ...")
            subprocess.run("sander -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
            print(" Minimization finished!")
    else:
        print(" %s.prmtop and %s.rst7 not found" %(Output, Output), end=" ")
        sys.exit(" Nothing to minimize. Something went wrong in the preceeding steps!")

################################################################################################################################################

def PreparedStructureSelection():
    for x in os.listdir():
        if x.endswith(".prmtop") or x.endswith(".rst7"):
            print(x)

    while True:
        Output = input(" Prefix used for previously generated parm/rst7 ")

        if (os.path.isfile(Output + ".prmtop") == True) and (os.path.isfile(Output + ".rst7") == True):
            print(" That pair of topology and coordinate files were found! ")
            return Output
            break
        elif (os.path.isfile(Output + ".prmtop") == False) or (os.path.isfile(Output + ".rst7") == False):
            sys.exit(""" 
 Either %s.prmtop or %s.rst7 is missing. 
 Please try running the Structure Preparation and Relaxation module 
 before proceeding""" %(Output, Output))
        else: 
            print(""" That pair of topology and coordinate files were NOT found! 
 Please try again.""")

################################################################################################################################################

def ReorderResByChain(Output):
    print("""
parm %s.prmtop
trajin min.rst7
fixatomorder parmout %s_reord.prmtop
trajout %s_reord.rst7 topresnum
run
quit
    """ %(Output, Output, Output), file=open("ReorderRes.in", "w"))

    if (os.path.isfile("%s_reord.prmtop" %(Output)) == True):
        print("\n Found the reordered topology: %s_reord.prmtop!" %(Output))
    if (os.path.isfile("%s_reord.prmtop" %(Output)) == False):
        subprocess.run("cpptraj -i ReorderRes.in > ReorderRes.log 2> /dev/null", shell=True)
#       subprocess.run("ambpdb -p %s_reord.prmtop -c min.rst7 > min.pdb" %(Output), shell=True)

    Output = Output+"_reord"

    return Output

################################################################################################################################################

def LinearizeHemeSequence():

    if (os.path.isfile("LinearizedHemeSequence.txt") == True):
        print(" Found LinearizedHemeSequence.txt!")

    if (os.path.isfile("LinearizedHemeSequence.txt") == False):
        New=list(map(int, input(" Linear Sequence: ").strip().split()))
        x = len(New)
        
        for idx in range(0, x):
            if (idx == 0):
                print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'w'))
            else:
                print(idx, New[idx], file=open('LinearizedHemeSequence.txt', 'a'))

#       idx1 = 0; idx2 = 0;
#       with open("ResIndexing.txt") as fp:
#           x = len(New) #len(fp.readlines())
#           Org = [0]*x
#           Lin = [0]*x

#           fp.seek(0)
#           Lines = fp.readlines()
#           for line in Lines:
#               Org[idx1] = int(line.strip().split(" ")[5])
#               idx1 +=1

#           for idx1 in range(0, x):
#               for idx2 in range(0, x):
#                   if (Org[idx1] == New[idx2]):
#                       Lin[idx2] = Org[idx1]

#           for idx2 in range(0, x):
#               if (idx2 == 0):
#                   print(idx2, New[idx2], file=open('LinearizedHemeSequence.txt', 'w'))
#               else:
#                   print(idx2, New[idx2], file=open('LinearizedHemeSequence.txt', 'a'))

#       print(" Original sequence of hemes: "+str(Org))
#       print(" Linear sequence of hemes: "+str(New)+"\n")

################################################################################################################################################

def LambdaFromSASA(Output):

    if (os.path.isfile(Output+".prmtop") == True and os.path.isfile("min.rst7")):
        print(" Found %s.prmtop and min.rst7" %(Output))
        print("""
 To estimate the reorganization energy (lambda) from the 
 solvent accessible surface area, two steps will be take:
    (1) Convert min.rst7 to a PDB-formatted file using ambpdb""", end=" ")

        if (PolySel == 'Yes') or (PolySel == "yes") or (PolySel == "Y") or (PolySel == "y"):
            print("""
        Note: This step is skipped because you indicated you 
        have a polymeric structute, When the topology was
        re-ordered to conform to AMBER conventions for 
        multi-chain structures, min.pdb was already 
        created.
            """, end=" ")
        print("""
    (2) Write and submit a TCL script to VMD
        """)

        if (PolySel == 'No') or (PolySel == "no") or (PolySel == "N") or (PolySel == "n"):
            subprocess.run("ambpdb -p "+Output+".prmtop -c min.rst7 > min.pdb", shell=True)

        idx = 0
        with open("LinearizedHemeSequence.txt") as fp:
            x = len(fp.readlines())
            HEH = [0]*x

            fp.seek(0)
            Lines = fp.readlines()
            for line in Lines:
                HEH[idx] = int(line.strip().split(" ")[1])
                idx += 1

        for idx in range(len(HEH)-1):
            print("""
 mol new min.pdb
 set HEH1 {%0d}; set HEH2 {%0d}""" %(HEH[idx], HEH[idx+1]), file=open('SASACalc.tcl', 'w'))  
            print("""
 set output [open "${HEH1},${HEH2}_SASAanalysis.dat" a]

 set allsel     [atomselect top "all and not water and not ions"]

 set donor      [atomselect top "resname HEH and resid $HEH1"]
 set Dfe        [atomselect top "resname HEH and resid $HEH1 and name FE"]
 set DfeIDX     [$Dfe get index]
 set DfeID      [$Dfe get resid]

 set acceptor   [atomselect top "resname HEH and resid $HEH2"]
 set Afe        [atomselect top "resname HEH and resid $HEH2 and name FE"]
 set AfeIDX     [$Afe get index]
 set AfeID      [$Afe get resid]
 set DA         [list $DfeIDX $AfeIDX]

 set nf [molinfo top get numframes]
 puts "There are $nf frames"
 for {set frame 0} {$frame < $nf} {incr frame} {
   set time [expr ($frame * 1 * 0.0020)]

   puts "  Analyzing Frame $frame..."

   $allsel   frame $frame; $allsel    update
   $donor    frame $frame; $donor     update
   $acceptor frame $frame; $acceptor  update

   set dsasa [measure sasa 1.4 $allsel -restrict $donor]
   set asasa [measure sasa 1.4 $allsel -restrict $acceptor]
   set Rfefe [measure bond $DA frame $frame]

   puts $output "$time $dsasa $asasa $Rfefe"
 }
 exit
            """, file=open('SASACalc.tcl', 'a'))

            if (os.path.isfile("%0d,%0d_SASAanalysis.dat" %(HEH[idx], HEH[idx+1])) == True):
                print(""" Found %0d,%0d_SASAanalysis.dat from a prior execution.
 This prior output will be used for the analysis.
  """ %(HEH[idx], HEH[idx+1])) 

            if (os.path.isfile("%0d,%0d_SASAanalysis.dat" %(HEH[idx], HEH[idx+1])) == False):
                print(" Now using VMD to compute SASA Donor = %0d & Acceptor = %0d..." %(HEH[idx], HEH[idx+1]))
                subprocess.run("vmd -e SASACalc.tcl > SASACalc.log", shell=True)

        print(" Computing Reorganization Energy from Solvent Accessibility...")
        alpha = 5.18; beta = 0.016;
        Rd=4.6; Ra=4.6;
        Eopt=1.84 

        Dsasa = [0]*(len(HEH)-1)
        Asasa = [0]*(len(HEH)-1)
        TotalSASA = [0]*(len(HEH)-1)
        Es = [0]*(len(HEH)-1)
        M = [0]*(len(HEH)-1)
        Rda = [0]*(len(HEH)-1)
        R = [0]*(len(HEH)-1)
        Lambda = [0]*(len(HEH)-1)
        for idx in range(len(HEH)-1):
            with open(str(HEH[idx])+","+str(HEH[idx+1])+"_SASAanalysis.dat") as fp:
               #print(str(HEH[idx])+","+str(HEH[idx+1])+"_SASAanalysis.dat")
                Lines = fp.readlines()
                for line in Lines:
                   #print(line)
                    Dsasa[idx] = float(line.strip().split(" ")[1])
                    Asasa[idx] = float(line.strip().split(" ")[2])
                    Rda[idx] = float(line.strip().split(" ")[3])
                   #print(Dsasa, Asasa, Rda)

            TotalSASA[idx] = Dsasa[idx] + Asasa[idx]
            Es[idx] = alpha + (beta * TotalSASA[idx])
            M[idx] = (1/Eopt) - (1/Es[idx])
            R[idx] = (1/((2*Rd)/0.53)) + (1/((2*Ra)/0.53)) - (1/(Rda[idx]/0.53))
            Lambda[idx] = ((-1)**2) * (M[idx]) * (R[idx]) * (27.2114)
           #print(TotalSASA[idx], Es[idx], M[idx], R[idx], Lambda[idx])

        for idx in range(len(HEH)-1):
            if (idx == 0):
                print(""" 
 HEH-%0d -> HEH-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" %(HEH[idx], HEH[idx+1], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'w'))

            if (idx != 0):
                print(""" 
 HEH-%0d -> HEH-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" %(HEH[idx], HEH[idx+1], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'a'))

    print(" Done!")

    return Lambda

################################################################################################################################################

def DeltaGFromPBSA(Output, Solv):

    print("""
 Single point PB calculation
 &cntrl
  IPB=2, INP=2, ntx=1, imin=1,
 /

 &pb
  epsin=5.19, epsout=78.2, smoothopt=1, istrng=100, pbtemp=300, radiopt=0, dprob=1.4, iprob=2.0, sasopt=0, saopt=1,
  npbopt=0, solvopt=1, accept=0.001, maxitn=100, fillratio=1.5, space=0.5, nfocus=2, fscale=8, npbgrid=1,
  bcopt=5, eneopt=2, frcopt=2, scalec=0, cutfd=5, cutnb=0,
  !isurfchg=1, npbverb=1,
  !#phiout=1, phiform=2, outlvlset=true,
 /
    """, file=open('pbsa.key', 'w'))

    if (Solv == "explicit") or (Solv == "Explicit") or (Solv == "e") or (Solv == "E"):
        print("""
 parm    %s.prmtop
 trajin  min.rst7
 strip   :WAT,Na+,Cl- outprefix StrucForPBSA
 trajout StrucForPBSA.rst7
 run
        """ %(Output), file=open("GenerateCoordForPBSA.in", 'w'))
        subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)
    else:
        print("""
 parm    %s.prmtop
 trajin  min.rst7
 parmwrite out StrucForPBSA.%s.prmtop
 trajout StrucForPBSA.rst7
 run
        """ %(Output, Output), file=open("GenerateCoordForPBSA.in", 'w'))
        subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)

    idx = 0; idxc = 0;
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x
        command = [' ']*(2*x)
        RedoxState = [0]*2

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[1])

            for state in ["ox", "red"]:

                if (state == "ox"):
                    RedoxState[0] = "o"
                    RedoxState[1] = int(0)
                if (state == "red"):
                    RedoxState[0] = "r"
                    RedoxState[1] = int(1)

                if (os.path.isfile("%s%0d.prmtop" %(RedoxState[0], HEH[idx])) == False):
                    print("""
 parm StrucForPBSA.%s.prmtop
 changeRedoxState :%0d %0d
 outparm %s%0d.prmtop
 quit
                    """ %(Output, HEH[idx], RedoxState[1], RedoxState[0], HEH[idx]), file=open(str(RedoxState[0])+str(HEH[idx])+".inp", 'w'))

                    print("\n Generating topology for %sHEH-%0d ..." %(RedoxState[0], HEH[idx]))
                    subprocess.run("parmed -i "+str(RedoxState[0])+str(HEH[idx])+".inp > "+str(RedoxState[0])+str(HEH[idx])+".log", shell=True)

                if (os.path.isfile("pbsa_"+str(RedoxState[0])+str(HEH[idx])) == True):
                    print(""" 
 Found pbsa_%s%0d from a prior execution.
 This prior output will be used for the analysis.""" %(RedoxState[0], HEH[idx]))
                else: 
                    #print(" Running PBSA calculation for oxidized HEH-%0d ..." %(HEH[idx]))
                    #subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(RedoxState[0])+str(HEH[idx])+" -p "+str(RedoxState[0])+str(HEH[idx])+".prmtop -c StrucForPBSA.rst7", shell=True)
                    command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(RedoxState[0])+str(HEH[idx])+" -p "+str(RedoxState[0])+str(HEH[idx])+".prmtop -c StrucForPBSA.rst7"
                    #print(idx, idxc, command[idx])
                    idxc += 1

            idx += 1

        if (idxc != 0):
            print("\n Submitting "+str(len(command))+" PBSA calculations in parallel")

            procs = [ subprocess.Popen(i, shell=True) for i in command ]
            for p in procs:
                p.wait()
                print("  Finished: "+str(p))

    DEtot = [0]*(len(HEH))
    DEelec = [0]*(len(HEH))
    DG = [0]*(len(HEH)-1)
    for idx in range(len(HEH)):

        idx1 = 0; idx2 = 0;
        with open('pbsa_o'+str(HEH[idx]), 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                word1 = 'Etot'
                word2 = 'EELEC'

                if (line.find(word1) != -1) and (idx1 == 0):
                    EtotOx = float(line.strip().split()[2]) * 0.043 
                    idx1 += 1

                if (line.find(word2) != -1) and (idx2 == 0):
                    EelecOx = float(line.strip().split()[2]) * 0.043
                    idx2 += 1

        idx3 = 0; idx4 = 0;
        with open('pbsa_r'+str(HEH[idx]), 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                word1 = 'Etot'
                word2 = 'EELEC'

                if (line.find(word1) != -1) and (idx3 == 0):
                    EtotRed = float(line.strip().split()[2]) * 0.043
                    idx3 += 1

                if (line.find(word2) != -1) and (idx4 == 0):
                    EelecRed = float(line.strip().split()[2]) * 0.043
                    idx4 += 1

        DEtot[idx] = (EtotOx - EtotRed) 
        DEelec[idx] = (EelecOx - EelecRed) 

        if (idx == 0):
            print("\n Result:")

        print("""  step=%0d HEH-%0d EtotOx=%.3f eV EtotRed=%.3f eV DEtot=%.3f eV EelecOx=%.3f eV EelecRed=%.3f eV DEelec=%.3f eV""" %(idx, HEH[idx], EtotOx, EtotRed, DEtot[idx], EelecOx, EelecRed, DEelec[idx]))

    for idx in range(len(HEH)-1):
        DG[idx] = -1 * ((-1 * DEtot[idx]) + (DEtot[idx+1])) 

        if (idx == 0):
            print("(HEH-%0d = %.3f eV) -> (HEH-%0d = %.3f eV); DG = %10.3f eV" %(HEH[idx],  DEtot[idx], HEH[idx+1], DEtot[idx+1], DG[idx]), file=open('DG.txt', 'w'))
        else:
            print("(HEH-%0d = %.3f eV) -> (HEH-%0d = %.3f eV); DG = %10.3f eV" %(HEH[idx],  DEtot[idx], HEH[idx+1], DEtot[idx+1], DG[idx]), file=open('DG.txt', 'a'))

    return DG

################################################################################################################################################

def DeltaGFromLIE(Output):

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[1])

            for state in ["ox", "red"]:
                if (state == "ox"):
                    print("""
 parm %s.prmtop
 changeRedoxState :%0d 0
 outparm o%0d.prmtop
 quit
                    """ %(Output, HEH[idx],  HEH[idx]), file=open("o"+str(HEH[idx])+".inp", 'w'))

                    if (os.path.isfile("o%0d.prmtop" %(HEH[idx])) == False):
                        print("\n Generating oxidized topology for HEH-%0d ..." %(HEH[idx]))
                        subprocess.run("parmed -i o"+str(HEH[idx])+".inp > r"+str(HEH[idx])+".log", shell=True)

                elif (state == "red"):
                    print("""
 parm %s.prmtop
 changeRedoxState :%0d 1
 outparm r%0d.prmtop
 quit
                    """ %(Output, HEH[idx],  HEH[idx]), file=open("r"+str(HEH[idx])+".inp", 'w'))

                    if (os.path.isfile("r%0d.prmtop" %(HEH[idx])) == False):
                        print(" Generating reduced topology for HEH-%0d ..." %(HEH[idx]))
                        subprocess.run("parmed -i r"+str(HEH[idx])+".inp > r"+str(HEH[idx])+".log", shell=True)

            if (os.path.isfile("o%0d.prmtop" %(HEH[idx])) == True) and (os.path.isfile("r%0d.prmtop" %(HEH[idx])) == True): 
                print("""
 parm   o%(HEH)0d.prmtop
 trajin min.pdb parmindex 0

 lie OHSys      :%(HEH)0d      !:%(HEH)0d                       out  OHSys.dat
 lie OHP        :%(HEH)0d      !(:%(HEH)0d,WAT,Na+,Cl-)         out  OHP.dat
 lie OHS        :%(HEH)0d      :WAT,Na+,Cl-                     out  OHS.dat
 lie OHNonPolar :%(HEH)0d      :GLY,ALA,VAL,LEU,ILE,MET,PRO     out  OHNonPolar.dat
 lie OHAromatic :%(HEH)0d      :PHE,TYR,TRP                     out  OHAromatic.dat
 lie OHPolar    :%(HEH)0d      :SER,THR,CYS,CYO,ASN,GLN         out  OHPolar.dat
 lie OHAcidic   :%(HEH)0d      :ASP,AS4,GLU,GL4                 out  OHAcidic.dat 
 lie OHBasic    :%(HEH)0d      :HIS,HIP,HID,HIE,HIO,LYS,ARG     out  OHBasic.dat 
 lie OHHEH      :%(HEH)0d      :HEH&!:%(HEH)0d                  out  OHHEH.dat 
 lie OHPRN      :%(HEH)0d      :PRN                             out  OHPRN.dat 
 run

 clear trajin

 parm   r%(HEH)0d.prmtop
 trajin min.pdb parmindex 1

 lie RHSys      :%(HEH)0d      !:%(HEH)0d                       out  RHSys.dat
 lie RHP        :%(HEH)0d      !(:%(HEH)0d,WAT,Na+,Cl-)         out  RHP.dat
 lie RHS        :%(HEH)0d      :WAT,Na+,Cl-                     out  RHS.dat
 lie RHNonPolar :%(HEH)0d      :GLY,ALA,VAL,LEU,ILE,MET,PRO     out  RHNonPolar.dat
 lie RHAromatic :%(HEH)0d      :PHE,TYR,TRP                     out  RHAromatic.dat
 lie RHPolar    :%(HEH)0d      :SER,THR,CYS,CYO,ASN,GLN         out  RHPolar.dat
 lie RHAcidic   :%(HEH)0d      :ASP,AS4,GLU,GL4                 out  RHAcidic.dat 
 lie RHBasic    :%(HEH)0d      :HIS,HIP,HID,HIE,HIO,LYS,ARG     out  RHBasic.dat 
 lie RHHEH      :%(HEH)0d      :HEH&!:%(HEH)0d                  out  RHHEH.dat 
 lie RHPRN      :%(HEH)0d      :PRN                             out  RHPRN.dat 
 run
#-----------------------------------------------------------------------------------

 ODHSys      = avg(OHSys[EELEC])      * 0.043 
 ODHP        = avg(OHP[EELEC])        * 0.043 
 ODHS        = avg(OHS[EELEC])        * 0.043
 ODHNonPolar = avg(OHNonPolar[EELEC]) * 0.043
 ODHAromatic = avg(OHAromatic[EELEC]) * 0.043
 ODHPolar    = avg(OHPolar[EELEC])    * 0.043
 ODHAcidic   = avg(OHAcidic[EELEC])   * 0.043
 ODHBasic    = avg(OHBasic[EELEC])    * 0.043
 ODHHEH      = avg(OHHEH[EELEC])      * 0.043
 ODHPRN      = avg(OHPRN[EELEC])      * 0.043

 RDHSys      = avg(RHSys[EELEC])      * 0.043 
 RDHP        = avg(RHP[EELEC])        * 0.043 
 RDHS        = avg(RHS[EELEC])        * 0.043
 RDHNonPolar = avg(RHNonPolar[EELEC]) * 0.043
 RDHAromatic = avg(RHAromatic[EELEC]) * 0.043
 RDHPolar    = avg(RHPolar[EELEC])    * 0.043
 RDHAcidic   = avg(RHAcidic[EELEC])   * 0.043
 RDHBasic    = avg(RHBasic[EELEC])    * 0.043
 RDHHEH      = avg(RHHEH[EELEC])      * 0.043
 RDHPRN      = avg(RHPRN[EELEC])      * 0.043

 DHSys          = (avg(OHSys[EELEC])      - avg(RHSys[EELEC]))      * 0.043    
 DHP            = (avg(OHP[EELEC])        - avg(RHP[EELEC]))        * 0.043 
 DHS            = (avg(OHS[EELEC])        - avg(RHS[EELEC]))        * 0.043 
 DHNonPolar     = (avg(OHNonPolar[EELEC]) - avg(RHNonPolar[EELEC])) * 0.043    
 DHAromatic     = (avg(OHAromatic[EELEC]) - avg(RHAromatic[EELEC])) * 0.043 
 DHPolar        = (avg(OHPolar[EELEC])    - avg(RHPolar[EELEC]))    * 0.043    
 DHAcidic       = (avg(OHAcidic[EELEC])   - avg(RHAcidic[EELEC]))   * 0.043 
 DHBasic        = (avg(OHBasic[EELEC])    - avg(RHBasic[EELEC]))    * 0.043
 DHHEH          = (avg(OHHEH[EELEC])      - avg(RHHEH[EELEC]))      * 0.043
 DHPRN          = (avg(OHPRN[EELEC])      - avg(RHPRN[EELEC]))      * 0.043

 #DHSysstd       = stdev(DHSys)       * 0.043    
 #DHPstd         = stdev(DHP)         * 0.043 
 #DHSstd         = stdev(DHS)         * 0.043 
 #DHNonPolarstd  = stdev(DHNonPolar)  * 0.043
 #DHAromaticstd  = stdev(DHAromatic)  * 0.043 
 #DHPolarstd     = stdev(DHPolar)     * 0.043 
 #DHAcidicstd    = stdev(DHAcidic)    * 0.043 
 #DHBasicstd     = stdev(DHBasic)     * 0.043 
 #DHHEHstd       = stdev(DHHEH)       * 0.043
 #DHPRNstd       = stdev(DHPRN)       * 0.043
 #writedata %(HEH)0dphsicochemStd.dat DHSysstd DHPstd DHSstd DHNonPolarstd DHAromaticstd DHPolarstd DHAcidicstd DHBasicstd DHHEHstd DHPRNstd 
 #printdata DHSysstd DHPstd DHSstd DHNonPolarstd DHAromaticstd DHPolarstd DHAcidicstd DHBasicstd DHHEHstd DHPRNstd 

 printdata ODHSys RDHSys DHSys ODHP RDHP DHP ODHS RDHS DHS ODHNonPolar RDHNonPolar DHNonPolar ODHAromatic RDHAromatic DHAromatic ODHPolar RDHPolar DHPolar ODHAcidic RDHAcidic DHAcidic ODHBasic RDHBasic DHBasic ODHHEH RDHHEH DHHEH ODHPRN RDHPRN DHPRN 
 writedata %(HEH)0dphsicochemAvg.dat ODHSys RDHSys DHSys ODHP RDHP DHP ODHS RDHS DHS ODHNonPolar RDHNonPolar DHNonPolar ODHAromatic RDHAromatic DHAromatic ODHPolar RDHPolar DHPolar ODHAcidic RDHAcidic DHAcidic ODHBasic RDHBasic DHBasic ODHHEH RDHHEH DHHEH ODHPRN RDHPRN DHPRN
                """ %{'HEH': HEH[idx]}, file=open("%0dredox.in" %(HEH[idx]), 'w'))

                if (os.path.isfile("%0dphsicochemAvg.dat" %(HEH[idx])) == True):
                    print(""" 
 Found %0dphsicochemAvg.dat from a prior execution.
 This prior output will be used for the analysis.""" %(HEH[idx])) 

                if (os.path.isfile("%0dphsicochemAvg.dat" %(HEH[idx])) == False):
                    print(" Running LIE calculation for Oxidation of HEH-%0d ..." %(HEH[idx]))
                    subprocess.run("cpptraj -i "+str(HEH[idx])+"redox.in > "+str(HEH[idx])+"redox.log", shell=True)

            else:
                print(""" o%0d.prmtop and/or r%0d.prmtop are missing and we can not proceed. Something went wrong.""" %(HEH[idx], HEH[idx]))

            idx += 1

    ODHSys = [0]*(len(HEH))
    RDHSys = [0]*(len(HEH))
    DHSys = [0]*(len(HEH))

    ODHP = [0]*(len(HEH))
    RDHP = [0]*(len(HEH))
    DHP = [0]*(len(HEH))

    ODHS = [0]*(len(HEH))
    RDHS = [0]*(len(HEH))
    DHS = [0]*(len(HEH))

    ODHNonPolar = [0]*(len(HEH))
    RDHNonPolar = [0]*(len(HEH))
    DHNonPolar = [0]*(len(HEH))

    ODHAromatic = [0]*(len(HEH))
    RDHAromatic = [0]*(len(HEH))
    DHAromatic = [0]*(len(HEH))

    ODHPolar = [0]*(len(HEH))
    RDHPolar = [0]*(len(HEH))
    DHPolar = [0]*(len(HEH))

    ODHAcidic = [0]*(len(HEH))
    RDHAcidic = [0]*(len(HEH))
    DHAcidic = [0]*(len(HEH))

    ODHBasic = [0]*(len(HEH))
    RDHBasic = [0]*(len(HEH))
    DHBasic = [0]*(len(HEH))

    ODHHEH = [0]*(len(HEH))
    RDHHEH = [0]*(len(HEH))
    DHHEH = [0]*(len(HEH))

    ODHPRN = [0]*(len(HEH))
    RDHPRN = [0]*(len(HEH))
    DHPRN = [0]*(len(HEH))

    DG = [0]*(len(HEH)-1)

    for idx in range(len(HEH)):
        print("\n Result for HEH-%0d:" %(HEH[idx]))

        with open("%0dphsicochemAvg.dat" %(HEH[idx]), 'r') as fp:
            LinesNumToRead = [1]
            LinesToRead = []*len(LinesNumToRead)
            for i, line in enumerate(fp):
                if i in LinesNumToRead:
                    LinesToRead.append(line.strip())

                    ODHSys[idx] = float(LinesToRead[0].split()[1])
                    RDHSys[idx] = float(LinesToRead[0].split()[2])
                    DHSys[idx] = float(LinesToRead[0].split()[3])
                    print("  Full System: %8.3f %8.3f %8.3f" %(ODHSys[idx], RDHSys[idx], DHSys[idx]))
            
                    ODHP[idx] = float(LinesToRead[0].split()[4])
                    RDHP[idx] = float(LinesToRead[0].split()[5])
                    DHP[idx] = float(LinesToRead[0].split()[6])
                    print("      Protein: %8.3f %8.3f %8.3f" %(ODHP[idx], RDHP[idx], DHP[idx]))

                    ODHS[idx] = float(LinesToRead[0].split()[7])
                    RDHS[idx] = float(LinesToRead[0].split()[8])
                    DHS[idx] = float(LinesToRead[0].split()[9])
                    print("      Solvent: %8.3f %8.3f %8.3f" %(ODHS[idx], RDHS[idx], DHS[idx]))

                    ODHNonPolar[idx] = float(LinesToRead[0].split()[10])
                    RDHNonPolar[idx] = float(LinesToRead[0].split()[11])
                    DHNonPolar[idx] = float(LinesToRead[0].split()[12])
                    print("    Non-Polar: %8.3f %8.3f %8.3f" %(ODHNonPolar[idx], RDHNonPolar[idx], DHNonPolar[idx]))

                    ODHAromatic[idx] = float(LinesToRead[0].split()[13])
                    RDHAromatic[idx] = float(LinesToRead[0].split()[14])
                    DHAromatic[idx] = float(LinesToRead[0].split()[15])
                    print("     Aromatic: %8.3f %8.3f %8.3f" %(ODHAromatic[idx], RDHAromatic[idx], DHAromatic[idx]))

                    ODHPolar[idx] = float(LinesToRead[0].split()[16])
                    RDHPolar[idx] = float(LinesToRead[0].split()[17])
                    DHPolar[idx] = float(LinesToRead[0].split()[18])
                    print("        Polar: %8.3f %8.3f %8.3f" %(ODHPolar[idx], RDHPolar[idx], DHPolar[idx]))

                    ODHAcidic[idx] = float(LinesToRead[0].split()[19])
                    RDHAcidic[idx] = float(LinesToRead[0].split()[20])
                    DHAcidic[idx] = float(LinesToRead[0].split()[21])
                    print("       Acidic: %8.3f %8.3f %8.3f" %(ODHAcidic[idx], RDHAcidic[idx], DHAcidic[idx]))

                    ODHBasic[idx] = float(LinesToRead[0].split()[22])
                    RDHBasic[idx] = float(LinesToRead[0].split()[23])
                    DHBasic[idx] = float(LinesToRead[0].split()[24])
                    print("        Basic: %8.3f %8.3f %8.3f" %(ODHBasic[idx], RDHBasic[idx], DHBasic[idx]))

                    ODHHEH[idx] = float(LinesToRead[0].split()[25])
                    RDHHEH[idx] = float(LinesToRead[0].split()[26])
                    DHHEH[idx] = float(LinesToRead[0].split()[27])
                    print("  Other Hemes: %8.3f %8.3f %8.3f" %(ODHHEH[idx], RDHHEH[idx], DHHEH[idx]))

                    ODHPRN[idx] = float(LinesToRead[0].split()[28])
                    RDHPRN[idx] = float(LinesToRead[0].split()[29])
                    DHPRN[idx] = float(LinesToRead[0].split()[30])
                    print("  Propionates: %8.3f %8.3f %8.3f" %(ODHPRN[idx], RDHPRN[idx], DHPRN[idx]))
                elif i > 1:
                    break

    for idx in range(len(HEH)-1):
        DG[idx] = -1 * ((-1 * DHSys[idx]) + (DHSys[idx+1])) 

        if (idx == 0):
            print("(HEH-%0d = %.3f eV) -> (HEH-%0d = %.3f eV); DG = %10.3f eV" %(HEH[idx],  DHSys[idx], HEH[idx+1], DHSys[idx+1], DG[idx]), file=open('DG.txt', 'w'))
        else:
            print("(HEH-%0d = %.3f eV) -> (HEH-%0d = %.3f eV); DG = %10.3f eV" %(HEH[idx],  DHSys[idx], HEH[idx+1], DHSys[idx+1], DG[idx]), file=open('DG.txt', 'a'))

    return DG

################################################################################################################################################

def AssignCouplingFromGeom(Output):

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[1])
            idx += 1

    ang = [0]*(len(HEH)-1)
    Hda = [0]*(len(HEH)-1)
    for idx in range(len(HEH)-1):
        print("""
 parm    %s.prmtop
 trajin  min.rst7
 
 vector Hem%0d  corrplane :%0d&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 vector Hem%0d  corrplane :%0d&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 run
 
 vectormath vec1 Hem%0d vec2 Hem%0d dotangle out hemplaneorient_%0d-%0d.dat
 
 run
 quit
        """ %(Output, HEH[idx], HEH[idx], HEH[idx+1], HEH[idx+1], HEH[idx], HEH[idx+1], HEH[idx], HEH[idx+1]), file=open("CalcOrientation.in", "w"))

        if (os.path.isfile("hemplaneorient_%0d-%0d.dat" %(HEH[idx], HEH[idx+1])) == True):
            print(""" 
 Found hemplaneorient_%0d-%0d.dat from a prior execution.
 This prior output will be used for the analysis.""" %(HEH[idx], HEH[idx+1])) 

        if (os.path.isfile("hemplaneorient_%0d-%0d.dat" %(HEH[idx], HEH[idx+1])) == False):
            subprocess.run("cpptraj -i CalcOrientation.in > CalcOrientation.log", shell=True)

        lc = 0
        with open("hemplaneorient_%0d-%0d.dat" %(HEH[idx], HEH[idx+1])) as fp:
            Lines = fp.readlines()
            for line in Lines:
                if (lc == 1):
                    ang[idx] = float(line.strip().split()[1])
                    if (ang[idx] < 45):
                        Hda[idx] = 8.0
                    elif (ang[idx] >= 45):
                        Hda[idx] = 2.0
                    else:
                        Hda[idx] = "?"
                lc += 1

    print("\n Assiging coupling values based on inter-macrocycle planar anlge ... ")
    for idx in range(len(HEH)-1):

        if (idx == 0):
            print("Hda(HEH-%0d <-> HEH-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" %(HEH[idx],  HEH[idx+1], ang[idx], Hda[idx]), file=open('Hda.txt', 'w'))
        else:
            print("Hda(HEH-%0d <-> HEH-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" %(HEH[idx],  HEH[idx+1], ang[idx], Hda[idx]), file=open('Hda.txt', 'a'))

    print(" Done!")

    return Hda

################################################################################################################################################

def ComputeMarcusRates(Hda, Lambda, DG):

    NumSteps = len(DG)
    prefactor = [0]*(NumSteps)
    Eactf = [0]*(NumSteps)
    ketf = [0]*(NumSteps)
    Eactb = [0]*(NumSteps)
    ketb = [0]*(NumSteps)
    for idx in range(NumSteps):
        T = 300.0
        pi = 3.141592654
        kB = 8.6173304E-5
        hbar = 6.582119514E-16

        prefactor[idx] = (2 * pi * (Hda[idx]/1000)**2) / (hbar * math.sqrt(4 * pi * Lambda[idx] * kB* T))

        Eactf[idx] = ((DG[idx] + Lambda[idx])**2)/(4 * Lambda[idx])
        ketf[idx] = prefactor[idx] * math.exp((-1 * Eactf[idx])/(kB * T))

        Eactb[idx] = (((-1 * DG[idx]) + Lambda[idx])**2)/(4 * Lambda[idx])
        ketb[idx] = prefactor[idx] * math.exp((-1 * Eactb[idx])/(kB * T))

    for idx in range(NumSteps):
        print("""
 Step #%0d:
    Activation Energy:
        Forward: %.3E eV
        Reverse: %.3E eV
    Rates:
        Forward: %.3E 
        Reverse: %.3E\n""" %(idx, Eactf[idx], Eactb[idx], ketf[idx], ketb[idx]), end=" ")

        if (idx == 0):
            print("ketf,ketb", file=open("rates.txt", "w"))
            print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "a"))
        else:
            print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "a"))

################################################################################################################################################

def ComputeDiffusionCoefficient():
#   subprocess.run("make 2> /dev/null", shell=True)
#   subprocess.run("./derrida 2> /dev/null", shell=True)

    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()
    V,D = derrida.VD(ketf, ketb)
    print(" Diffusion constant = %E cm^2/S" % (2.5e-15 * D))
    print("Diffusion constant = %E (cm^2/S" % (2.5e-15 * D), file=open('D.txt', 'w'))
################################################################################################################################################

def ComputeFlux():
    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()

    Jf,Jb = blumberger.flux(ketf, ketb)
    print("Forward Flux: %.2E" %(Jf)) 
    print("Reverse Flux: %.2E" %(Jb))
################################################################################################################################################

def MeasureSubunitLength():

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[1])
            idx += 1
    print("""
 mol new min.pdb
 set output [open  "SubunitLength.txt" w]

 set FirstHeme [[atomselect top "resname HEH and resid %0d and name FE"] get index]
 set LastHeme  [[atomselect top "resname HEH and resid %0d and name FE"] get index]

 set dista  [measure bond [list $FirstHeme $LastHeme]]
 set distcm [expr {$dista * 1E-8}]

 puts $output "Proposed Subunit Length (cm) = $distcm"
 exit
    """ %(HEH[0], HEH[x-1]), file=open("MeasureSubunitLength.tcl", "w"))
    subprocess.run("vmd -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)

    with open("SubunitLength.txt") as fp:
        SubunitLength = float(fp.readline().strip().split()[5])

    return SubunitLength

################################################################################################################################################

def ComputeRedoxCurrent():

    r=0.75E-7      #cm
    e=1.602E-19    #c
    kb=1.38E-23    #J/k

    print(" Please provide the following parmaeters: ")
    T = float(input("  Temperature (K)? "))
    cps = int(input("  Number of Charges per subunit? "))
    SubunitLength = MeasureSubunitLength()
#   lsub = float(input("  Length of subunit (cm)? "))
    lsub = SubunitLength
    print(""" 
 The length of a subunit of the cytochrome polymer is needed.
  The subunit length measured between the first and the last heme
  specified in LinearizedHemeSequence.txt is %.2E \n""" %(SubunitLength))

    lw = float(input("  Length of wire (cm)? "))
    Gexp = float(input("""  Experimental Conductance (S/cm) ? \n   (Enter "0" if not known)         """))
    
    if (os.path.isfile("D.txt") == True):
        with open("D.txt") as fp:
            Dcalc = float(fp.readline().strip().split()[3])

    cpsul = ((cps)/(lsub))
    csa = (math.pi * (r)**2) 
    crgden = (cpsul)/(csa)
    Dexp = ((Gexp * kb * T * lw) / (csa * crgden * (e)**2))

    print("""
 Charge per Subunit Length          = %e q/cm 
 Cross-Sectional Area               = %e cm^2
 Charge Density                     = %e q/cm^2
 Experimental Diffusion Constant    = %e cm^2/s
 Computed Diffusion Constant        = %e cm^2/s
    """ %(cpsul, csa, crgden, Dexp, Dcalc))

    print(" %6s %8s %8s" %("Voltage (V)", "Exp. Current (pA)", "Computed Current (pA)"))
    for V in np.arange (-0.5, 0.5, 0.05):
        Iexp = ((Gexp * V) * 1E12)
        Icmp = (((csa * crgden * (e)**2 * Dcalc) / (kb * T * lw)) * V) * 1E12 
        print(" %11.3f %17.3f %21.3E" %(V, Iexp, Icmp))

################################################################################################################################################

print("""
 ================================================================== 
                      Weclome to BioDC
           A program that automates and accelerates
             the computation of redox currents in
              (polymeric) multi-heme cytochormes 

            Written by Matthew J. Guberman-Pfeffer
                Last Updated: 08/30/2023

 ================================================================== 

 BioDC presents a highly modular workflow that has three 
 large divisions: 
    (1) Structure Preparaiton & Relaxation
    (2) Energetic Estimation
    (3) Redox Current Prediction
""", end=" ")

while True:
    try:
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

if (DivSel == 0) or (DivSel == 1):
    print("""
 ===================================================================  
 First: Structure Preparation & Relaxation
 ===================================================================  
    """, end=" ")

    CheckProgInPath()
    PDB = InitialStructureSelection()
    CreateResIndexing(PDB)
    ProcessPDB(PDB)
    Output, Solv = ReBuildStructure(PDB)
    StructRelax(Output, Solv)

if (DivSel == 0) or (DivSel == 2):
    print("""
 ===================================================================  
 Second: Energetic Estimation
 ===================================================================
""")

    if (DivSel == 2):
        Output = PreparedStructureSelection()

        if (os.path.isfile(Output+".prmtop") == True) and (os.path.isfile(Output+".rst7") == True):
            if (os.path.isfile("min.rst7") == True):
                print("""
 Found %s.prmtop and %s.rst7, and the minimized structure (min.rst7).
 We are all set to proceed!""" %(Output, Output))
            else:
                while True:
                    MinSel = input("""
 Found %s.prmtop and %s.rst7, but min.rst7.
 Would you like to relax the structure (yes/no)? """ %(Output, Output))

                    if (MinSel == "Yes") or (MinSel == "yes") or (MinSel == "Y") or (MinSel == "y"):
                        StructRelax(Output, Solv)
                        break
                    elif (MinSel == "No") or (MinSel == "no") or (MinSel == "N") or (MinSel == "n"):
                        sys.exit("""
 Use of a structually relaxed geometry is hard-corded into the 
 program. If you wish to override this best-practice, please
 rename your prepared %s.rst7 to min.rst7 and re-run this
 module of BioDC.       """ %(Output))
                    else:
                        print(" Sorry, I didn't understand your response.")

        if (os.path.isfile(Output+".prmtop") == False) or (os.path.isfile(Output+".rst7") == False):
            sys.exit("""
 Sorry, but %s.prmtop and/or %s.rst7 is/are missing! 
 Please run the Structure Preparation and Relaxation 
 module, or if you already ran it, please check the 
 output from TLEaP for errors when building the 
 topology.

 If you simply made a typo in your response, please
 start again.""" %(Output, Output))

        while True:
            Solv = input("\n Is an explicit or implicit solvent present (explicit/implicit)? ")

            if (Solv == "Explicit") or (Solv == "Implicit") or (Solv == "explicit") or (Solv == "implicit") or (Solv == "E") or (Solv == "I") or (Solv == "e") or (Solv == "i"):
                pass
                break
            else:
                print("""\n Sorry, I didn't understand your selection for the type of solvent used.""")

    while True:
        PolySel = input("""
 Is your structure polymeric (yes/no)? """)

        if (PolySel == 'Yes') or (PolySel == "yes") or (PolySel == "Y") or (PolySel == "y"):
            print("""
 The structure preparation stage required you to place all the heme
 residues from all the chains at the end of the PDB with sequential
 numbering. The programs in AmberTools that will be used to estimate 
 the charge transfer energetics want instead the residues to come in 
 the order of the connectiviity; that is, the hemes of chain A should 
 come both any residue in chain B. 

 To oblige this different numbering convention, we'll use the 
 CPPTRAJ of the AmberTools package to re-order the residues. 
 This process will write a new topology and coordinate file, 
 where the latter is of the structure you previously minimized.""")
            Output = ReorderResByChain(Output)
            break

        elif (PolySel == 'No') or (PolySel == "n") or (PolySel == "N") or (PolySel == "n"):
            print("""
 Great! Thanks for the clarification.""")
            break

        else:
            print(" Sorry, I didn't understand your response.")

    print("""
 To compute the energetics for heme-to-heme electron transfer. We 
 need to know the linear sequence of hemes that will serve as 
 charge hopping sites. Typically the linear sequence is NOT the 
 sequence of residue IDs in the PDB. We therefore need to specify
 the linear sequence to compute the right electron transfer steps. 
    """)
    LinearizeHemeSequence()

    print(" ------------------------------------------------------------------- ")

    while True:
        CompLambda = input("""
 Should we compute the reorganization energy (yes/no)? """)

        if (CompLambda == 'Yes') or (CompLambda == "yes") or (CompLambda == "Y") or (CompLambda == "y"):
            Lambda = LambdaFromSASA(Output)
            break
        elif (CompLambda == 'No') or (CompLambda == "no") or (CompLambda == "N") or (CompLambda == "n"):
            print("""
 An array where each eleemnt is the lambda for a charge transfer is needed.""")

            NumSteps = int(input(" Enter the total number of charge transfer step? "))
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

        if (CompDG == 'Yes') or (CompDG == "yes") or (CompDG == "Y") or (CompDG == "y"):

            if (Solv == "explicit") or (Solv == "Explicit") or (Solv == "e") or (Solv == "E"):
                print("""
 Two different methods are implemented to estimate heme redox 
 potentials and thereby reaction free energies when an explict
 solvent is presnet:

    (1) Compute the change in electrostatic interaciton energy upon 
        heme oxidaiton using the Linear Interaction Energy method 
        in CPPTRAJ.

    (2) Compute the change in electrostatic interaciton energy upon 
        heme oxidaiton using the Poisson–Boltzmann Surface Area
        method implemented in AmberTools (essentially a Delphi-type
        calculation). In this case, the explicit solvent 
        prepared with the sturcutre is discarded.

    Note that method #1 has two advantages:
        (1) It is considerably faster than method #2

        (2) The overall change in electorstatic energy is decomposed
            into contributions from different groups of residues.
                """)

                while True:
                    DGmethod = input(""" Should we use the LIE or PBSA method (lie/pbsa)? """)

                    if (DGmethod == "LIE") or (DGmethod == "lie") or (DGmethod == str(1)):
                        DG = DeltaGFromLIE(Output)
                        break

                    elif (DGmethod == "PBSA") or (DGmethod == "pbsa") or (DGmethod == str(2)):
                        DG = DeltaGFromPBSA(Output, Solv)
                        break
                    else:
                        print(" Sorry, I didn't understand your response.\n")
                break

            if (Solv == "implicit") or (Solv == "Implicit") or (Solv == "i") or (Solv == "I"):
                print("""
 The method implemented to estimate heme redox potentials and thereby
 reaction free energies when an explicit solvent is NOT presnet is 
 to compute the change in electrostatic interaciton energy upon 
 heme oxidaiton using the Poisson–Boltzmann Surface Area
 module in the AmberTools package (essentially a Delphi-type
 calculation).""")
                DG = DeltaGFromPBSA(Output, Solv)
                break
        elif (CompDG == 'No') or (CompDG == "no") or (CompDG == "N") or (CompDG == "n"):
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
        CompHda = input("""
 Should we estimate the electronic coupling from the geometry (yes/no)? """)

        if (CompHda == 'Yes') or (CompHda == "yes") or (CompHda == "Y") or (CompHda == "y"):
            Hda = AssignCouplingFromGeom(Output)
            break
        elif (CompHda == 'No') or (CompHda == "no") or (CompHda == "N") or (CompHda == "n"):
            print("""
 An array where each eleemnt is the Hda for a charge transfer is needed.""")

            NumSteps = int(input(" Enter the total number of charge transfer step? "))
            print("""
 Enter Hda for each charge transfer step followed by return """)

            Hda = [0]*(NumSteps)
            for idx in range(0, NumSteps):
                Hda[idx] = float(input(""" """))
            break
        else:
            print(" Sorry, I didn't understand your response.")

    print(" ------------------------------------------------------------------- ")

    while True:
        CompKet = input("""
 Should we compute the non-adiabatic Marcus-theory rates (yes/no)? """)

        if (CompKet == 'Yes') or (CompKet == "yes") or (CompKet == "Y") or (CompKet == "y"):
            ComputeMarcusRates(Hda, Lambda, DG)
            break
        elif (CompKet == 'No') or (CompKet == "no") or (CompKet == "N") or (CompKet == "n"):
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

if (DivSel == 0) or (DivSel == 3):
    print("""
 ===================================================================  
 Third: Redox Current Prediction
 ===================================================================  
    """, end=" ")

    print("""
 This division of the BioDC workflow computes, via the analytical
 Derrida formula [Ref #1], the charge diffuction coefficient 
 based on the non-adiabatic Marcus theory rates. The diffusion 
 coefficient is then related to the electrical resistance and used 
 in Ohm's law to compute the current as a function of applied bias
 [Ref #2]. Note that this approach is only rigorously correct in 
 the limit of zero bias and single (or non-interacting) mobile 
 charges. The implementaiton of the Derrida formula was kindly 
 provided by Fredrik Jansson [Refs #3 & #4].

 To relax the single-particle condition, the multi-particle 
 steady-state flux is computed according to an analytical 
 expression derived by Blumberger and co-workers. The original code 
 was kindly provided by Jochen Blumberger and Xiuyun Jiang.

 References
   [1] B. Derrida, 
       "Velocity and diffusion constant of a periodic one-dimensional hopping model" 
       J. Stat. Phys. 31, 433 (1983).

   [2] M. J. Guberman-Pfeffer
       "Assessing Thermal Response of Redox Conduction for Anti-Arrhenius Kinetics 
       in a Microbial Cytochrome Nanowire"
       J. Phys. Chem. B 2022, 126, 48, 10083–10097

   [3] Thesis, F. Jansson, Charge transport in disordered materials -
       simulations, theory, and numerical modeling of hopping transport and
       electron-hole recombination. Åbo Akademi University, 2011
       https://urn.fi/URN:NBN:fi-fe201311277464

       Implementation details in section 3.6, especially a method to evaluate
       the expressions in linear time. Application in section 6.2.

   [4] Effect of Electric Field on Diffusion in Disordered Materials I.
       One-dimensional Hopping Transport, A. V. Nenashev, F. Jansson,
       S. D. Baranovskii, R. Österbacka, A. V. Dvurechenskii, F. Gebhard,
       Phys. Rev. B 81, 115203 (2010)
       http://dx.doi.org/10.1103/PhysRevB.81.115203

   [5] M. Breuer, K. M. Rosso, and J. Blumberger, 
       “Electron flow in multi-heme bacterial cytochromes is a balancing act 
       between heme electronic interaction and redox potentials,”
       Proc. Nat. Acad. Sci. USA, vol. 111, p. 611, 2014.

   [6] X. Jiang, Z. Futera, M. E. Ali, F. Gajdos, G. F. von Rudorff, A. Carof, M. Breuer, and J. Blumberger, 
       “Cysteine linkages accelerate electron flow through tetra-heme protein STC,” 
       J. Am. Chem. Soc., vol. 139, p. 17237–17240, 2017.

   [7] X. Jiang, J. H. van Wonderen, J. N. Butt, M. J. Edwards, T. A. Clarke, and J. Blumberger, 
       “Which multi-heme protein complex transfers electrons more efficiently? Comparing MtrCAB from Shewanella 
       with OmcS from Geobacter,” 
       J. Phys. Chem. Lett., vol. 11, pp. 9421-9425, 2020.
    """, end=" ")

    if (os.path.isfile("rates.txt") == True):
        print("""
 Found rates.txt, which is needed to proceed!

 We will now compute the multi-particle, stead-state flux.
 """)
        ComputeFlux()

        print("""
 We will now compute the single-particle diffusion coefficient. 
 """)
        ComputeDiffusionCoefficient()
    else:
        sys.exit("""
 rates.txt not found! Please run division #2 of the BioDC 
 program to generate rates.txt, or create the file by hand.
 In the latter case, the syntax on each line should be: 
    [forward rate],[revers rate]
 where each line is for a charge transfer step. The lines
 are assumed to be in the linear sequence of charge 
 transfer steps. 
        """)

    if (os.path.isfile("D.txt") == True):
        print("""
 Found D.txt, which is needed to proceed!
 We will at last compute the redox current. 
 To do this, some system-specific information is needed.
        """)
        ComputeRedoxCurrent()
        print(" Done!")
    else:
        sys.exit("""
 D.txt not found! Please re-run division #3 of BioDC or
 create the file by hand. In the latter case, the syntax
 should be:
    Diffusion constant = [D] cm^2/S
 where "[D]" is to be replaced with the diffusion coefficient.
        """)



