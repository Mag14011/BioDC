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

def Initialization():

    while True:
        ProgInPath = input("""
 This program requires VMD and the Amber Molecular Dynamics Suite 
 to be in your system's PATH variable. Are they? (yes/no)? """)

        if (ProgInPath == 'Yes') or (ProgInPath == "yes") or (ProgInPath == "Y") or (ProgInPath == "y"):
            print(""" 
 Good! Now, here are the PDBs in the present direcotry:\n""")

            for x in os.listdir():
                if x.endswith(".pdb"):
                    print(x)

            while True:
                OriginalPDB = input("""
 Which PDB would you like to setup 
 (omit the .pdb file extension)? """)

                if (os.path.isfile(OriginalPDB + ".pdb") == True):
                    print("\n That PDB was found! ")
                    return OriginalPDB
                    break
                else: 
                    print("\n That PDB does not exist unfortunately")
            break
        elif (ProgInPath == "No") or (ProgInPath == "no") or (ProgInPath == "N") or (ProgInPath == "n"):
            sys.exit("""
 Please make VMD and the Amber Molecular Dynamics Suite 
 findalbe in your system PATH variable. Then, please 
 re-run this program \n""")
        else: 
            print("\n Sorry, I didn't understand your response.")

################################################################################################################################################

def Mutate(PDB):

    OrigResName = input(" What is the three-letter code residue name before the mutaiton? ")
    MutResID = int(input(" What is the residue ID of the point mutation? "))
    MutResName = input(" What is the three-letter code residue name after the mutaiton? ")

    if (os.path.isfile(f"{OrigResName}-{MutResID}-{MutResName}.tcl") == True):
        subprocess.run(f"vmd -e {OrigResName}-{MutResID}-{MutResName}.tcl > {OrigResName}-{MutResID}-{MutResName}.log", shell=True)
    elif (os.path.isfile(f"{OrigResName}-{MutResID}-{MutResName}.tcl") == False):

        print(f"""
mol new {PDB}.pdb

set res [atomselect top "resname {OrigResName} and resid {MutResID}"]
$res set resname {MutResName}

set mut [atomselect top "(all and not resid {MutResID}) or (not sidechain and resid {MutResID})"]
$mut writepdb {OrigResName}-{MutResID}-{MutResName}.pdb
exit
        """, file=open(f"{OrigResName}-{MutResID}-{MutResName}.tcl", 'w'))

        print(f"\n Generating PDB for {OrigResName}-{MutResID}-{MutResName} mutant ...")
        subprocess.run(f"vmd -e {OrigResName}-{MutResID}-{MutResName}.tcl > {OrigResName}-{MutResID}-{MutResName}.log", shell=True)
    PDB = f"{OrigResName}-{MutResID}-{MutResName}"

    return PDB

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

 The distance cutoff is incremented by 0.1 Å from the a minimum 
 value until two His and two Cys are found near each heme. The
 distance thresholds are incremented independently for the His
 residues and the Cys residues.

 A recommended minimum distance is 2.0 Å.

 What minimum distance threshold would you liek to use? """)

                print("""
 Please make sure that the correct residues are identified by, 
 for example, creating representations in VMD with the residue IDs 
 given on each line of the ResIndexing.txt file. 

 If the wrong residues are identified, the setup later with TLEaP 
 will fail because the bond definitions will be wrong. In this case, 
 please correct the residue IDs and save the changes to 
 CorrectedResIndexing.txt. When you re-run BioDC the
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

    set NumHIS {0}  
    set DistThreshHIS $DistThresh
    while {$NumHIS != 2} {
      set HIS [lsort -integer [[atomselect top "resname HIS and name NE2 and within $DistThreshHIS of resname HEC HEM and resid $ShiftHResID"] get resid]]
      set NumHIS [llength $HIS]
      set DistThreshHIS [expr {$DistThreshHIS+0.1}]
    }
    puts "$ShiftHResID | $DistThreshHIS | $NumHIS | $HIS"

    set NumCYS {0}  
    set DistThreshCYS $DistThresh
    while {$NumCYS != 2} {
      set CYS [lsort -integer [[atomselect top "resname CYS and name SG  and within $DistThreshCYS of resname HEC HEM and resid $ShiftHResID"] get resid]]
      set NumCYS [llength $CYS]
      set DistThreshCYS [expr {$DistThreshCYS+0.1}]
    }
    puts "$ShiftHResID | $DistThreshCYS | $NumCYS | $CYS"

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
            elif (CreateResIndexingMethod == "man") or (CreateResIndexingMethod == "Man") or (CreateResIndexingMethod == "manual") or (CreateResIndexingMethod == "Manual"): 
                print("""
 To create ResIndexing.txt by hand:
    Create a txt file with an editor of your choosing (e.g. 
    vi ResIndexing.txt). In this file,there must be one line for 
    each heme cofactor in your structure. Each line should contain 
    six numbers separated by a single space. The first four 
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
    the final and properly formatted PDB for use with TLEaP. This 
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
 everything from here on out will be, put politely, junk!
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
 also want to inspect the generated PDBs for the protein, each heme, 
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
 # Load parameters 
 source leaprc.constph
 source leaprc.conste
 source leaprc.water.tip3p
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
        SolvEnv = input(""" 
 Should the structure be prepared with 
  an explicit or implicit solvent (explicit/implicit)? """)

        if (SolvEnv == "explicit") or (SolvEnv == "Explicit") or (SolvEnv == "e") or (SolvEnv == "E"):
            while True:
                BoxShape = input(" Using a rectangular or an octahedral box (rec/octahed)? ")
                BufferSize = int(input(" With how much of a solvent buffer (in angstroms)? "))
            
                if (BoxShape == "rectangular") or (BoxShape == "rec"):
                    print("""

 #Solvate
 solvateBox %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                    break
                elif (BoxShape == "octahedral") or (BoxShape == "octahed") or (BoxShape == "oct"):
                    print("""

 #Solvate
 solvateOct %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                    break
                else:
                    print(" Sorry, I didn't understand your response.")

            NaCount = int(input(" And how many Na+ ions; 0 = enough for charge neutrality? "))
            ClCount = int(input(" And how many Cl- ions; 0 = enough for charge neutrality? "))

            print("""

 #Add ions
 addions %s Na+ %0d""" %(OutPrefix, NaCount), end=" ", file=open('tleap.in', 'a'))
            print("""
 addions %s Cl- %0d""" %(OutPrefix, ClCount), end=" ", file=open('tleap.in', 'a'))

            break
        elif (SolvEnv == "implicit") or (SolvEnv == "Implicit") or (SolvEnv == "i") or (SolvEnv == "I"):
            break
        else:
            print(" Sorry, I didn't understand your response.")

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

def SetUpHemeB(PDB):
    print("""
 We will write and run a TCL script that identifies the residue
 IDs of the heme groups and whcih His residues are coordinated
 to each. These residues will be re-named according to the AMBER
 parameterization for b-type hemes and a new PDB will be written.

 The TCL script also creates a BondDefinitions.txt file that will 
 be used to build the TLEaP input file. 
 
 To create the TLEaP input, we'll ask about the type of solvent, 
 unit cell shape, and numbers of cations/anions you wish to use.
    """, end=" ")
    
    OutPrefix = input(" \n Prefix for output parm/rst7 ")

    print("""
 #-------------------------------------------------------------
 #Input
 mol load pdb {}.pdb
 #-------------------------------------------------------------""".format(PDB), file=open('SetupStructure.tcl', 'w'))

    print("""
 set out [open "BondDefinitions.txt" w]
 
 set Fe [[atomselect top "name FE"] get resid]
 set NumFe [llength $Fe]
    
 for {set i 0} {$i < $NumFe} {incr i} {
   set SelFE [lindex $Fe $i]

   set NumHIS {0}; set d {2.0}
   while {$NumHIS != 2} {
     set HIS [[atomselect top "resname HIS and name NE2 and same residue as within $d of resid $SelFE"] get resid]
     set NumHIS [llength $HIS]
     set d [expr {$d+0.1}]
   }
   puts "$SelFE | $d | $NumHIS | $HIS"

   set SelFEid [atomselect top "resid $SelFE"]
   $SelFEid set resname HEB

   set HIM [lindex $HIS 0] 
   set HIMid [atomselect top "resname HIS and resid $HIM"]
   $HIMid set resname HIM

   set HIN [lindex $HIS 1] 
   set HINid [atomselect top "resname HIS and resid $HIN"]
   $HINid set resname HIN

   puts $out "bond %s.${HIM}.NE2 %s.${SelFE}.FE"
   puts $out "bond %s.${HIN}.NE2 %s.${SelFE}.FE"
 }""" %(OutPrefix, OutPrefix, OutPrefix, OutPrefix), file=open('SetupStructure.tcl', 'a'))

    print("""
 set all [atomselect top "all"]
 $all writepdb {}-{}.pdb

 exit
    """.format(OutPrefix, PDB), file=open('SetupStructure.tcl', 'a'))

   #subprocess.run("vmd -e SetupStructure.tcl > SetupStructure.log", shell=True)
    subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e SetupStructure.tcl > SetupStructure.log", shell=True)

    print("""
 VMD finished. Please check SetupStructure.log for any erros. You may 
 also want to inspect the generated PDBs for the protein,each heme, 
 and each heme propionic acid group.
    """, end=" ")

    Format = "sed -i '/OXT/a TER' "+OutPrefix+"-"+PDB+".pdb"
    subprocess.run(Format, shell=True)

    print("""
 Now Setting up TLEaP input""")

    print("""
# Reference:
# Yang, Longhua, Åge A. Skjevik, Wen-Ge Han Du, Louis Noodleman, Ross C. Walker, and Andreas W. Götz.
# Data for molecular dynamics simulations of B-type cytochrome c oxidase with the Amber force field.
# Data in brief 8 (2016): 1209-1214.
#
#
#----- leaprc for loading the Cytochrome c Oxidase of ba3 type from Thermus thermophilus. -----
#      Charges for the DNC derived using a cluster model for state 6 of the reaction cycle
#      (ref. L.Noodleman et al. Inorg. Chem., 53 (2014) 6458;
#            J.A.Fee et al. J.Am.Chem.Soc., 130 (2008) 15002.)
#

 source leaprc.protein.ff14SB
 source leaprc.water.tip3p
 loadAmberParams frcmod.ionsjc_tip3p

 addAtomTypes {
         {"FE" "Fe" "sp3"}
         { "NO" "N" "sp2" }    ## Modified by George to define NO and NP atoms as sp2 hybridised.
         { "NP" "N" "sp2" }    ## Prevents sp0 errors in leap.
 }

# Load heme b parameters:
 loadamberparams heme.frcmod
 loadoff hemeb.lib

# Load pdb and bond iron to ligating histidine nitrogens:
 %s = loadpdb %s-%s.pdb

# Connecting HEH to the protein""" %(OutPrefix, OutPrefix, PDB), file=open('tleap.in', 'w'))

    with open('BondDefinitions.txt','r') as firstfile, open('tleap.in','a') as secondfile:
        for line in firstfile:
            secondfile.write(line)

    while True:
        SolvEnv = input(""" 
 Should the structure be prepared with 
  an explicit or implicit solvent (explicit/implicit)? """)

        if (SolvEnv == "explicit") or (SolvEnv == "Explicit") or (SolvEnv == "e") or (SolvEnv == "E"):
            while True:
                BoxShape = input(" Using a rectangular or an octahedral box (rec/octahed)? ")
                BufferSize = int(input(" With how much of a solvent buffer (in angstroms)? "))
            
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
                    print(" Sorry, I didn't understand your response.")

            NaCount = int(input(" And how many Na+ ions; 0 = enough for charge neutrality? "))
            ClCount = int(input(" And how many Cl- ions; 0 = enough for charge neutrality? "))

            print("""

#Add ions
 addions %s Na+ %0d""" %(OutPrefix, NaCount), end=" ", file=open('tleap.in', 'a'))
            print("""
 addions %s Cl- %0d""" %(OutPrefix, ClCount), end=" ", file=open('tleap.in', 'a'))

            break
        elif (SolvEnv == "implicit") or (SolvEnv == "Implicit") or (SolvEnv == "i") or (SolvEnv == "I"):
            break
        else:
            print(" Sorry, I didn't understand your response.")

    print("""

#Save topology and coordinate files
 saveamberparm %s %s.prmtop %s.rst7""" %(OutPrefix, OutPrefix, OutPrefix), end=" ", file=open('tleap.in', 'a'))

    print("""

 quit""", end=" ", file=open('tleap.in', 'a'))

    print("""
 The prepared structure will now be processed with TLEaP.
    """, end=" ")
    subprocess.run("tleap -s -f tleap.in > tleap.log", shell=True)
    print("""
 TLEaP finished! 
 Please inspect the structure to make sure it is correct.
    """)
    return OutPrefix, SolvEnv 

################################################################################################################################################

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
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
 /
            """, file=open('min.in', 'w'))


            while True:
                CompChoice = input("\n Run the minimization using SANDER (S) or PMEMD (P)? ")

                if (CompChoice == "S") or (CompChoice == "s") or (CompChoice == "SANDER") or (CompChoice == "sander"):
                    subprocess.run("                 sander -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
                    print(" Minimization finished!")
                    break
                elif (CompChoice == "P") or (CompChoice == "p") or (CompChoice == "PMEMD") or (CompChoice == "pmemd"):
                    NProc = input(" parallelize the minimization over how many CPUs? ")
                    subprocess.run("mpirun -np"+NProc+"pmemd.MPI -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
                    print(" Running minimization ...")
                    break
                else:
                    print(" Sorry, I didn't understand your response. Please try again.")
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
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
/
            """, file=open('min.in', 'w'))

            while True:
                CompChoice = input("\n Run the minimization using SANDER (S) or PMEMD (P)? ")

                if (CompChoice == "S") or (CompChoice == "s") or (CompChoice == "SANDER") or (CompChoice == "sander"):
                    subprocess.run("                 sander -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
                    print(" Minimization finished!")
                    break
                elif (CompChoice == "P") or (CompChoice == "p") or (CompChoice == "PMEMD") or (CompChoice == "pmemd"):
                    NProc = input(" parallelize the minimization over how many CPUs? ")
                    subprocess.run("mpirun -np"+NProc+"pmemd.MPI -O -i min.in -o min.out -p "+Output+".prmtop -c "+Output+".rst7 -inf min.mdinfo -r min.rst7 -ref "+Output+".rst7", shell=True)
                    print(" Running minimization ...")
                    break
                else:
                    print(" Sorry, I didn't understand your response. Please try again.")
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
    
    subprocess.run("ambpdb -p "+Output+"_reord.prmtop -c "+Output+"_reord.rst7 > min.pdb", shell=True)

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
                #subprocess.run("vmd -e SASACalc.tcl > SASACalc.log", shell=True)
                subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e SASACalc.tcl > SASACalc.log", shell=True)

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

################################################################################################################################################

def LambdaFromSASA_HemeB(Output):

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

        HisID = []; HebID = []; HEM = []
        with open("BondDefinitions.txt") as f:
            for line in f:
                HisID.append(int(line.strip().split('.')[1]))
                HebID.append(int(line.strip().split('.')[3]))
 
        for i in range(0, len(HisID), 2):
            HEM.append([HebID[i], HisID[i], HisID[i+1]])

        for idx in range(len(HEM)-1):
           #print(HEM[idx][0], HEM[idx][1], HEM[idx][2], HEM[idx+1][0], HEM[idx+1][1], HEM[idx+1][2])

            print("""
 mol new min.pdb

 set HEB1 %0d; set HIM1 %0d; set HIN1 %0d
 set HEB2 %0d; set HIM2 %0d; set HIN2 %0d""" %(HEM[idx][0], HEM[idx][1], HEM[idx][2], HEM[idx+1][0], HEM[idx+1][1], HEM[idx+1][2]), file=open('SASACalc.tcl', 'w'))  

            print("""
 set output [open "${HEB1},${HEB2}_SASAanalysis.dat" a]

 set allsel     [atomselect top "all and not water and not ions"]

 set donor      [atomselect top "(resname HEB and resid $HEB1) or (resname HIM and resid $HIM1) or (resname HIN and resid $HIN1)"]
 set Dfe        [atomselect top "(resname HEB and resid $HEB1 and name FE)"]
 set DfeIDX     [$Dfe get index]
 set DfeID      [$Dfe get resid]

 set acceptor   [atomselect top "(resname HEB and resid $HEB2) or (resname HIM and resid $HIM2) or (resname HIN and resid $HIN2)"]
 set Afe        [atomselect top "(resname HEB and resid $HEB2 and name FE)"]
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
            
            if (os.path.isfile("%0d,%0d_SASAanalysis.dat" %(HEM[idx][0], HEM[idx+1][0])) == True):
                print(""" Found %0d,%0d_SASAanalysis.dat from a prior execution.
 This prior output will be used for the analysis.
  """ %(HEM[idx][0], HEM[idx+1][0])) 

            if (os.path.isfile("%0d,%0d_SASAanalysis.dat" %(HEM[idx][0], HEM[idx+1][0])) == False):
                print(" Now using VMD to compute SASA Donor = %0d & Acceptor = %0d..." %(HEM[idx][0], HEM[idx+1][0]))
                #subprocess.run("vmd -e SASACalc.tcl > SASACalc.log", shell=True)
                subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e SASACalc.tcl > SASACalc.log", shell=True)

        print(" Computing Reorganization Energy from Solvent Accessibility...")
        alpha = 5.18; beta = 0.016;
        Rd=4.6; Ra=4.6;
        Eopt=1.84 

        Dsasa = [0]*(len(HEM)-1)
        Asasa = [0]*(len(HEM)-1)
        TotalSASA = [0]*(len(HEM)-1)
        Es = [0]*(len(HEM)-1)
        M = [0]*(len(HEM)-1)
        Rda = [0]*(len(HEM)-1)
        R = [0]*(len(HEM)-1)
        Lambda = [0]*(len(HEM)-1)
        for idx in range(len(HEM)-1):
            with open(str(HEM[idx][0])+","+str(HEM[idx+1][0])+"_SASAanalysis.dat") as fp:
                #print(str(HEM[idx][0])+","+str(HEM[idx+1][0])+"_SASAanalysis.dat")
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

        for idx in range(len(HEM)-1):
            if (idx == 0):
                print(""" 
 HEM-%0d -> HEM-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" %(HEM[idx][0], HEM[idx+1][0], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'w'))

            if (idx != 0):
                print(""" 
 HEM-%0d -> HEM-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" %(HEM[idx][0], HEM[idx+1][0], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'a'))

        print(" Done!")

    return Lambda

################################################################################################################################################

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

    while True:
        CompChoice = input(" Do you wish to run any needed computaitons in serial or parallel (s/p)? ")

        if (CompChoice == "s") or (CompChoice == "p") or (CompChoice == "serial") or (CompChoice == "parallel"):
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

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
                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                        print(" Running PBSA calculation for %s. HEH-%0d ..." %(state, HEH[idx]))
                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(RedoxState[0])+str(HEH[idx])+" -p "+str(RedoxState[0])+str(HEH[idx])+".prmtop -c StrucForPBSA.rst7", shell=True)
                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(RedoxState[0])+str(HEH[idx])+" -p "+str(RedoxState[0])+str(HEH[idx])+".prmtop -c StrucForPBSA.rst7"
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

################################################################################################################################################

def PairedChargeAssignment(Output, HEH, i, j, k, l):
    if (k == "o"): 
        kstate = 0
    if (k == "r"):
        kstate = 1

    if (l == "o"): 
        lstate = 0
    if (l == "r"):
        lstate = 1

    print(k, kstate, l, lstate)

    if (os.path.isfile(f"{k}{HEH[i]}-{l}{HEH[j]}.inp") == True):
        print(f"\n Generating topology for {k}{HEH[i]}-{l}{HEH[j]}  ...")
        subprocess.run(f"parmed -i {k}{HEH[i]}-{l}{HEH[j]}.inp > {k}{HEH[i]}-{l}{HEH[j]}.log", shell=True)
    if (os.path.isfile(f"{k}{HEH[i]}-{l}{HEH[j]}.inp") == False):
        print(f"""
 parm StrucForPBSA.{Output}.prmtop
 changeRedoxState :{HEH[i]} {kstate}
 changeRedoxState :{HEH[j]} {lstate}
 outparm {k}{HEH[i]}-{l}{HEH[j]}.prmtop
 quit
        """, file=open(f"{k}{HEH[i]}-{l}{HEH[j]}.inp", 'w'))

        print(f"\n Generating topology for {k}{HEH[i]}-{l}{HEH[j]}  ...")
        subprocess.run(f"parmed -i {k}{HEH[i]}-{l}{HEH[j]}.inp > {k}{HEH[i]}-{l}{HEH[j]}.log", shell=True)

################################################################################################################################################

################################################################################################################################################

def HemeHemeInt(Output, Solv):

    if (os.path.isfile("pbsa.key") == True):
        print(" Found pbsa.key")
    elif (os.path.isfile("pbsa.key") == False):
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
        if (os.path.isfile("GenerateCoordForPBSA.in") == True):
            if (os.path.isfile("StrucForPBSA."+Output+".prmtop") == True) and (os.path.isfile("StrucForPBSA.rst7") == True):
                print(f" Found StrucForPBSA.{Output}.prmtop & StrucForPBSA.rst7")
        elif (os.path.isfile("GenerateCoordForPBSA.in") == False):
            print("""
 parm    %s.prmtop
 trajin  min.rst7
 strip   :WAT,Na+,Cl- outprefix StrucForPBSA
 trajout StrucForPBSA.rst7
 run
            """ %(Output), file=open("GenerateCoordForPBSA.in", 'w'))
            subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)
    elif (Solv == "implicit") or (Solv == "Implicit") or (Solv == "i") or (Solv == "I"):
        if (os.path.isfile("GenerateCoordForPBSA.in") == True):
            if (os.path.isfile("StrucForPBSA."+Output+".prmtop") == True) and (os.path.isfile("StrucForPBSA.rst7") == True):
                print(f" Found StrucForPBSA.{Output}.prmtop & StrucForPBSA.rst7")
        elif (os.path.isfile("GenerateCoordForPBSA.in") == False):
            print("""
 parm    %s.prmtop
 trajin  min.rst7
 parmwrite out StrucForPBSA.%s.prmtop
 trajout StrucForPBSA.rst7
 run
            """ %(Output, Output), file=open("GenerateCoordForPBSA.in", 'w'))
            subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)

    while True:
        CompChoice = input(" Do you wish to run any needed computaitons in serial or parallel (s/p)? ")

        if (CompChoice == "s") or (CompChoice == "p") or (CompChoice == "serial") or (CompChoice == "parallel"):
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEH = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEH[idx] = int(line.strip().split(" ")[1])
            idx += 1

    M = [[0 for column in range(len(HEH))] for row in range(len(HEH))]
    N = [[0 for column in range(4)] for row in range(4)]

    print("State Energies", file=open('CheckValues.txt', 'w'))
    for i in range(0, len(HEH)):
        for j in range(0, len(HEH)):
            if (j == i):
                pass
            if (j != i):
                idxc = 0
                command = ['']*(4)
                print(f"""
 ================================================
 Analyzing pair {HEH[i]} {HEH[j]}...
 ================================================""")

                for k in ("o", "r"):
                    for l in ("o", "r"):

                        if (k == "o") and (l == "o"):
                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == False):
                                PairedChargeAssignment(Output, HEH, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEH[i]}-{l}{HEH[j]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "r") and (l == "o"):
                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == False):
                                PairedChargeAssignment(Output, HEH, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEH[i]}-{l}{HEH[j]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "o") and (l == "r"):
                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == False):
                                PairedChargeAssignment(Output, HEH, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEH[i]}-{l}{HEH[j]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "r") and (l == "r"):
                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == False):
                                PairedChargeAssignment(Output, HEH, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEH[i]}-{l}{HEH[j]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEH[i]}-{l}{HEH[j]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+" -p "+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1
            
                if (idxc != 0):
                    commandrev = list(filter(None, command))
                    print("\n Submitting "+str(len(commandrev))+" PBSA calculations in parallel")
                    procs = [ subprocess.Popen(i, shell=True) for i in commandrev ]

                    for p in procs:
                        p.wait()
                        print("  Finished: "+str(p))

                chk = 4
                if (os.path.isfile("pbsa_o"+str(HEH[i])+"-o"+str(HEH[j])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_o{HEH[i]}-o{HEH[j]} is missing.""")
                if (os.path.isfile("pbsa_o"+str(HEH[i])+"-r"+str(HEH[j])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_o{HEH[i]}-r{HEH[j]} is missing.""")
                if (os.path.isfile("pbsa_r"+str(HEH[i])+"-o"+str(HEH[j])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_r{HEH[i]}-o{HEH[j]} is missing.""")
                if (os.path.isfile("pbsa_r"+str(HEH[i])+"-r"+str(HEH[j])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_r{HEH[i]}-r{HEH[j]} is missing.""")

                if (chk == 4):
                   #print(" All four files found") 
                     
                    for k in ("o", "r"):
                        for l in ("o", "r"):

                            if (k == "o") and (l == "o"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[0][0] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[0][0], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                            if (k == "o") and (l == "r"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[0][1] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[0][1], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                            if (k == "r") and (l == "o"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[1][0] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[1][0], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1
                                    
                            if (k == "r") and (l == "r"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEH[i])+"-"+str(l)+str(HEH[j]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[1][1] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[1][1], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                M[i][i] = round(((N[0][1] - N[1][1]))*1000, 3)
                M[i][j] = round(((N[0][0] - N[1][0]) - (N[0][1] - N[1][1]))*1000, 3)
                M[j][i] = round(((N[0][0] - N[0][1]) - (N[1][0] - N[1][1]))*1000, 3)
                M[j][j] = round(((N[1][0] - N[1][1]))*1000, 3)
    print(*M, sep='\n')
    print(*M, sep='\n', file=open('CheckValues.txt', 'a'))

################################################################################################################################################
   
################################################################################################################################################
def PairedChargeAssignment_HemeB(Output, HEM, i, j, k, l):
    if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp") == True):

        print(f"\n Generating topology for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
        subprocess.run("parmed -i "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp > "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".log", shell=True)

    elif (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp") == False):
        print(f"""
parm StrucForPBSA.{Output}.prmtop""", file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'w'))

        if (k == "o"):
            print(f"""
netcharge :{HEM[i][0]},{HEM[i][1]},{HEM[i][2]}
change CHARGE :{HEM[i][0]}&@FE         0.6083
change CHARGE :{HEM[i][0]}&@NA         0.1780
change CHARGE :{HEM[i][0]}&@C1A       -0.2990
change CHARGE :{HEM[i][0]}&@C2A        0.2176
change CHARGE :{HEM[i][0]}&@C3A       -0.0578
change CHARGE :{HEM[i][0]}&@CMA       -0.0653
change CHARGE :{HEM[i][0]}&@HMA1       0.0838
change CHARGE :{HEM[i][0]}&@HMA2       0.0416
change CHARGE :{HEM[i][0]}&@HMA3       0.1007
change CHARGE :{HEM[i][0]}&@C4A       -0.2008
change CHARGE :{HEM[i][0]}&@CHB       -0.1385
change CHARGE :{HEM[i][0]}&@HHB        0.3521
change CHARGE :{HEM[i][0]}&@C1B       -0.1639
change CHARGE :{HEM[i][0]}&@NB         0.1343
change CHARGE :{HEM[i][0]}&@C2B       -0.1095
change CHARGE :{HEM[i][0]}&@CMB        0.0161
change CHARGE :{HEM[i][0]}&@HMB1       0.0095
change CHARGE :{HEM[i][0]}&@HMB2       0.0688
change CHARGE :{HEM[i][0]}&@HMB3       0.0593
change CHARGE :{HEM[i][0]}&@C3B        0.1107
change CHARGE :{HEM[i][0]}&@CAB       -0.2074
change CHARGE :{HEM[i][0]}&@HAB        0.1943
change CHARGE :{HEM[i][0]}&@CBB       -0.4054
change CHARGE :{HEM[i][0]}&@HBB1       0.2192
change CHARGE :{HEM[i][0]}&@HBB2       0.1998
change CHARGE :{HEM[i][0]}&@C4B       -0.0043
change CHARGE :{HEM[i][0]}&@CHC       -0.2657
change CHARGE :{HEM[i][0]}&@HHC        0.2630
change CHARGE :{HEM[i][0]}&@C1C        0.1417
change CHARGE :{HEM[i][0]}&@NC         0.0897
change CHARGE :{HEM[i][0]}&@C2C       -0.0225
change CHARGE :{HEM[i][0]}&@CMC       -0.1841
change CHARGE :{HEM[i][0]}&@HMC1       0.1221
change CHARGE :{HEM[i][0]}&@HMC2       0.1268
change CHARGE :{HEM[i][0]}&@HMC3       0.0719
change CHARGE :{HEM[i][0]}&@C3C        0.0485
change CHARGE :{HEM[i][0]}&@CAC        0.0012
change CHARGE :{HEM[i][0]}&@HAC        0.1508
change CHARGE :{HEM[i][0]}&@CBC       -0.5097
change CHARGE :{HEM[i][0]}&@HBC1       0.2648
change CHARGE :{HEM[i][0]}&@HBC2       0.2618
change CHARGE :{HEM[i][0]}&@C4C       -0.2418
change CHARGE :{HEM[i][0]}&@CHD       -0.1909
change CHARGE :{HEM[i][0]}&@HHD        0.3192
change CHARGE :{HEM[i][0]}&@C1D       -0.0856
change CHARGE :{HEM[i][0]}&@ND         0.1145
change CHARGE :{HEM[i][0]}&@C2D       -0.0592
change CHARGE :{HEM[i][0]}&@CMD       -0.0045
change CHARGE :{HEM[i][0]}&@HMD1       0.0235
change CHARGE :{HEM[i][0]}&@HMD2       0.0942
change CHARGE :{HEM[i][0]}&@HMD3       0.0215
change CHARGE :{HEM[i][0]}&@C3D       -0.1936
change CHARGE :{HEM[i][0]}&@C4D        0.0685
change CHARGE :{HEM[i][0]}&@CHA       -0.1387
change CHARGE :{HEM[i][0]}&@HHA        0.3203
change CHARGE :{HEM[i][1]}&@CB         0.0590
change CHARGE :{HEM[i][1]}&@HB2        0.0431
change CHARGE :{HEM[i][1]}&@HB3        0.0431
change CHARGE :{HEM[i][1]}&@CG         0.0368
change CHARGE :{HEM[i][1]}&@ND1       -0.1361
change CHARGE :{HEM[i][1]}&@HD1        0.4020
change CHARGE :{HEM[i][1]}&@CE1       -0.1891
change CHARGE :{HEM[i][1]}&@HE1        0.2281
change CHARGE :{HEM[i][1]}&@NE2       -0.0977
change CHARGE :{HEM[i][1]}&@CD2       -0.2862
change CHARGE :{HEM[i][1]}&@HD2        0.2051
change CHARGE :{HEM[i][2]}&@CB        -0.0209
change CHARGE :{HEM[i][2]}&@HB2        0.0667
change CHARGE :{HEM[i][2]}&@HB3        0.0667
change CHARGE :{HEM[i][2]}&@CG        -0.0517
change CHARGE :{HEM[i][2]}&@ND1       -0.0005
change CHARGE :{HEM[i][2]}&@HD1        0.3810
change CHARGE :{HEM[i][2]}&@CE1       -0.3022
change CHARGE :{HEM[i][2]}&@HE1        0.2830
change CHARGE :{HEM[i][2]}&@NE2       -0.1741
change CHARGE :{HEM[i][2]}&@CD2       -0.2124
change CHARGE :{HEM[i][2]}&@HD2        0.2337
change CHARGE :{HEM[i][0]}&@CAA        0.1105
change CHARGE :{HEM[i][0]}&@HAA1       0.0269
change CHARGE :{HEM[i][0]}&@HAA2      -0.0203
change CHARGE :{HEM[i][0]}&@CBA       -0.0733
change CHARGE :{HEM[i][0]}&@HBA1       0.0006
change CHARGE :{HEM[i][0]}&@HBA2       0.0585
change CHARGE :{HEM[i][0]}&@CGA        0.5918
change CHARGE :{HEM[i][0]}&@O1A       -0.6537
change CHARGE :{HEM[i][0]}&@O2A       -0.6853
change CHARGE :{HEM[i][0]}&@CAD       -0.0291
change CHARGE :{HEM[i][0]}&@HAD1       0.0078
change CHARGE :{HEM[i][0]}&@HAD2       0.0658
change CHARGE :{HEM[i][0]}&@CBD        0.2827
change CHARGE :{HEM[i][0]}&@HBD1      -0.0293
change CHARGE :{HEM[i][0]}&@HBD2      -0.0436
change CHARGE :{HEM[i][0]}&@CGD        0.4783
change CHARGE :{HEM[i][0]}&@O1D       -0.6025
change CHARGE :{HEM[i][0]}&@O2D       -0.6098
netcharge :{HEM[i][0]},{HEM[i][1]},{HEM[i][2]}""", file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'a'))

        if (k == "r"):
            print(f"""
netcharge :{HEM[i][0]},{HEM[i][1]},{HEM[i][2]}
change CHARGE :{HEM[i][0]}&@FE         0.4223
change CHARGE :{HEM[i][0]}&@NA         0.1973
change CHARGE :{HEM[i][0]}&@C1A       -0.3489
change CHARGE :{HEM[i][0]}&@C2A        0.2187
change CHARGE :{HEM[i][0]}&@C3A       -0.0877
change CHARGE :{HEM[i][0]}&@CMA       -0.0598
change CHARGE :{HEM[i][0]}&@HMA1       0.0703
change CHARGE :{HEM[i][0]}&@HMA2       0.0281
change CHARGE :{HEM[i][0]}&@HMA3       0.0872
change CHARGE :{HEM[i][0]}&@C4A       -0.2107
change CHARGE :{HEM[i][0]}&@CHB       -0.1614
change CHARGE :{HEM[i][0]}&@HHB        0.3441
change CHARGE :{HEM[i][0]}&@C1B       -0.1858
change CHARGE :{HEM[i][0]}&@NB         0.1413
change CHARGE :{HEM[i][0]}&@C2B       -0.0884   
change CHARGE :{HEM[i][0]}&@CMB        0.0466
change CHARGE :{HEM[i][0]}&@HMB1      -0.0146
change CHARGE :{HEM[i][0]}&@HMB2       0.0447
change CHARGE :{HEM[i][0]}&@HMB3       0.0352
change CHARGE :{HEM[i][0]}&@C3B        0.0348
change CHARGE :{HEM[i][0]}&@CAB       -0.2072
change CHARGE :{HEM[i][0]}&@HAB        0.1898
change CHARGE :{HEM[i][0]}&@CBB       -0.4301
change CHARGE :{HEM[i][0]}&@HBB1       0.2165
change CHARGE :{HEM[i][0]}&@HBB2       0.1971
change CHARGE :{HEM[i][0]}&@C4B       -0.0142
change CHARGE :{HEM[i][0]}&@CHC       -0.2576
change CHARGE :{HEM[i][0]}&@HHC        0.2520
change CHARGE :{HEM[i][0]}&@C1C        0.0528
change CHARGE :{HEM[i][0]}&@NC         0.1307
change CHARGE :{HEM[i][0]}&@C2C        0.0116
change CHARGE :{HEM[i][0]}&@CMC       -0.1501
change CHARGE :{HEM[i][0]}&@HMC1       0.0991
change CHARGE :{HEM[i][0]}&@HMC2       0.1038
change CHARGE :{HEM[i][0]}&@HMC3       0.0489
change CHARGE :{HEM[i][0]}&@C3C        0.0046
change CHARGE :{HEM[i][0]}&@CAC       -0.0156
change CHARGE :{HEM[i][0]}&@HAC        0.1293
change CHARGE :{HEM[i][0]}&@CBC       -0.5384
change CHARGE :{HEM[i][0]}&@HBC1       0.2611
change CHARGE :{HEM[i][0]}&@HBC2       0.2581
change CHARGE :{HEM[i][0]}&@C4C       -0.2427
change CHARGE :{HEM[i][0]}&@CHD       -0.2138
change CHARGE :{HEM[i][0]}&@HHD        0.3122
change CHARGE :{HEM[i][0]}&@C1D       -0.0975
change CHARGE :{HEM[i][0]}&@ND         0.1255
change CHARGE :{HEM[i][0]}&@C2D       -0.0861
change CHARGE :{HEM[i][0]}&@CMD        0.0005
change CHARGE :{HEM[i][0]}&@HMD1       0.0098
change CHARGE :{HEM[i][0]}&@HMD2       0.0805
change CHARGE :{HEM[i][0]}&@HMD3       0.0078
change CHARGE :{HEM[i][0]}&@C3D       -0.1965
change CHARGE :{HEM[i][0]}&@C4D        0.0236
change CHARGE :{HEM[i][0]}&@CHA       -0.1216
change CHARGE :{HEM[i][0]}&@HHA        0.3063
change CHARGE :{HEM[i][1]}&@CB         0.0426
change CHARGE :{HEM[i][1]}&@HB2        0.0349
change CHARGE :{HEM[i][1]}&@HB3        0.0349
change CHARGE :{HEM[i][1]}&@CG         0.0408
change CHARGE :{HEM[i][1]}&@ND1       -0.1881
change CHARGE :{HEM[i][1]}&@HD1        0.3920
change CHARGE :{HEM[i][1]}&@CE1       -0.1901
change CHARGE :{HEM[i][1]}&@HE1        0.2191
change CHARGE :{HEM[i][1]}&@NE2       -0.0847
change CHARGE :{HEM[i][1]}&@CD2       -0.3012
change CHARGE :{HEM[i][1]}&@HD2        0.2011
change CHARGE :{HEM[i][2]}&@CB        -0.0373
change CHARGE :{HEM[i][2]}&@HB2        0.0585
change CHARGE :{HEM[i][2]}&@HB3        0.0585
change CHARGE :{HEM[i][2]}&@CG        -0.0477
change CHARGE :{HEM[i][2]}&@ND1       -0.0525
change CHARGE :{HEM[i][2]}&@HD1        0.3710
change CHARGE :{HEM[i][2]}&@CE1       -0.3032
change CHARGE :{HEM[i][2]}&@HE1        0.2740
change CHARGE :{HEM[i][2]}&@NE2       -0.1611
change CHARGE :{HEM[i][2]}&@CD2       -0.2274
change CHARGE :{HEM[i][2]}&@HD2        0.2297
change CHARGE :{HEM[i][0]}&@CAA        0.1105
change CHARGE :{HEM[i][0]}&@HAA1       0.0269
change CHARGE :{HEM[i][0]}&@HAA2      -0.0203
change CHARGE :{HEM[i][0]}&@CBA       -0.0733
change CHARGE :{HEM[i][0]}&@HBA1       0.0006
change CHARGE :{HEM[i][0]}&@HBA2       0.0585
change CHARGE :{HEM[i][0]}&@CGA        0.5918
change CHARGE :{HEM[i][0]}&@O1A       -0.6537
change CHARGE :{HEM[i][0]}&@O2A       -0.6853
change CHARGE :{HEM[i][0]}&@CAD       -0.0291
change CHARGE :{HEM[i][0]}&@HAD1       0.0078
change CHARGE :{HEM[i][0]}&@HAD2       0.0658
change CHARGE :{HEM[i][0]}&@CBD        0.2827
change CHARGE :{HEM[i][0]}&@HBD1      -0.0293
change CHARGE :{HEM[i][0]}&@HBD2      -0.0436
change CHARGE :{HEM[i][0]}&@CGD        0.4783
change CHARGE :{HEM[i][0]}&@O1D       -0.6025
change CHARGE :{HEM[i][0]}&@O2D       -0.6098
netcharge :{HEM[j][0]},{HEM[j][1]},{HEM[j][2]}""", file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'a'))

        if (l == "o"):
            print(f"""
netcharge :{HEM[j][0]},{HEM[j][1]},{HEM[j][2]}
change CHARGE :{HEM[j][0]}&@FE         0.6083
change CHARGE :{HEM[j][0]}&@NA         0.1780
change CHARGE :{HEM[j][0]}&@C1A       -0.2990
change CHARGE :{HEM[j][0]}&@C2A        0.2176
change CHARGE :{HEM[j][0]}&@C3A       -0.0578
change CHARGE :{HEM[j][0]}&@CMA       -0.0653
change CHARGE :{HEM[j][0]}&@HMA1       0.0838
change CHARGE :{HEM[j][0]}&@HMA2       0.0416
change CHARGE :{HEM[j][0]}&@HMA3       0.1007
change CHARGE :{HEM[j][0]}&@C4A       -0.2008
change CHARGE :{HEM[j][0]}&@CHB       -0.1385
change CHARGE :{HEM[j][0]}&@HHB        0.3521
change CHARGE :{HEM[j][0]}&@C1B       -0.1639
change CHARGE :{HEM[j][0]}&@NB         0.1343
change CHARGE :{HEM[j][0]}&@C2B       -0.1095
change CHARGE :{HEM[j][0]}&@CMB        0.0161
change CHARGE :{HEM[j][0]}&@HMB1       0.0095
change CHARGE :{HEM[j][0]}&@HMB2       0.0688
change CHARGE :{HEM[j][0]}&@HMB3       0.0593
change CHARGE :{HEM[j][0]}&@C3B        0.1107
change CHARGE :{HEM[j][0]}&@CAB       -0.2074
change CHARGE :{HEM[j][0]}&@HAB        0.1943
change CHARGE :{HEM[j][0]}&@CBB       -0.4054
change CHARGE :{HEM[j][0]}&@HBB1       0.2192
change CHARGE :{HEM[j][0]}&@HBB2       0.1998
change CHARGE :{HEM[j][0]}&@C4B       -0.0043
change CHARGE :{HEM[j][0]}&@CHC       -0.2657
change CHARGE :{HEM[j][0]}&@HHC        0.2630
change CHARGE :{HEM[j][0]}&@C1C        0.1417
change CHARGE :{HEM[j][0]}&@NC         0.0897
change CHARGE :{HEM[j][0]}&@C2C       -0.0225
change CHARGE :{HEM[j][0]}&@CMC       -0.1841
change CHARGE :{HEM[j][0]}&@HMC1       0.1221
change CHARGE :{HEM[j][0]}&@HMC2       0.1268
change CHARGE :{HEM[j][0]}&@HMC3       0.0719
change CHARGE :{HEM[j][0]}&@C3C        0.0485
change CHARGE :{HEM[j][0]}&@CAC        0.0012
change CHARGE :{HEM[j][0]}&@HAC        0.1508
change CHARGE :{HEM[j][0]}&@CBC       -0.5097
change CHARGE :{HEM[j][0]}&@HBC1       0.2648
change CHARGE :{HEM[j][0]}&@HBC2       0.2618
change CHARGE :{HEM[j][0]}&@C4C       -0.2418
change CHARGE :{HEM[j][0]}&@CHD       -0.1909
change CHARGE :{HEM[j][0]}&@HHD        0.3192
change CHARGE :{HEM[j][0]}&@C1D       -0.0856
change CHARGE :{HEM[j][0]}&@ND         0.1145
change CHARGE :{HEM[j][0]}&@C2D       -0.0592
change CHARGE :{HEM[j][0]}&@CMD       -0.0045
change CHARGE :{HEM[j][0]}&@HMD1       0.0235
change CHARGE :{HEM[j][0]}&@HMD2       0.0942
change CHARGE :{HEM[j][0]}&@HMD3       0.0215
change CHARGE :{HEM[j][0]}&@C3D       -0.1936
change CHARGE :{HEM[j][0]}&@C4D        0.0685
change CHARGE :{HEM[j][0]}&@CHA       -0.1387
change CHARGE :{HEM[j][0]}&@HHA        0.3203
change CHARGE :{HEM[j][1]}&@CB         0.0590
change CHARGE :{HEM[j][1]}&@HB2        0.0431
change CHARGE :{HEM[j][1]}&@HB3        0.0431
change CHARGE :{HEM[j][1]}&@CG         0.0368
change CHARGE :{HEM[j][1]}&@ND1       -0.1361
change CHARGE :{HEM[j][1]}&@HD1        0.4020
change CHARGE :{HEM[j][1]}&@CE1       -0.1891
change CHARGE :{HEM[j][1]}&@HE1        0.2281
change CHARGE :{HEM[j][1]}&@NE2       -0.0977
change CHARGE :{HEM[j][1]}&@CD2       -0.2862
change CHARGE :{HEM[j][1]}&@HD2        0.2051
change CHARGE :{HEM[j][2]}&@CB        -0.0209
change CHARGE :{HEM[j][2]}&@HB2        0.0667
change CHARGE :{HEM[j][2]}&@HB3        0.0667
change CHARGE :{HEM[j][2]}&@CG        -0.0517
change CHARGE :{HEM[j][2]}&@ND1       -0.0005
change CHARGE :{HEM[j][2]}&@HD1        0.3810
change CHARGE :{HEM[j][2]}&@CE1       -0.3022
change CHARGE :{HEM[j][2]}&@HE1        0.2830
change CHARGE :{HEM[j][2]}&@NE2       -0.1741
change CHARGE :{HEM[j][2]}&@CD2       -0.2124
change CHARGE :{HEM[j][2]}&@HD2        0.2337
change CHARGE :{HEM[j][0]}&@CAA        0.1105
change CHARGE :{HEM[j][0]}&@HAA1       0.0269
change CHARGE :{HEM[j][0]}&@HAA2      -0.0203
change CHARGE :{HEM[j][0]}&@CBA       -0.0733
change CHARGE :{HEM[j][0]}&@HBA1       0.0006
change CHARGE :{HEM[j][0]}&@HBA2       0.0585
change CHARGE :{HEM[j][0]}&@CGA        0.5918
change CHARGE :{HEM[j][0]}&@O1A       -0.6537
change CHARGE :{HEM[j][0]}&@O2A       -0.6853
change CHARGE :{HEM[j][0]}&@CAD       -0.0291
change CHARGE :{HEM[j][0]}&@HAD1       0.0078
change CHARGE :{HEM[j][0]}&@HAD2       0.0658
change CHARGE :{HEM[j][0]}&@CBD        0.2827
change CHARGE :{HEM[j][0]}&@HBD1      -0.0293
change CHARGE :{HEM[j][0]}&@HBD2      -0.0436
change CHARGE :{HEM[j][0]}&@CGD        0.4783
change CHARGE :{HEM[j][0]}&@O1D       -0.6025
change CHARGE :{HEM[j][0]}&@O2D       -0.6098
netcharge :{HEM[j][0]},{HEM[j][1]},{HEM[j][2]}""", file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'a'))
 
        if (l == "r"):
            print(f"""
netcharge :{HEM[j][0]},{HEM[j][1]},{HEM[j][2]}
change CHARGE :{HEM[j][0]}&@FE         0.4223
change CHARGE :{HEM[j][0]}&@NA         0.1973
change CHARGE :{HEM[j][0]}&@C1A       -0.3489
change CHARGE :{HEM[j][0]}&@C2A        0.2187
change CHARGE :{HEM[j][0]}&@C3A       -0.0877
change CHARGE :{HEM[j][0]}&@CMA       -0.0598
change CHARGE :{HEM[j][0]}&@HMA1       0.0703
change CHARGE :{HEM[j][0]}&@HMA2       0.0281
change CHARGE :{HEM[j][0]}&@HMA3       0.0872
change CHARGE :{HEM[j][0]}&@C4A       -0.2107
change CHARGE :{HEM[j][0]}&@CHB       -0.1614
change CHARGE :{HEM[j][0]}&@HHB        0.3441
change CHARGE :{HEM[j][0]}&@C1B       -0.1858
change CHARGE :{HEM[j][0]}&@NB         0.1413
change CHARGE :{HEM[j][0]}&@C2B       -0.0884   
change CHARGE :{HEM[j][0]}&@CMB        0.0466
change CHARGE :{HEM[j][0]}&@HMB1      -0.0146
change CHARGE :{HEM[j][0]}&@HMB2       0.0447
change CHARGE :{HEM[j][0]}&@HMB3       0.0352
change CHARGE :{HEM[j][0]}&@C3B        0.0348
change CHARGE :{HEM[j][0]}&@CAB       -0.2072
change CHARGE :{HEM[j][0]}&@HAB        0.1898
change CHARGE :{HEM[j][0]}&@CBB       -0.4301
change CHARGE :{HEM[j][0]}&@HBB1       0.2165
change CHARGE :{HEM[j][0]}&@HBB2       0.1971
change CHARGE :{HEM[j][0]}&@C4B       -0.0142
change CHARGE :{HEM[j][0]}&@CHC       -0.2576
change CHARGE :{HEM[j][0]}&@HHC        0.2520
change CHARGE :{HEM[j][0]}&@C1C        0.0528
change CHARGE :{HEM[j][0]}&@NC         0.1307
change CHARGE :{HEM[j][0]}&@C2C        0.0116
change CHARGE :{HEM[j][0]}&@CMC       -0.1501
change CHARGE :{HEM[j][0]}&@HMC1       0.0991
change CHARGE :{HEM[j][0]}&@HMC2       0.1038
change CHARGE :{HEM[j][0]}&@HMC3       0.0489
change CHARGE :{HEM[j][0]}&@C3C        0.0046
change CHARGE :{HEM[j][0]}&@CAC       -0.0156
change CHARGE :{HEM[j][0]}&@HAC        0.1293
change CHARGE :{HEM[j][0]}&@CBC       -0.5384
change CHARGE :{HEM[j][0]}&@HBC1       0.2611
change CHARGE :{HEM[j][0]}&@HBC2       0.2581
change CHARGE :{HEM[j][0]}&@C4C       -0.2427
change CHARGE :{HEM[j][0]}&@CHD       -0.2138
change CHARGE :{HEM[j][0]}&@HHD        0.3122
change CHARGE :{HEM[j][0]}&@C1D       -0.0975
change CHARGE :{HEM[j][0]}&@ND         0.1255
change CHARGE :{HEM[j][0]}&@C2D       -0.0861
change CHARGE :{HEM[j][0]}&@CMD        0.0005
change CHARGE :{HEM[j][0]}&@HMD1       0.0098
change CHARGE :{HEM[j][0]}&@HMD2       0.0805
change CHARGE :{HEM[j][0]}&@HMD3       0.0078
change CHARGE :{HEM[j][0]}&@C3D       -0.1965
change CHARGE :{HEM[j][0]}&@C4D        0.0236
change CHARGE :{HEM[j][0]}&@CHA       -0.1216
change CHARGE :{HEM[j][0]}&@HHA        0.3063
change CHARGE :{HEM[j][1]}&@CB         0.0426
change CHARGE :{HEM[j][1]}&@HB2        0.0349
change CHARGE :{HEM[j][1]}&@HB3        0.0349
change CHARGE :{HEM[j][1]}&@CG         0.0408
change CHARGE :{HEM[j][1]}&@ND1       -0.1881
change CHARGE :{HEM[j][1]}&@HD1        0.3920
change CHARGE :{HEM[j][1]}&@CE1       -0.1901
change CHARGE :{HEM[j][1]}&@HE1        0.2191
change CHARGE :{HEM[j][1]}&@NE2       -0.0847
change CHARGE :{HEM[j][1]}&@CD2       -0.3012
change CHARGE :{HEM[j][1]}&@HD2        0.2011
change CHARGE :{HEM[j][2]}&@CB        -0.0373
change CHARGE :{HEM[j][2]}&@HB2        0.0585
change CHARGE :{HEM[j][2]}&@HB3        0.0585
change CHARGE :{HEM[j][2]}&@CG        -0.0477
change CHARGE :{HEM[j][2]}&@ND1       -0.0525
change CHARGE :{HEM[j][2]}&@HD1        0.3710
change CHARGE :{HEM[j][2]}&@CE1       -0.3032
change CHARGE :{HEM[j][2]}&@HE1        0.2740
change CHARGE :{HEM[j][2]}&@NE2       -0.1611
change CHARGE :{HEM[j][2]}&@CD2       -0.2274
change CHARGE :{HEM[j][2]}&@HD2        0.2297
change CHARGE :{HEM[j][0]}&@CAA        0.1105
change CHARGE :{HEM[j][0]}&@HAA1       0.0269
change CHARGE :{HEM[j][0]}&@HAA2      -0.0203
change CHARGE :{HEM[j][0]}&@CBA       -0.0733
change CHARGE :{HEM[j][0]}&@HBA1       0.0006
change CHARGE :{HEM[j][0]}&@HBA2       0.0585
change CHARGE :{HEM[j][0]}&@CGA        0.5918
change CHARGE :{HEM[j][0]}&@O1A       -0.6537
change CHARGE :{HEM[j][0]}&@O2A       -0.6853
change CHARGE :{HEM[j][0]}&@CAD       -0.0291
change CHARGE :{HEM[j][0]}&@HAD1       0.0078
change CHARGE :{HEM[j][0]}&@HAD2       0.0658
change CHARGE :{HEM[j][0]}&@CBD        0.2827
change CHARGE :{HEM[j][0]}&@HBD1      -0.0293
change CHARGE :{HEM[j][0]}&@HBD2      -0.0436
change CHARGE :{HEM[j][0]}&@CGD        0.4783
change CHARGE :{HEM[j][0]}&@O1D       -0.6025
change CHARGE :{HEM[j][0]}&@O2D       -0.6098
netcharge :{HEM[j][0]},{HEM[j][1]},{HEM[j][2]}""", file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'a'))

        print(f"""
outparm {k}{HEM[i][0]}-{l}{HEM[j][0]}.prmtop
quit
        """, file=open(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp", 'a'))

        print(f"\n Generating topology for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
        subprocess.run("parmed -i "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".inp > "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".log", shell=True)
################################################################################################################################################

################################################################################################################################################

def HemeHemeInt_HemeB(Output, Solv):

    if (os.path.isfile("pbsa.key") == True):
        print(" Found pbsa.key")
    elif (os.path.isfile("pbsa.key") == False):
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
        if (os.path.isfile("GenerateCoordForPBSA.in") == True):
            if (os.path.isfile("StrucForPBSA."+Output+".prmtop") == True) and (os.path.isfile("StrucForPBSA.rst7") == True):
                print(f" Found StrucForPBSA.{Output}.prmtop & StrucForPBSA.rst7")
        elif (os.path.isfile("GenerateCoordForPBSA.in") == False):
            print("""
 parm    %s.prmtop
 trajin  min.rst7
 strip   :WAT,Na+,Cl- outprefix StrucForPBSA
 trajout StrucForPBSA.rst7
 run
            """ %(Output), file=open("GenerateCoordForPBSA.in", 'w'))
            subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)
    elif (Solv == "implicit") or (Solv == "Implicit") or (Solv == "i") or (Solv == "I"):
        if (os.path.isfile("GenerateCoordForPBSA.in") == True):
            if (os.path.isfile("StrucForPBSA."+Output+".prmtop") == True) and (os.path.isfile("StrucForPBSA.rst7") == True):
                print(f" Found StrucForPBSA.{Output}.prmtop & StrucForPBSA.rst7")
        elif (os.path.isfile("GenerateCoordForPBSA.in") == False):
            print("""
 parm    %s.prmtop
 trajin  min.rst7
 parmwrite out StrucForPBSA.%s.prmtop
 trajout StrucForPBSA.rst7
 run
            """ %(Output, Output), file=open("GenerateCoordForPBSA.in", 'w'))
            subprocess.run("cpptraj -i GenerateCoordForPBSA.in > GenerateCoordForPBSA.log", shell=True)

    while True:
        CompChoice = input(" Do you wish to run any needed computaitons in serial or parallel (s/p)? ")

        if (CompChoice == "s") or (CompChoice == "p") or (CompChoice == "serial") or (CompChoice == "parallel"):
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    HisID = []; HebID = []; HEM = []
    with open("BondDefinitions.txt") as f:
        for line in f:
            HisID.append(int(line.strip().split('.')[1]))
            HebID.append(int(line.strip().split('.')[3]))

    for i in range(0, len(HisID), 2):
        HEM.append([HebID[i], HisID[i], HisID[i+1]])

    M = [[0 for column in range(len(HEM))] for row in range(len(HEM))]
    N = [[0 for column in range(4)] for row in range(4)]

    print("State Energies", file=open('CheckValues.txt', 'w'))
    for i in range(0, len(HEM)):
        for j in range(0, len(HEM)):
            if (j == i):
                pass
            if (j != i):
                idxc = 0
                command = ['']*(4)
                print(f"""
 ================================================
 Analyzing pair {HEM[i][0]} {HEM[j][0]}...
 ================================================""")

                for k in ("o", "r"):
                    for l in ("o", "r"):

                        if (k == "o") and (l == "o"):
                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == False):
                                PairedChargeAssignment_HemeB(Output, HEM, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "r") and (l == "o"):
                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == False):
                                PairedChargeAssignment_HemeB(Output, HEM, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "o") and (l == "r"):
                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == False):
                                PairedChargeAssignment_HemeB(Output, HEM, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1

                        if (k == "r") and (l == "r"):
                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == False):
                                PairedChargeAssignment_HemeB(Output, HEM, i, j, k, l)

                            if (os.path.isfile(str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop") == True):
                                if (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == True):
                                    print(f""" 
 Found pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 This prior output will be used for the analysis.""")
                                elif (os.path.isfile("pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])) == False):
                                    print(f""" 
 Did not find pbsa_{k}{HEM[i][0]}-{l}{HEM[j][0]} from a prior execution.
 The calculation will be submitted.""")
                                    if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                                        print(f" Running PBSA calculation for {k}{HEM[i][0]}-{l}{HEM[j][0]} ...")
                                        subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
                                    if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                                        command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+" -p "+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0])+".prmtop -c StrucForPBSA.rst7"
                                        idxc += 1
            
                if (idxc != 0):
                    commandrev = list(filter(None, command))
                    print("\n Submitting "+str(len(commandrev))+" PBSA calculations in parallel")
                    procs = [ subprocess.Popen(i, shell=True) for i in commandrev ]

                    for p in procs:
                        p.wait()
                        print("  Finished: "+str(p))

                chk = 4
                if (os.path.isfile("pbsa_o"+str(HEM[i][0])+"-o"+str(HEM[j][0])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_o{HEM[i][0]}-o{HEM[j][0]} is missing.""")
                if (os.path.isfile("pbsa_o"+str(HEM[i][0])+"-r"+str(HEM[j][0])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_o{HEM[i][0]}-r{HEM[j][0]} is missing.""")
                if (os.path.isfile("pbsa_r"+str(HEM[i][0])+"-o"+str(HEM[j][0])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_r{HEM[i][0]}-o{HEM[j][0]} is missing.""")
                if (os.path.isfile("pbsa_r"+str(HEM[i][0])+"-r"+str(HEM[j][0])) == False):
                    chk -= 1
                    print(f" Something went wrong! pbsa_r{HEM[i][0]}-r{HEM[j][0]} is missing.""")

                if (chk == 4):
                    #print(" All four files found") 
                     
                    for k in ("o", "r"):
                        for l in ("o", "r"):

                            if (k == "o") and (l == "o"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[0][0] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[0][0], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                            if (k == "o") and (l == "r"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[0][1] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[0][1], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                            if (k == "r") and (l == "o"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[1][0] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[1][0], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1
                                    
                            if (k == "r") and (l == "r"):
                                idx = 0
                                with open('pbsa_'+str(k)+str(HEM[i][0])+"-"+str(l)+str(HEM[j][0]), 'r') as fp:
                                    lines = fp.readlines()
                                    for line in lines:
                                        word1 = 'Etot'

                                        if (line.find(word1) != -1) and (idx == 0):
                                            N[1][1] = float(line.strip().split()[2]) * 0.043 
                                            print(k, i, l, j, round(N[1][1], 3), file=open('CheckValues.txt', 'a'))
                                            idx += 1

                M[i][i] = round(((N[0][1] - N[1][1]))*1000, 3)
                M[i][j] = round(((N[0][0] - N[1][0]) - (N[0][1] - N[1][1]))*1000, 3)
                M[j][i] = round(((N[0][0] - N[0][1]) - (N[1][0] - N[1][1]))*1000, 3)
                M[j][j] = round(((N[1][0] - N[1][1]))*1000, 3)
    print(*M, sep='\n', file=open('CheckValues.txt', 'a'))

################################################################################################################################################

def DeltaGFromPBSA_HemeB(Output, Solv):

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

    while True:
        CompChoice = input(" Do you wish to run any needed computaitons in serial or parallel (s/p)? ")

        if (CompChoice == "s") or (CompChoice == "p") or (CompChoice == "serial") or (CompChoice == "parallel"):
            break
        else:
            print(" Sorry, I didn't understand your choice. Please try again.")

    HisID = []; HebID = []; HEM = []
    with open("BondDefinitions.txt") as f:
        for line in f:
            HisID.append(int(line.strip().split('.')[1]))
            HebID.append(int(line.strip().split('.')[3]))

    for i in range(0, len(HisID), 2):
        HEM.append([HebID[i], HisID[i], HisID[i+1]])
    
    idxc = 0
    command = [' ']*(2*len(HEM))
    for idx in range(len(HEM)):
        print(HEM[idx][0], HEM[idx][1], HEM[idx][2])

        if (os.path.isfile("r%0d.prmtop" %(HEM[idx][0])) == False):
            print(f"""
parm StrucForPBSA.{Output}.prmtop

netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
change CHARGE :{HEM[idx][0]}&@FE         0.4223
change CHARGE :{HEM[idx][0]}&@NA         0.1973
change CHARGE :{HEM[idx][0]}&@C1A       -0.3489
change CHARGE :{HEM[idx][0]}&@C2A        0.2187
change CHARGE :{HEM[idx][0]}&@C3A       -0.0877
change CHARGE :{HEM[idx][0]}&@CMA       -0.0598
change CHARGE :{HEM[idx][0]}&@HMA1       0.0703
change CHARGE :{HEM[idx][0]}&@HMA2       0.0281
change CHARGE :{HEM[idx][0]}&@HMA3       0.0872
change CHARGE :{HEM[idx][0]}&@C4A       -0.2107
change CHARGE :{HEM[idx][0]}&@CHB       -0.1614
change CHARGE :{HEM[idx][0]}&@HHB        0.3441
change CHARGE :{HEM[idx][0]}&@C1B       -0.1858
change CHARGE :{HEM[idx][0]}&@NB         0.1413
change CHARGE :{HEM[idx][0]}&@C2B       -0.0884   
change CHARGE :{HEM[idx][0]}&@CMB        0.0466
change CHARGE :{HEM[idx][0]}&@HMB1      -0.0146
change CHARGE :{HEM[idx][0]}&@HMB2       0.0447
change CHARGE :{HEM[idx][0]}&@HMB3       0.0352
change CHARGE :{HEM[idx][0]}&@C3B        0.0348
change CHARGE :{HEM[idx][0]}&@CAB       -0.2072
change CHARGE :{HEM[idx][0]}&@HAB        0.1898
change CHARGE :{HEM[idx][0]}&@CBB       -0.4301
change CHARGE :{HEM[idx][0]}&@HBB1       0.2165
change CHARGE :{HEM[idx][0]}&@HBB2       0.1971
change CHARGE :{HEM[idx][0]}&@C4B       -0.0142
change CHARGE :{HEM[idx][0]}&@CHC       -0.2576
change CHARGE :{HEM[idx][0]}&@HHC        0.2520
change CHARGE :{HEM[idx][0]}&@C1C        0.0528
change CHARGE :{HEM[idx][0]}&@NC         0.1307
change CHARGE :{HEM[idx][0]}&@C2C        0.0116
change CHARGE :{HEM[idx][0]}&@CMC       -0.1501
change CHARGE :{HEM[idx][0]}&@HMC1       0.0991
change CHARGE :{HEM[idx][0]}&@HMC2       0.1038
change CHARGE :{HEM[idx][0]}&@HMC3       0.0489
change CHARGE :{HEM[idx][0]}&@C3C        0.0046
change CHARGE :{HEM[idx][0]}&@CAC       -0.0156
change CHARGE :{HEM[idx][0]}&@HAC        0.1293
change CHARGE :{HEM[idx][0]}&@CBC       -0.5384
change CHARGE :{HEM[idx][0]}&@HBC1       0.2611
change CHARGE :{HEM[idx][0]}&@HBC2       0.2581
change CHARGE :{HEM[idx][0]}&@C4C       -0.2427
change CHARGE :{HEM[idx][0]}&@CHD       -0.2138
change CHARGE :{HEM[idx][0]}&@HHD        0.3122
change CHARGE :{HEM[idx][0]}&@C1D       -0.0975
change CHARGE :{HEM[idx][0]}&@ND         0.1255
change CHARGE :{HEM[idx][0]}&@C2D       -0.0861
change CHARGE :{HEM[idx][0]}&@CMD        0.0005
change CHARGE :{HEM[idx][0]}&@HMD1       0.0098
change CHARGE :{HEM[idx][0]}&@HMD2       0.0805
change CHARGE :{HEM[idx][0]}&@HMD3       0.0078
change CHARGE :{HEM[idx][0]}&@C3D       -0.1965
change CHARGE :{HEM[idx][0]}&@C4D        0.0236
change CHARGE :{HEM[idx][0]}&@CHA       -0.1216
change CHARGE :{HEM[idx][0]}&@HHA        0.3063
change CHARGE :{HEM[idx][1]}&@CB         0.0426
change CHARGE :{HEM[idx][1]}&@HB2        0.0349
change CHARGE :{HEM[idx][1]}&@HB3        0.0349
change CHARGE :{HEM[idx][1]}&@CG         0.0408
change CHARGE :{HEM[idx][1]}&@ND1       -0.1881
change CHARGE :{HEM[idx][1]}&@HD1        0.3920
change CHARGE :{HEM[idx][1]}&@CE1       -0.1901
change CHARGE :{HEM[idx][1]}&@HE1        0.2191
change CHARGE :{HEM[idx][1]}&@NE2       -0.0847
change CHARGE :{HEM[idx][1]}&@CD2       -0.3012
change CHARGE :{HEM[idx][1]}&@HD2        0.2011
change CHARGE :{HEM[idx][2]}&@CB        -0.0373
change CHARGE :{HEM[idx][2]}&@HB2        0.0585
change CHARGE :{HEM[idx][2]}&@HB3        0.0585
change CHARGE :{HEM[idx][2]}&@CG        -0.0477
change CHARGE :{HEM[idx][2]}&@ND1       -0.0525
change CHARGE :{HEM[idx][2]}&@HD1        0.3710
change CHARGE :{HEM[idx][2]}&@CE1       -0.3032
change CHARGE :{HEM[idx][2]}&@HE1        0.2740
change CHARGE :{HEM[idx][2]}&@NE2       -0.1611
change CHARGE :{HEM[idx][2]}&@CD2       -0.2274
change CHARGE :{HEM[idx][2]}&@HD2        0.2297
change CHARGE :{HEM[idx][0]}&@CAA        0.1105
change CHARGE :{HEM[idx][0]}&@HAA1       0.0269
change CHARGE :{HEM[idx][0]}&@HAA2      -0.0203
change CHARGE :{HEM[idx][0]}&@CBA       -0.0733
change CHARGE :{HEM[idx][0]}&@HBA1       0.0006
change CHARGE :{HEM[idx][0]}&@HBA2       0.0585
change CHARGE :{HEM[idx][0]}&@CGA        0.5918
change CHARGE :{HEM[idx][0]}&@O1A       -0.6537
change CHARGE :{HEM[idx][0]}&@O2A       -0.6853
change CHARGE :{HEM[idx][0]}&@CAD       -0.0291
change CHARGE :{HEM[idx][0]}&@HAD1       0.0078
change CHARGE :{HEM[idx][0]}&@HAD2       0.0658
change CHARGE :{HEM[idx][0]}&@CBD        0.2827
change CHARGE :{HEM[idx][0]}&@HBD1      -0.0293
change CHARGE :{HEM[idx][0]}&@HBD2      -0.0436
change CHARGE :{HEM[idx][0]}&@CGD        0.4783
change CHARGE :{HEM[idx][0]}&@O1D       -0.6025
change CHARGE :{HEM[idx][0]}&@O2D       -0.6098
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
outparm r{HEM[idx][0]}.prmtop
""", file=open("r"+str(HEM[idx][0])+".inp", 'w'))

            print("\n Generating topology for rHEH-%0d ..." %(HEM[idx][0]))
            subprocess.run("parmed -i "+str("r")+str(HEM[idx][0])+".inp > "+str("r")+str(HEM[idx][0])+".log", shell=True)

        if (os.path.isfile("pbsa_r"+str(HEM[idx][0])) == True):
            print(""" 
 Found pbsa_r%0d from a prior execution.
 This prior output will be used for the analysis.""" %(HEM[idx][0]))
        else: 
            if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                print(" Running PBSA calculation for red. HEH-%0d ..." %(HEM[idx][0]))
                subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str("r")+str(HEM[idx][0])+" -p "+str("r")+str(HEM[idx][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
            if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str("r")+str(HEM[idx][0])+" -p "+str("r")+str(HEM[idx][0])+".prmtop -c StrucForPBSA.rst7"
                idxc += 1
            
        if (os.path.isfile("o%0d.prmtop" %(HEM[idx][0])) == False):
            print(f"""
parm StrucForPBSA.{Output}.prmtop

netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
change CHARGE :{HEM[idx][0]}&@FE         0.6083
change CHARGE :{HEM[idx][0]}&@NA         0.1780
change CHARGE :{HEM[idx][0]}&@C1A       -0.2990
change CHARGE :{HEM[idx][0]}&@C2A        0.2176
change CHARGE :{HEM[idx][0]}&@C3A       -0.0578
change CHARGE :{HEM[idx][0]}&@CMA       -0.0653
change CHARGE :{HEM[idx][0]}&@HMA1       0.0838
change CHARGE :{HEM[idx][0]}&@HMA2       0.0416
change CHARGE :{HEM[idx][0]}&@HMA3       0.1007
change CHARGE :{HEM[idx][0]}&@C4A       -0.2008
change CHARGE :{HEM[idx][0]}&@CHB       -0.1385
change CHARGE :{HEM[idx][0]}&@HHB        0.3521
change CHARGE :{HEM[idx][0]}&@C1B       -0.1639
change CHARGE :{HEM[idx][0]}&@NB         0.1343
change CHARGE :{HEM[idx][0]}&@C2B       -0.1095
change CHARGE :{HEM[idx][0]}&@CMB        0.0161
change CHARGE :{HEM[idx][0]}&@HMB1       0.0095
change CHARGE :{HEM[idx][0]}&@HMB2       0.0688
change CHARGE :{HEM[idx][0]}&@HMB3       0.0593
change CHARGE :{HEM[idx][0]}&@C3B        0.1107
change CHARGE :{HEM[idx][0]}&@CAB       -0.2074
change CHARGE :{HEM[idx][0]}&@HAB        0.1943
change CHARGE :{HEM[idx][0]}&@CBB       -0.4054
change CHARGE :{HEM[idx][0]}&@HBB1       0.2192
change CHARGE :{HEM[idx][0]}&@HBB2       0.1998
change CHARGE :{HEM[idx][0]}&@C4B       -0.0043
change CHARGE :{HEM[idx][0]}&@CHC       -0.2657
change CHARGE :{HEM[idx][0]}&@HHC        0.2630
change CHARGE :{HEM[idx][0]}&@C1C        0.1417
change CHARGE :{HEM[idx][0]}&@NC         0.0897
change CHARGE :{HEM[idx][0]}&@C2C       -0.0225
change CHARGE :{HEM[idx][0]}&@CMC       -0.1841
change CHARGE :{HEM[idx][0]}&@HMC1       0.1221
change CHARGE :{HEM[idx][0]}&@HMC2       0.1268
change CHARGE :{HEM[idx][0]}&@HMC3       0.0719
change CHARGE :{HEM[idx][0]}&@C3C        0.0485
change CHARGE :{HEM[idx][0]}&@CAC        0.0012
change CHARGE :{HEM[idx][0]}&@HAC        0.1508
change CHARGE :{HEM[idx][0]}&@CBC       -0.5097
change CHARGE :{HEM[idx][0]}&@HBC1       0.2648
change CHARGE :{HEM[idx][0]}&@HBC2       0.2618
change CHARGE :{HEM[idx][0]}&@C4C       -0.2418
change CHARGE :{HEM[idx][0]}&@CHD       -0.1909
change CHARGE :{HEM[idx][0]}&@HHD        0.3192
change CHARGE :{HEM[idx][0]}&@C1D       -0.0856
change CHARGE :{HEM[idx][0]}&@ND         0.1145
change CHARGE :{HEM[idx][0]}&@C2D       -0.0592
change CHARGE :{HEM[idx][0]}&@CMD       -0.0045
change CHARGE :{HEM[idx][0]}&@HMD1       0.0235
change CHARGE :{HEM[idx][0]}&@HMD2       0.0942
change CHARGE :{HEM[idx][0]}&@HMD3       0.0215
change CHARGE :{HEM[idx][0]}&@C3D       -0.1936
change CHARGE :{HEM[idx][0]}&@C4D        0.0685
change CHARGE :{HEM[idx][0]}&@CHA       -0.1387
change CHARGE :{HEM[idx][0]}&@HHA        0.3203
change CHARGE :{HEM[idx][1]}&@CB         0.0590
change CHARGE :{HEM[idx][1]}&@HB2        0.0431
change CHARGE :{HEM[idx][1]}&@HB3        0.0431
change CHARGE :{HEM[idx][1]}&@CG         0.0368
change CHARGE :{HEM[idx][1]}&@ND1       -0.1361
change CHARGE :{HEM[idx][1]}&@HD1        0.4020
change CHARGE :{HEM[idx][1]}&@CE1       -0.1891
change CHARGE :{HEM[idx][1]}&@HE1        0.2281
change CHARGE :{HEM[idx][1]}&@NE2       -0.0977
change CHARGE :{HEM[idx][1]}&@CD2       -0.2862
change CHARGE :{HEM[idx][1]}&@HD2        0.2051
change CHARGE :{HEM[idx][2]}&@CB        -0.0209
change CHARGE :{HEM[idx][2]}&@HB2        0.0667
change CHARGE :{HEM[idx][2]}&@HB3        0.0667
change CHARGE :{HEM[idx][2]}&@CG        -0.0517
change CHARGE :{HEM[idx][2]}&@ND1       -0.0005
change CHARGE :{HEM[idx][2]}&@HD1        0.3810
change CHARGE :{HEM[idx][2]}&@CE1       -0.3022
change CHARGE :{HEM[idx][2]}&@HE1        0.2830
change CHARGE :{HEM[idx][2]}&@NE2       -0.1741
change CHARGE :{HEM[idx][2]}&@CD2       -0.2124
change CHARGE :{HEM[idx][2]}&@HD2        0.2337
change CHARGE :{HEM[idx][0]}&@CAA        0.1105
change CHARGE :{HEM[idx][0]}&@HAA1       0.0269
change CHARGE :{HEM[idx][0]}&@HAA2      -0.0203
change CHARGE :{HEM[idx][0]}&@CBA       -0.0733
change CHARGE :{HEM[idx][0]}&@HBA1       0.0006
change CHARGE :{HEM[idx][0]}&@HBA2       0.0585
change CHARGE :{HEM[idx][0]}&@CGA        0.5918
change CHARGE :{HEM[idx][0]}&@O1A       -0.6537
change CHARGE :{HEM[idx][0]}&@O2A       -0.6853
change CHARGE :{HEM[idx][0]}&@CAD       -0.0291
change CHARGE :{HEM[idx][0]}&@HAD1       0.0078
change CHARGE :{HEM[idx][0]}&@HAD2       0.0658
change CHARGE :{HEM[idx][0]}&@CBD        0.2827
change CHARGE :{HEM[idx][0]}&@HBD1      -0.0293
change CHARGE :{HEM[idx][0]}&@HBD2      -0.0436
change CHARGE :{HEM[idx][0]}&@CGD        0.4783
change CHARGE :{HEM[idx][0]}&@O1D       -0.6025
change CHARGE :{HEM[idx][0]}&@O2D       -0.6098
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
outparm o{HEM[idx][0]}.prmtop
""", file=open("o"+str(HEM[idx][0])+".inp", 'w'))

            print("\n Generating topology for oHEH-%0d ..." %(HEM[idx][0]))
            subprocess.run("parmed -i "+str("o")+str(HEM[idx][0])+".inp > "+str("o")+str(HEM[idx][0])+".log", shell=True)

        if (os.path.isfile("pbsa_o"+str(HEM[idx][0])) == True):
            print(""" 
 Found pbsa_o%0d from a prior execution.
 This prior output will be used for the analysis.""" %(HEM[idx][0]))
        else: 
            if (CompChoice == "serial") or (CompChoice == "s") or (CompChoice == "S"):
                print(" Running PBSA calculation for ox. HEH-%0d ..." %(HEM[idx][0]))
                subprocess.run("pbsa -O -i pbsa.key -o pbsa_"+str("o")+str(HEM[idx][0])+" -p "+str("o")+str(HEM[idx][0])+".prmtop -c StrucForPBSA.rst7", shell=True)
            if (CompChoice == "parallel") or (CompChoice == "p") or (CompChoice == "P"):
                command[idxc] = "pbsa -O -i pbsa.key -o pbsa_"+str("o")+str(HEM[idx][0])+" -p "+str("o")+str(HEM[idx][0])+".prmtop -c StrucForPBSA.rst7"
                idxc += 1

    if (idxc != 0):
        print("\n Submitting "+str(len(command))+" PBSA calculations in parallel")

        procs = [ subprocess.Popen(i, shell=True) for i in command ]
        for p in procs:
            p.wait()
            print("  Finished: "+str(p))

    DEtot = [0]*(len(HEM))
    DEelec = [0]*(len(HEM))
    DG = [0]*(len(HEM)-1)
    for idx in range(len(HEM)):

        idx1 = 0; idx2 = 0;
        with open('pbsa_o'+str(HEM[idx][0]), 'r') as fp:
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
        with open('pbsa_r'+str(HEM[idx][0]), 'r') as fp:
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

        print("""  step=%0d HEM-%0d EtotOx=%.3f eV EtotRed=%.3f eV DEtot=%.3f eV EelecOx=%.3f eV EelecRed=%.3f eV DEelec=%.3f eV""" %(idx, HEM[idx][0], EtotOx, EtotRed, DEtot[idx], EelecOx, EelecRed, DEelec[idx]))

    for idx in range(len(HEM)-1):
        DG[idx] = -1 * ((-1 * DEtot[idx]) + (DEtot[idx+1])) 

        if (idx == 0):
            print("(HEM-%0d = %.3f eV) -> (HEM-%0d = %.3f eV); DG = %10.3f eV" %(HEM[idx][0],  DEtot[idx], HEM[idx+1][0], DEtot[idx+1], DG[idx]), file=open('DG.txt', 'w'))
        else:
            print("(HEM-%0d = %.3f eV) -> (HEM-%0d = %.3f eV); DG = %10.3f eV" %(HEM[idx][0],  DEtot[idx], HEM[idx+1][0], DEtot[idx+1], DG[idx]), file=open('DG.txt', 'a'))

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

################################################################################################################################################

def DeltaGFromLIE_HemeB(Output):

    idx = 0
    HisID = []; HebID = []; HEM = []
    with open("BondDefinitions.txt") as f:
        for line in f:
            HisID.append(int(line.strip().split('.')[1]))
            HebID.append(int(line.strip().split('.')[3]))

    for i in range(0, len(HisID), 2):
        HEM.append([HebID[i], HisID[i], HisID[i+1]])

    for idx in range(len(HEM)):
       #print(HEM[idx][0], HEM[idx][1], HEM[idx][2])
        HEB = str(HEB[idx][0])+","+str(HEM[idx][1])+","+str(HEM[idx][2])

        if (os.path.isfile("o%0d.prmtop" %(HEM[idx][0])) == False):
            print("""
parm %s.prmtop""".format(Output), file=open("o"+str(HEM[idx][0])+".inp", 'w'))

            print(f"""
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
change CHARGE :{HEM[idx][0]}&@FE         0.6083
change CHARGE :{HEM[idx][0]}&@NA         0.1780
change CHARGE :{HEM[idx][0]}&@C1A       -0.2990
change CHARGE :{HEM[idx][0]}&@C2A        0.2176
change CHARGE :{HEM[idx][0]}&@C3A       -0.0578
change CHARGE :{HEM[idx][0]}&@CMA       -0.0653
change CHARGE :{HEM[idx][0]}&@HMA1       0.0838
change CHARGE :{HEM[idx][0]}&@HMA2       0.0416
change CHARGE :{HEM[idx][0]}&@HMA3       0.1007
change CHARGE :{HEM[idx][0]}&@C4A       -0.2008
change CHARGE :{HEM[idx][0]}&@CHB       -0.1385
change CHARGE :{HEM[idx][0]}&@HHB        0.3521
change CHARGE :{HEM[idx][0]}&@C1B       -0.1639
change CHARGE :{HEM[idx][0]}&@NB         0.1343
change CHARGE :{HEM[idx][0]}&@C2B       -0.1095
change CHARGE :{HEM[idx][0]}&@CMB        0.0161
change CHARGE :{HEM[idx][0]}&@HMB1       0.0095
change CHARGE :{HEM[idx][0]}&@HMB2       0.0688
change CHARGE :{HEM[idx][0]}&@HMB3       0.0593
change CHARGE :{HEM[idx][0]}&@C3B        0.1107
change CHARGE :{HEM[idx][0]}&@CAB       -0.2074
change CHARGE :{HEM[idx][0]}&@HAB        0.1943
change CHARGE :{HEM[idx][0]}&@CBB       -0.4054
change CHARGE :{HEM[idx][0]}&@HBB1       0.2192
change CHARGE :{HEM[idx][0]}&@HBB2       0.1998
change CHARGE :{HEM[idx][0]}&@C4B       -0.0043
change CHARGE :{HEM[idx][0]}&@CHC       -0.2657
change CHARGE :{HEM[idx][0]}&@HHC        0.2630
change CHARGE :{HEM[idx][0]}&@C1C        0.1417
change CHARGE :{HEM[idx][0]}&@NC         0.0897
change CHARGE :{HEM[idx][0]}&@C2C       -0.0225
change CHARGE :{HEM[idx][0]}&@CMC       -0.1841
change CHARGE :{HEM[idx][0]}&@HMC1       0.1221
change CHARGE :{HEM[idx][0]}&@HMC2       0.1268
change CHARGE :{HEM[idx][0]}&@HMC3       0.0719
change CHARGE :{HEM[idx][0]}&@C3C        0.0485
change CHARGE :{HEM[idx][0]}&@CAC        0.0012
change CHARGE :{HEM[idx][0]}&@HAC        0.1508
change CHARGE :{HEM[idx][0]}&@CBC       -0.5097
change CHARGE :{HEM[idx][0]}&@HBC1       0.2648
change CHARGE :{HEM[idx][0]}&@HBC2       0.2618
change CHARGE :{HEM[idx][0]}&@C4C       -0.2418
change CHARGE :{HEM[idx][0]}&@CHD       -0.1909
change CHARGE :{HEM[idx][0]}&@HHD        0.3192
change CHARGE :{HEM[idx][0]}&@C1D       -0.0856
change CHARGE :{HEM[idx][0]}&@ND         0.1145
change CHARGE :{HEM[idx][0]}&@C2D       -0.0592
change CHARGE :{HEM[idx][0]}&@CMD       -0.0045
change CHARGE :{HEM[idx][0]}&@HMD1       0.0235
change CHARGE :{HEM[idx][0]}&@HMD2       0.0942
change CHARGE :{HEM[idx][0]}&@HMD3       0.0215
change CHARGE :{HEM[idx][0]}&@C3D       -0.1936
change CHARGE :{HEM[idx][0]}&@C4D        0.0685
change CHARGE :{HEM[idx][0]}&@CHA       -0.1387
change CHARGE :{HEM[idx][0]}&@HHA        0.3203
change CHARGE :{HEM[idx][1]}&@CB         0.0590
change CHARGE :{HEM[idx][1]}&@HB2        0.0431
change CHARGE :{HEM[idx][1]}&@HB3        0.0431
change CHARGE :{HEM[idx][1]}&@CG         0.0368
change CHARGE :{HEM[idx][1]}&@ND1       -0.1361
change CHARGE :{HEM[idx][1]}&@HD1        0.4020
change CHARGE :{HEM[idx][1]}&@CE1       -0.1891
change CHARGE :{HEM[idx][1]}&@HE1        0.2281
change CHARGE :{HEM[idx][1]}&@NE2       -0.0977
change CHARGE :{HEM[idx][1]}&@CD2       -0.2862
change CHARGE :{HEM[idx][1]}&@HD2        0.2051
change CHARGE :{HEM[idx][2]}&@CB        -0.0209
change CHARGE :{HEM[idx][2]}&@HB2        0.0667
change CHARGE :{HEM[idx][2]}&@HB3        0.0667
change CHARGE :{HEM[idx][2]}&@CG        -0.0517
change CHARGE :{HEM[idx][2]}&@ND1       -0.0005
change CHARGE :{HEM[idx][2]}&@HD1        0.3810
change CHARGE :{HEM[idx][2]}&@CE1       -0.3022
change CHARGE :{HEM[idx][2]}&@HE1        0.2830
change CHARGE :{HEM[idx][2]}&@NE2       -0.1741
change CHARGE :{HEM[idx][2]}&@CD2       -0.2124
change CHARGE :{HEM[idx][2]}&@HD2        0.2337
change CHARGE :{HEM[idx][0]}&@CAA        0.1105
change CHARGE :{HEM[idx][0]}&@HAA1       0.0269
change CHARGE :{HEM[idx][0]}&@HAA2      -0.0203
change CHARGE :{HEM[idx][0]}&@CBA       -0.0733
change CHARGE :{HEM[idx][0]}&@HBA1       0.0006
change CHARGE :{HEM[idx][0]}&@HBA2       0.0585
change CHARGE :{HEM[idx][0]}&@CGA        0.5918
change CHARGE :{HEM[idx][0]}&@O1A       -0.6537
change CHARGE :{HEM[idx][0]}&@O2A       -0.6853
change CHARGE :{HEM[idx][0]}&@CAD       -0.0291
change CHARGE :{HEM[idx][0]}&@HAD1       0.0078
change CHARGE :{HEM[idx][0]}&@HAD2       0.0658
change CHARGE :{HEM[idx][0]}&@CBD        0.2827
change CHARGE :{HEM[idx][0]}&@HBD1      -0.0293
change CHARGE :{HEM[idx][0]}&@HBD2      -0.0436
change CHARGE :{HEM[idx][0]}&@CGD        0.4783
change CHARGE :{HEM[idx][0]}&@O1D       -0.6025
change CHARGE :{HEM[idx][0]}&@O2D       -0.6098
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
outparm o{HEM[idx][0]}.prmtop
""", file=open("o"+str(HEM[idx][0])+".inp", 'a'))

            print("\n Generating oxidized topology for HEB-%0d ..." %(HEM[idx][0]))
            subprocess.run("parmed -i o"+str(HEM[idx][0])+".inp > r"+str(HEM[idx][0])+".log", shell=True)

        if (os.path.isfile("r%0d.prmtop" %(HEM[idx][0])) == False):
            print("""
parm %s.prmtop""".format(Output), file=open("r"+str(HEM[idx][0])+".inp", 'w'))

            print(f"""
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
change CHARGE :{HEM[idx][0]}&@FE         0.4223
change CHARGE :{HEM[idx][0]}&@NA         0.1973
change CHARGE :{HEM[idx][0]}&@C1A       -0.3489
change CHARGE :{HEM[idx][0]}&@C2A        0.2187
change CHARGE :{HEM[idx][0]}&@C3A       -0.0877
change CHARGE :{HEM[idx][0]}&@CMA       -0.0598
change CHARGE :{HEM[idx][0]}&@HMA1       0.0703
change CHARGE :{HEM[idx][0]}&@HMA2       0.0281
change CHARGE :{HEM[idx][0]}&@HMA3       0.0872
change CHARGE :{HEM[idx][0]}&@C4A       -0.2107
change CHARGE :{HEM[idx][0]}&@CHB       -0.1614
change CHARGE :{HEM[idx][0]}&@HHB        0.3441
change CHARGE :{HEM[idx][0]}&@C1B       -0.1858
change CHARGE :{HEM[idx][0]}&@NB         0.1413
change CHARGE :{HEM[idx][0]}&@C2B       -0.0884   
change CHARGE :{HEM[idx][0]}&@CMB        0.0466
change CHARGE :{HEM[idx][0]}&@HMB1      -0.0146
change CHARGE :{HEM[idx][0]}&@HMB2       0.0447
change CHARGE :{HEM[idx][0]}&@HMB3       0.0352
change CHARGE :{HEM[idx][0]}&@C3B        0.0348
change CHARGE :{HEM[idx][0]}&@CAB       -0.2072
change CHARGE :{HEM[idx][0]}&@HAB        0.1898
change CHARGE :{HEM[idx][0]}&@CBB       -0.4301
change CHARGE :{HEM[idx][0]}&@HBB1       0.2165
change CHARGE :{HEM[idx][0]}&@HBB2       0.1971
change CHARGE :{HEM[idx][0]}&@C4B       -0.0142
change CHARGE :{HEM[idx][0]}&@CHC       -0.2576
change CHARGE :{HEM[idx][0]}&@HHC        0.2520
change CHARGE :{HEM[idx][0]}&@C1C        0.0528
change CHARGE :{HEM[idx][0]}&@NC         0.1307
change CHARGE :{HEM[idx][0]}&@C2C        0.0116
change CHARGE :{HEM[idx][0]}&@CMC       -0.1501
change CHARGE :{HEM[idx][0]}&@HMC1       0.0991
change CHARGE :{HEM[idx][0]}&@HMC2       0.1038
change CHARGE :{HEM[idx][0]}&@HMC3       0.0489
change CHARGE :{HEM[idx][0]}&@C3C        0.0046
change CHARGE :{HEM[idx][0]}&@CAC       -0.0156
change CHARGE :{HEM[idx][0]}&@HAC        0.1293
change CHARGE :{HEM[idx][0]}&@CBC       -0.5384
change CHARGE :{HEM[idx][0]}&@HBC1       0.2611
change CHARGE :{HEM[idx][0]}&@HBC2       0.2581
change CHARGE :{HEM[idx][0]}&@C4C       -0.2427
change CHARGE :{HEM[idx][0]}&@CHD       -0.2138
change CHARGE :{HEM[idx][0]}&@HHD        0.3122
change CHARGE :{HEM[idx][0]}&@C1D       -0.0975
change CHARGE :{HEM[idx][0]}&@ND         0.1255
change CHARGE :{HEM[idx][0]}&@C2D       -0.0861
change CHARGE :{HEM[idx][0]}&@CMD        0.0005
change CHARGE :{HEM[idx][0]}&@HMD1       0.0098
change CHARGE :{HEM[idx][0]}&@HMD2       0.0805
change CHARGE :{HEM[idx][0]}&@HMD3       0.0078
change CHARGE :{HEM[idx][0]}&@C3D       -0.1965
change CHARGE :{HEM[idx][0]}&@C4D        0.0236
change CHARGE :{HEM[idx][0]}&@CHA       -0.1216
change CHARGE :{HEM[idx][0]}&@HHA        0.3063
change CHARGE :{HEM[idx][1]}&@CB         0.0426
change CHARGE :{HEM[idx][1]}&@HB2        0.0349
change CHARGE :{HEM[idx][1]}&@HB3        0.0349
change CHARGE :{HEM[idx][1]}&@CG         0.0408
change CHARGE :{HEM[idx][1]}&@ND1       -0.1881
change CHARGE :{HEM[idx][1]}&@HD1        0.3920
change CHARGE :{HEM[idx][1]}&@CE1       -0.1901
change CHARGE :{HEM[idx][1]}&@HE1        0.2191
change CHARGE :{HEM[idx][1]}&@NE2       -0.0847
change CHARGE :{HEM[idx][1]}&@CD2       -0.3012
change CHARGE :{HEM[idx][1]}&@HD2        0.2011
change CHARGE :{HEM[idx][2]}&@CB        -0.0373
change CHARGE :{HEM[idx][2]}&@HB2        0.0585
change CHARGE :{HEM[idx][2]}&@HB3        0.0585
change CHARGE :{HEM[idx][2]}&@CG        -0.0477
change CHARGE :{HEM[idx][2]}&@ND1       -0.0525
change CHARGE :{HEM[idx][2]}&@HD1        0.3710
change CHARGE :{HEM[idx][2]}&@CE1       -0.3032
change CHARGE :{HEM[idx][2]}&@HE1        0.2740
change CHARGE :{HEM[idx][2]}&@NE2       -0.1611
change CHARGE :{HEM[idx][2]}&@CD2       -0.2274
change CHARGE :{HEM[idx][2]}&@HD2        0.2297
change CHARGE :{HEM[idx][0]}&@CAA        0.1105
change CHARGE :{HEM[idx][0]}&@HAA1       0.0269
change CHARGE :{HEM[idx][0]}&@HAA2      -0.0203
change CHARGE :{HEM[idx][0]}&@CBA       -0.0733
change CHARGE :{HEM[idx][0]}&@HBA1       0.0006
change CHARGE :{HEM[idx][0]}&@HBA2       0.0585
change CHARGE :{HEM[idx][0]}&@CGA        0.5918
change CHARGE :{HEM[idx][0]}&@O1A       -0.6537
change CHARGE :{HEM[idx][0]}&@O2A       -0.6853
change CHARGE :{HEM[idx][0]}&@CAD       -0.0291
change CHARGE :{HEM[idx][0]}&@HAD1       0.0078
change CHARGE :{HEM[idx][0]}&@HAD2       0.0658
change CHARGE :{HEM[idx][0]}&@CBD        0.2827
change CHARGE :{HEM[idx][0]}&@HBD1      -0.0293
change CHARGE :{HEM[idx][0]}&@HBD2      -0.0436
change CHARGE :{HEM[idx][0]}&@CGD        0.4783
change CHARGE :{HEM[idx][0]}&@O1D       -0.6025
change CHARGE :{HEM[idx][0]}&@O2D       -0.6098
netcharge :{HEM[idx][0]},{HEM[idx][1]},{HEM[idx][2]}
outparm r{HEM[idx][0]}.prmtop
""", file=open("r"+str(HEM[idx][0])+".inp", 'a'))

            print(" Generating reduced topology for HEH-%0d ..." %(HEM[idx][0]))
            subprocess.run("parmed -i r"+str(HEM[idx][0])+".inp > r"+str(HEM[idx][0])+".log", shell=True)

            if (os.path.isfile("%0dphsicochemAvg.dat" %(HEM[idx][0])) == True):
                print(""" 
 Found %0dphsicochemAvg.dat from a prior execution.
 This prior output will be used for the analysis.""" %(HEH[idx])) 

            if (os.path.isfile("%0dphsicochemAvg.dat" %(HEM[idx][0])) == False):
                if (os.path.isfile("o%0d.prmtop" %(HEM[idx][0])) == True) and (os.path.isfile("r%0d.prmtop" %(HEM[idx][0])) == True): 
                    print("""
 parm   o%0d.prmtop""".format(HEM[idx][0]), file=open("%0dredox.in" %(HEM[idx][0]), 'w'))

                    print(f"""
 trajin min.pdb parmindex 0

 lie OHSys      :({HEB})        !:({HEB})                        out  OHSys.dat
 lie OHP        :({HEB})        !(:{HEB},WAT,Na+,Cl-)            out  OHP.dat
 lie OHS        :({HEB})        :WAT,Na+,Cl-                     out  OHS.dat
 lie OHNonPolar :({HEB})        :GLY,ALA,VAL,LEU,ILE,MET,PRO     out  OHNonPolar.dat
 lie OHAromatic :({HEB})        :PHE,TYR,TRP                     out  OHAromatic.dat
 lie OHPolar    :({HEB})        :SER,THR,CYS,CYO,ASN,GLN         out  OHPolar.dat
 lie OHAcidic   :({HEB})        :ASP,AS4,GLU,GL4                 out  OHAcidic.dat 
 lie OHBasic    :({HEB})        :HIS,HIP,HID,HIE,HIO,LYS,ARG     out  OHBasic.dat 
 lie OHHEB      :({HEB})        :HEM,HIM,HIN&!:{HEB}                  out  OHHEB.dat 
 run

 clear trajin
                    """ %{'HEB': HEB}, file=open("%0dredox.in" %(HEM[idx][0]), 'a'))

                    print("""
 parm   r%0d.prmtop""".format(HEM[idx][0]), file=open("%0dredox.in" %(HEM[idx][0]), 'a'))

                    print(f"""
 trajin min.pdb parmindex 1

 lie RHSys      :({HEB})        !:{HEB}                          out  RHSys.dat
 lie RHP        :({HEB})        !(:{HEB},WAT,Na+,Cl-)            out  RHP.dat
 lie RHS        :({HEB})        :WAT,Na+,Cl-                     out  RHS.dat
 lie RHNonPolar :({HEB})        :GLY,ALA,VAL,LEU,ILE,MET,PRO     out  RHNonPolar.dat
 lie RHAromatic :({HEB})        :PHE,TYR,TRP                     out  RHAromatic.dat
 lie RHPolar    :({HEB})        :SER,THR,CYS,CYO,ASN,GLN         out  RHPolar.dat
 lie RHAcidic   :({HEB})        :ASP,AS4,GLU,GL4                 out  RHAcidic.dat 
 lie RHBasic    :({HEB})        :HIS,HIP,HID,HIE,HIO,LYS,ARG     out  RHBasic.dat 
 lie RHHEB      :({HEB})        :HEB&!:{HEB}                     out  RHHEB.dat 
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
 ODHHEB      = avg(OHHEB[EELEC])      * 0.043

 RDHSys      = avg(RHSys[EELEC])      * 0.043 
 RDHP        = avg(RHP[EELEC])        * 0.043 
 RDHS        = avg(RHS[EELEC])        * 0.043
 RDHNonPolar = avg(RHNonPolar[EELEC]) * 0.043
 RDHAromatic = avg(RHAromatic[EELEC]) * 0.043
 RDHPolar    = avg(RHPolar[EELEC])    * 0.043
 RDHAcidic   = avg(RHAcidic[EELEC])   * 0.043
 RDHBasic    = avg(RHBasic[EELEC])    * 0.043
 RDHHEB      = avg(RHHEB[EELEC])      * 0.043

 DHSys          = (avg(OHSys[EELEC])      - avg(RHSys[EELEC]))      * 0.043    
 DHP            = (avg(OHP[EELEC])        - avg(RHP[EELEC]))        * 0.043 
 DHS            = (avg(OHS[EELEC])        - avg(RHS[EELEC]))        * 0.043 
 DHNonPolar     = (avg(OHNonPolar[EELEC]) - avg(RHNonPolar[EELEC])) * 0.043    
 DHAromatic     = (avg(OHAromatic[EELEC]) - avg(RHAromatic[EELEC])) * 0.043 
 DHPolar        = (avg(OHPolar[EELEC])    - avg(RHPolar[EELEC]))    * 0.043    
 DHAcidic       = (avg(OHAcidic[EELEC])   - avg(RHAcidic[EELEC]))   * 0.043 
 DHBasic        = (avg(OHBasic[EELEC])    - avg(RHBasic[EELEC]))    * 0.043
 DHHEB          = (avg(OHHEB[EELEC])      - avg(RHHEB[EELEC]))      * 0.043

 #DHSysstd       = stdev(DHSys)       * 0.043    
 #DHPstd         = stdev(DHP)         * 0.043 
 #DHSstd         = stdev(DHS)         * 0.043 
 #DHNonPolarstd  = stdev(DHNonPolar)  * 0.043
 #DHAromaticstd  = stdev(DHAromatic)  * 0.043 
 #DHPolarstd     = stdev(DHPolar)     * 0.043 
 #DHAcidicstd    = stdev(DHAcidic)    * 0.043 
 #DHBasicstd     = stdev(DHBasic)     * 0.043 
 #DHHEBstd       = stdev(DHHEB)       * 0.043
 #DHPRNstd       = stdev(DHPRN)       * 0.043
 #writedata %(HEB)0dphsicochemStd.dat DHSysstd DHPstd DHSstd DHNonPolarstd DHAromaticstd DHPolarstd DHAcidicstd DHBasicstd DHHEBstd DHPRNstd 
 #printdata DHSysstd DHPstd DHSstd DHNonPolarstd DHAromaticstd DHPolarstd DHAcidicstd DHBasicstd DHHEBstd DHPRNstd 

 printdata ODHSys RDHSys DHSys ODHP RDHP DHP ODHS RDHS DHS ODHNonPolar RDHNonPolar DHNonPolar ODHAromatic RDHAromatic DHAromatic ODHPolar RDHPolar DHPolar ODHAcidic RDHAcidic DHAcidic ODHBasic RDHBasic DHBasic ODHHEB RDHHEB DHHEB 
 writedata %(HEB)0dphsicochemAvg.dat ODHSys RDHSys DHSys ODHP RDHP DHP ODHS RDHS DHS ODHNonPolar RDHNonPolar DHNonPolar ODHAromatic RDHAromatic DHAromatic ODHPolar RDHPolar DHPolar ODHAcidic RDHAcidic DHAcidic ODHBasic RDHBasic DHBasic ODHHEB RDHHEB DHHEB 
                    """ %{'HEB': HEB}, file=open("%0dredox.in" %(HEM[idx][0]), 'a'))
                    print(" Running LIE calculation for Oxidation of HEB-%0d ..." %(HEM[idx][0]))
                    subprocess.run("cpptraj -i "+str(HEM[idx][0])+"redox.in > "+str(HEM[idx][0])+"redox.log", shell=True)

                else:
                    print(""" o%0d.prmtop and/or r%0d.prmtop are missing and we can not proceed. Something went wrong.""" %(HEM[idx][0], HEM[idx][0]))

        idx += 1

    ODHSys = [0]*(len(HEM))
    RDHSys = [0]*(len(HEM))
    DHSys = [0]*(len(HEM))

    ODHP = [0]*(len(HEM))
    RDHP = [0]*(len(HEM))
    DHP = [0]*(len(HEM))

    ODHS = [0]*(len(HEM))
    RDHS = [0]*(len(HEM))
    DHS = [0]*(len(HEM))

    ODHNonPolar = [0]*(len(HEM))
    RDHNonPolar = [0]*(len(HEM))
    DHNonPolar = [0]*(len(HEM))

    ODHAromatic = [0]*(len(HEM))
    RDHAromatic = [0]*(len(HEM))
    DHAromatic = [0]*(len(HEM))

    ODHPolar = [0]*(len(HEM))
    RDHPolar = [0]*(len(HEM))
    DHPolar = [0]*(len(HEM))

    ODHAcidic = [0]*(len(HEM))
    RDHAcidic = [0]*(len(HEM))
    DHAcidic = [0]*(len(HEM))

    ODHBasic = [0]*(len(HEM))
    RDHBasic = [0]*(len(HEM))
    DHBasic = [0]*(len(HEM))

    ODHHEM = [0]*(len(HEM))
    RDHHEM = [0]*(len(HEM))
    DHHEM = [0]*(len(HEM))

    DG = [0]*(len(HEH)-1)

    for idx in range(len(HEM)):
        print("\n Result for HEB-%0d:" %(HEM[idx][0]))

        with open("%0dphsicochemAvg.dat" %(HEM[idx][0]), 'r') as fp:
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

                    ODHHEM[idx] = float(LinesToRead[0].split()[25])
                    RDHHEM[idx] = float(LinesToRead[0].split()[26])
                    DHHEM[idx] = float(LinesToRead[0].split()[27])
                    print("  Other Hemes: %8.3f %8.3f %8.3f" %(ODHHEM[idx], RDHHEM[idx], DHHEM[idx]))

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

################################################################################################################################################

def AssignCouplingFromGeom_HemeB(Output):

    HisID = []; HebID = []; HEM = []
    with open("BondDefinitions.txt") as f:
        for line in f:
            HisID.append(int(line.strip().split('.')[1]))
            HebID.append(int(line.strip().split('.')[3]))

    for i in range(0, len(HisID), 2):
        HEM.append([HebID[i], HisID[i], HisID[i+1]])
    
    ang = [0]*(len(HEM)-1)
    Hda = [0]*(len(HEM)-1)
    for idx in range(len(HEM)-1):
        print(HEM[idx][0], HEM[idx][1], HEM[idx][2])

        print("""
 parm    %s.prmtop
 trajin  min.rst7
 
 vector Hem%0d  corrplane :%0d&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 vector Hem%0d  corrplane :%0d&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 run
 
 vectormath vec1 Hem%0d vec2 Hem%0d dotangle out hemplaneorient_%0d-%0d.dat
 
 run
 quit
        """ %(Output, HEM[idx][0], HEM[idx][0], HEM[idx+1][0], HEM[idx+1][0], HEM[idx][0], HEM[idx+1][0], HEM[idx][0], HEM[idx+1][0]), file=open("CalcOrientation.in", "w"))

        if (os.path.isfile("hemplaneorient_%0d-%0d.dat" %(HEM[idx][0], HEM[idx+1][0])) == True):
            print(""" 
 Found hemplaneorient_%0d-%0d.dat from a prior execution.
 This prior output will be used for the analysis.""" %(HEM[idx][0], HEM[idx+1][0])) 

        if (os.path.isfile("hemplaneorient_%0d-%0d.dat" %(HEM[idx][0], HEM[idx+1][0])) == False):
            subprocess.run("cpptraj -i CalcOrientation.in > CalcOrientation.log", shell=True)

        lc = 0
        with open("hemplaneorient_%0d-%0d.dat" %(HEM[idx][0], HEM[idx+1][0])) as fp:
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
    for idx in range(len(HEM)-1):

        if (idx == 0):
            print("Hda(HEM-%0d <-> HEM-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" %(HEM[idx][0],  HEM[idx+1][0], ang[idx], Hda[idx]), file=open('Hda.txt', 'w'))
        else:
            print("Hda(HEM-%0d <-> HEM-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" %(HEM[idx][0],  HEM[idx+1][0], ang[idx], Hda[idx]), file=open('Hda.txt', 'a'))

    print(" Done!")

    return Hda

################################################################################################################################################

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

def ComputeDiffusionCoefficient(AvgHemeSpacing):

    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()
    V,D = derrida.VD(ketf, ketb)
    print("  Diffusion constant = %E cm^2/S"  % (D * ((AvgHemeSpacing)**2))) 
    print("Diffusion constant = %E (cm^2/S)" % (D * ((AvgHemeSpacing)**2)), file=open('D.txt', 'w'))
################################################################################################################################################

def ComputeFlux():

    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()

    Jf,Jb = blumberger.flux(ketf, ketb)
    Javg = (Jf+Jb)/2

    print("Forward Flux: %.2E" %(Jf)) 
    print("Reverse Flux: %.2E" %(Jb))
    print("Average forward/backward Flux: %.2E" %(Javg))

    return Jf, Jb, Javg
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
    #subprocess.run("vmd -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)
    subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)

    with open("SubunitLength.txt") as fp:
        SubunitLength = float(fp.readline().strip().split()[5])

    return SubunitLength

################################################################################################################################################

################################################################################################################################################

def MeasureSubunitLength_HemeB():
    HisID = []; HebID = []; HEM = []
    with open("BondDefinitions.txt") as f:
        for line in f:
            HisID.append(int(line.strip().split('.')[1]))
            HebID.append(int(line.strip().split('.')[3]))
 
    for i in range(0, len(HisID), 2):
        HEM.append([HebID[i], HisID[i], HisID[i+1]])

    print("""
 mol new min.pdb
 set output [open  "SubunitLength.txt" w]

 set FirstHeme [[atomselect top "resname HEB and resid %0d and name FE"] get index]
 set LastHeme  [[atomselect top "resname HEB and resid %0d and name FE"] get index]

 set dista  [measure bond [list $FirstHeme $LastHeme]]
 set distcm [expr {$dista * 1E-8}]

 puts $output "Proposed Subunit Length (cm) = $distcm"
 exit
    """ %(HEM[0][0], HEM[len(HEM)-1][0]), file=open("MeasureSubunitLength.tcl", "w"))
    #subprocess.run("vmd -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True) 
    subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)

    with open("SubunitLength.txt") as fp:
        SubunitLength = float(fp.readline().strip().split()[5])

    return SubunitLength

################################################################################################################################################

################################################################################################################################################

def ComputeRedoxCurrent(HemeType):

    r=0.75E-7      #cm
    c=3.00E10      #cm/s
    e=1.602E-19    #c
    kb=1.38E-23    #J/k
    hbar=1.05E-34  #J/s 
    deltabar=3     #cm^-1

    print(" Please provide the following parmaeters: ")
    T    = float(input("  Temperature (K)? "))
    cps    = int(input("  Number of charges per subunit? "))
    ahs  = float(input("  Average heme spacing (cm)? "))
    Ncnt = float(input("  Fewest number of contacts at either protein-electrode interface? "))
    
    DiffusionCoefficient = ComputeDiffusionCoefficient(ahs)
    
    if (HemeType == "c") or (HemeType == "C"):
        SubunitLength = MeasureSubunitLength()
    if (HemeType == "b") or (HemeType == "B"):
        SubunitLength = MeasureSubunitLength_HemeB()

    print(""" 
 The length of a subunit of the cytochrome polymer is needed.
  The subunit length measured between the first and the last heme
  specified in LinearizedHemeSequence.txt is %.2E \n""" %(SubunitLength))

#   lsub = float(input("  Length of subunit (cm)? "))
    lsub = SubunitLength

    lw = float(input("  Length of wire (cm)? "))
    Gexp = float(input("""  Experimental Conductance (S) ? \n   (Enter "0" if not known)         """))

    if (os.path.isfile("D.txt") == True):
        with open("D.txt") as fp:
            Dcalc = float(fp.readline().strip().split()[3])

    cpsul = ((cps)/(lsub))
    csa = (math.pi * (r)**2) 
    crgden = (cpsul)/(csa)

    if (Gexp == "") or (Gexp == 0):
        V = 0.1
    else:
        V = ((1E12)*(e))/Gexp

    k_dif = (Dcalc)/((ahs)**2)
    u_dif = (Dcalc)*(e/(kb*T))
    s_dif = e*(crgden)*(u_dif)
    j_dif = ((s_dif*csa)/(e*Ncnt*lw))*(V)

    s_flx = ((Javg*e*Ncnt*lw)/(V*csa)) 
    u_flx = s_flx/((e)*(crgden))
    d_flx = ((kb*T)/(e))*(u_flx)
    k_flx = (d_flx)/((ahs)**2)

    u_bt = ((e)*(ahs**2))/(2*hbar)
    d_bt = ((kb*T)/(e))*(u_bt)
    k_bt = (d_bt)/((ahs)**2)
    s_bt = (e)*(crgden)*(u_bt)
    j_bt = ((s_bt*csa)/(e*Ncnt*lw))*(V)

    u_hp = (((2)*(math.pi)*(c)*(e)*(ahs**2))/((kb)*(T)))*(deltabar)
    d_hp = ((kb*T)/(e))*(u_hp)
    k_hp = (d_hp)/((ahs)**2)
    s_hp = (e)*(crgden)*(u_hp)
    j_hp = ((s_hp*csa)/(e*Ncnt*lw))*(V)

    if (Gexp != "") or (Gexp != "0"):
        s_exp = Gexp*(lw/csa)
        j_exp = ((s_exp*csa)/(e*Ncnt*lw))*(V)
        u_exp = s_exp/((e)*(crgden))
        d_exp = ((Gexp * kb * T * lw) / (csa * crgden * (e)**2))
        k_exp = (d_exp)/((ahs)**2)

    print("""
 Entered Quantities 
   Temperature                                       = %.1f K    
   Average Heme Spacing                              = %e cm
   Charge per Subunit Length                         = %e q/cm 
   Length of Filament                                = %e cm 
   Fewest number of protein-electrode contacts       = %e cm 
   Experimental Conductance                          = %s S

 Computed Quantities 
   1) Structure Properties:
        Cross-Sectional Area                         = %e cm^2
        Charge Density                               = %e q/cm^2

   2) Single Particle Diffusion Model
        Diffusion Constant                           = %e cm^2/s
        Homogenious Chain Hopping Rate               = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   3) Multi-Particle Steady-State Model
        Diffusion Constant                           = %e cm^2/s
        Homogenious Chain Hopping Rate               = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   4) Band Theory Minimum Requirements
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   5) Hopping Rate Maximum Limits
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at at %.1f V                            = %e electrons/s
    """ %(T, ahs, cps, lw, Ncnt, Gexp, csa, crgden, Dcalc, k_dif, u_dif, s_dif, V, j_dif, d_flx, k_flx, u_flx, s_flx, V, Javg, d_bt, k_bt, u_bt, s_bt, V, j_bt, d_hp, k_hp, u_hp, s_hp, V, j_hp)) 

    if (Gexp != 0):
        print("""
   6) Experiment-Based Quantities
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s
        """ %(d_exp, k_exp, u_exp, s_exp, V, j_exp)) 

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

    PDB = Initialization()

    while True:
        ChooseMut = input("\n Would you like to mutate a residue? ")

        if (ChooseMut == "y") or (ChooseMut == "Y") or (ChooseMut == "yes") or (ChooseMut == "Yes") or (ChooseMut == "YES"):
            PDB = Mutate(PDB)
            break
        elif (ChooseMut == "n") or (ChooseMut == "N") or (ChooseMut == "no") or (ChooseMut == "No") or (ChooseMut == "NO"):
            print(" Structure will not be mutated.")
            break
        else:
            print(" Sorry, I didn't understand your response.")


    print("""
 Different parameters need to be applied for b- and c-type hemes. 
 Presently, only structures with all b- or all c-type hemes can be analyzed.""")

    while True:
        HemeType = input("\n Which type (b or c) does your structure contain? ")

        if (HemeType == "c") or (HemeType == "C"):
            CreateResIndexing(PDB)
            ProcessPDB(PDB)
            Output, Solv = ReBuildStructure(PDB)
            StructRelax(Output, Solv)
            break
        elif (HemeType == "b") or (HemeType == "B"):
            Output, Solv = SetUpHemeB(PDB)
            StructRelax(Output, Solv)
            break
        else: 
            print("""\n Sorry, I didn't understand your selection for the type of solvent used.""")

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

 To oblige this different numbering convention, we'll use 
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

    print("""
 We also need to know if the hemes are of the b- or c-type
 to compute the energetics correctly.
    """)

    HemeType = input(" Are there all b- or all c-type hemes? ")

    print(" ------------------------------------------------------------------- ")

    while True:
        CompLambda = input("""
 Should we compute the reorganization energy (yes/no)? """)

        if (CompLambda == 'Yes') or (CompLambda == "yes") or (CompLambda == "Y") or (CompLambda == "y"):
            if (HemeType == "c") or (HemeType == "C"):
                Lambda = LambdaFromSASA(Output)
                break

            if (HemeType == "b") or (HemeType == "B"):
                Lambda = LambdaFromSASA_HemeB(Output)
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

    Recommendation:
        The LIE approach has the advantages that it is considerably 
        faster than the PBSA approach, and the overall change in 
        electorstatic energy is decomposed into contributions from 
        different groups of residues. 
	
	However, the LIE approach is only recommended for the analysis
	of multiple snapshots. Electrostatic interactions need to be 
	thermally averaged to reduce the strong dependence of fixed
	charge (non-polarizable) interactions on a given configuration.
	This averaging is effectively done in the PBSA approach, at 
	least for the solvent.

	PBSA should be used to analyse a single or multiple configurations.
	LIE should be used to analyze only multiple configurations.

	Both methods should give better results as the number of examiined
	configurations increases. 
                """)

                while True:
                    DGmethod = input(""" Should we use the LIE or PBSA method (lie/pbsa)? """)

                    if (DGmethod == "LIE") or (DGmethod == "lie") or (DGmethod == str(1)):
                        if (HemeType == "c") or (HemeType == "C"):
                            DG = DeltaGFromLIE(Output)
                        if (HemeType == "b") or (HemeType == "B"):
                            DG = DeltaGFromLIE_HemeB(Output)
                        break

                    elif (DGmethod == "PBSA") or (DGmethod == "pbsa") or (DGmethod == str(2)):
                        if (HemeType == "c") or (HemeType == "C"):
                            DG = DeltaGFromPBSA(Output, Solv)
                        if (HemeType == "b") or (HemeType == "B"):
                            DG = DeltaGFromPBSA_HemeB(Output, Solv)
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
                if (HemeType == "c") or (HemeType == "C"):
                	DG = DeltaGFromPBSA(Output, Solv)
                if (HemeType == "b") or (HemeType == "B"):
                	DG = DeltaGFromPBSA_HemeB(Output, Solv)
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
        CompIntEng = input("""
 Should we compute heme-heme interaction energies (yes/no)? """)

        if (CompIntEng == 'Yes') or (CompIntEng == "yes") or (CompIntEng == "Y") or (CompIntEng == "y"):
            if (HemeType == "c") or (HemeType == "C"):
                IntEng = HemeHemeInt(Output, Solv)
            if (HemeType == "b") or (HemeType == "B"):
                IntEng = HemeHemeInt_HemeB(Output, Solv)
            break 
        elif (CompIntEng == 'No') or (CompIntEng == "no") or (CompIntEng == "N") or (CompIntEng == "n"):
            print("Skipping the computation of heme-heme interaction energies.")
        else: 
            print("Sorry, I didn't understand your response.")

    while True:
        CompHda = input("""
 Should we estimate the electronic coupling from the geometry (yes/no)? """)

        if (CompHda == 'Yes') or (CompHda == "yes") or (CompHda == "Y") or (CompHda == "y"):
            if (HemeType == "c") or (HemeType == "C"):
                Hda = AssignCouplingFromGeom(Output)
            if (HemeType == "b") or (HemeType == "B"):
                Hda = AssignCouplingFromGeom_HemeB(Output)
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
 expression derived by Blumberger and co-workers (Refs. #5-#7).
 The original code was kindly provided by Jochen Blumberger 
 and Xiuyun Jiang.

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

    HemeType = input("\n Are there all b- or all c-type hemes? ")

    if (os.path.isfile("rates.txt") == True):
        print("""
 Found rates.txt, which is needed to proceed!

 We will now compute the multi-particle, stead-state flux.
 """)
        Jf,Jb,Javg = ComputeFlux()

#       print("""
#We will now compute the single-particle diffusion coefficient. 
#""")
#       ComputeDiffusionCoefficient()
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

    if (os.path.isfile("rates.txt") == True):
        print("""
 We will at last compute the redox current. 
 To do this, some system-specific information is needed.
        """)
        ComputeRedoxCurrent(HemeType)
        print(" Done!")
    else:
        sys.exit("""
 D.txt not found! Please re-run division #3 of BioDC or
 create the file by hand. In the latter case, the syntax
 should be:
    Diffusion constant = [D] cm^2/S
 where "[D]" is to be replaced with the diffusion coefficient.
        """)



