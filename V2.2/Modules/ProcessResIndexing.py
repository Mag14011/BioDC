################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################
def ProcessResIndexing(PDB, DisulfList, SelASPIDs, SelGLUIDs, SelHISIDs, SelLYSIDs, SelTYRIDs, SelPRNIDs, InputDict, LaunchDir):

    if (os.path.isfile(f"ResIndexing.txt") == False):
        sys.exit("""
 ResIndexing.txt was not found! 
 I, unfortunately, do not know how to proceed without this file.\n""")

    elif (os.path.isfile(f"ResIndexing.txt") == True):
        print("""
 Your ResIndexing file will now be used with VMD
 to re-label the residues according to the available
 AMBER FF10 parameterizations of b- and c-type, 
 His-His or His-Met ligated hemes. 

 Also, if you specified any disulfide residues 
 earlier, the participating Cys residues will 
 be re-labeled as CYX, in accordance with AMBER
 conventions. 

 We will write and submit a script called 
 SetupStructure.tcl to perform this magic.

 CAUTION: The magic of the script is only 
 as good as the information in ResIndexing.txt. 
 If the wrong residue IDs are specified, everything 
 from here on out will be, put politely, junk!
    """, end=" ")

        while True:
            if ("ProcessPDBChoice" in InputDict):
                ProcessPDBChoice = InputDict["ProcessPDBChoice"]
            else:
                ProcessPDBChoice = input(""" 
 Shall we venture forward with SetupStructure.tcl (yes/no)? """)

            print(f"ProcessPDBChoice = {ProcessPDBChoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if (ProcessPDBChoice == "YES") or (ProcessPDBChoice == "Yes") or (ProcessPDBChoice == "yes") or (ProcessPDBChoice == "Y") or (ProcessPDBChoice == "y"):
                print(f"""
 #-------------------------------------------------------------------------
 #Input
   mol load pdb {PDB}_renumd.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'w'))

                if (len(DisulfList) != ""):
                    print(f"""
 #-------------------------------------------------------------------------
 #Change residue names of disulfide-linked Cys residues
   set DisulfCys [atomselect top "resname CYS CYX and resid {DisulfList}"] 
   $DisulfCys set resname CYX
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a'))

                if (len(SelASPIDs) != "") or (len(SelGLUIDs) != "") or (len(SelHISIDs) != ""):
                    print(f"""
 #-------------------------------------------------------------------------
 #Change residue names of selected titratable residues""", file=open('SetupStructure.tcl', 'a'))

                if (len(SelASPIDs) != ""):
                    print(f"""
   set ASP [atomselect top "resname ASP and resid {SelASPIDs}"] 
   $ASP set resname AS4 """, file=open('SetupStructure.tcl', 'a'))

                if (len(SelGLUIDs) != ""):
                    print(f"""
   set GLU [atomselect top "resname GLU and resid {SelGLUIDs}"] 
   $GLU set resname GL4 """, file=open('SetupStructure.tcl', 'a'))

                if (len(SelHISIDs) != ""):
                    print(f"""
   set HIS [atomselect top "resname HIS HID HIE HIP and resid {SelHISIDs}"] 
   $HIS set resname HIP """, file=open('SetupStructure.tcl', 'a'))

                idx=1
                LineNumErr = ''
                print(f"RedoxState =", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                with open("ResIndexing.txt") as fp:
                    NumHEC = len(fp.readlines())
                    print(f" \n There are {NumHEC} hemes")
                    fp.seek(0)
 
                    Lines = fp.readlines()
                    for line in Lines:
                        EntryLength = len(line.strip().split(" "))
                        HemeType = line.strip().split(" ")[-2]
                        AxLigType = line.strip().split(" ")[-1]

                        if ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HH" ):
                            CYSb = line.strip().split(" ")[0]
                            CYSc = line.strip().split(" ")[1]
                            HISp = line.strip().split(" ")[2]
                            HISd = line.strip().split(" ")[3]
                            HECHHshft = line.strip().split(" ")[4]
                            HECHHfnl = line.strip().split(" ")[5]

                            print(f""" 
 Heme-{HECHHfnl} is a bis-His-ligated c-type heme.
   Residue ID of Thioether-linked Cys to heme b-ring: {CYSb}
   Residue ID of Thioether-linked Cys to heme c-ring: {CYSc}
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal His to   Fe-center: {HISd}""")

                            while True:
                                if ("RedoxState" in InputDict):
                                    RedoxState = InputDict["RedoxState"][idx-1]
                                else:
                                    RedoxState = input(" Would you like to model this heme as oxidized or reduced (O/R)? ")

                                print(f"{RedoxState}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                                if (RedoxState == "OX") or (RedoxState == "Ox") or (RedoxState == "ox") or (RedoxState == "O") or (RedoxState == "o"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HECHHfnl}: 

 #Define atom groups
   set CYSbb [atomselect top "resid {CYSb} {CYSc} and name N CA C O"]
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set HISd  [atomselect top "resname HIS HIE HID HIP and resid {HISd}"]
   set HEM   [atomselect top "(resname HEC HEM and resid {HECHHshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D) or (resname CYS and resid {CYSb} {CYSc} and name CB SG)"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HECHHshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HECHHshft} and name CAD CBD CGD O1D O2D"]

 #Change the names of the Cys sidechain atoms bonded to the heme
   set CBB2 [atomselect top "resname CYS and resid {CYSb} and name CB"]
   $CBB2 set name CBB2
   set SGB2 [atomselect top "resname CYS and resid {CYSb} and name SG"]
   $SGB2 set name SGB2
   set CBC1 [atomselect top "resname CYS and resid {CYSc} and name CB"]
   $CBC1 set name CBC1
   set SGC1 [atomselect top "resname CYS and resid {CYSc} and name SG"]
   $SGC1 set name SGC1

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CAA) or (resname HEC HEM and resid {HECHHshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CBA) or (resname HEC HEM and resid {HECHHshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CGA) or (resname HEC HEM and resid {HECHHshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name O1A) or (resname HEC HEM and resid {HECHHshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name O2A) or (resname HEC HEM and resid {HECHHshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The backbone of the bonded Cys residues
     $CYSbb set resname CYO
   #The proximal His residue 
     $HISp  set resname PHO; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHO; #Distal   His for oxidized His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCO; #c-type His-His ligated oxidized heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 
 
 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HECHHfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HECHHfnl)+1}
     $PRND set resid {int(HECHHfnl)+2} 

 #Write PDBs
   $HEM  writepdb HCO{HECHHfnl}.pdb
   $PRNA writepdb PRN{int(HECHHfnl)+1}.pdb
   $PRND writepdb PRN{int(HECHHfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break

                                if (RedoxState == "RED") or (RedoxState == "Red") or (RedoxState == "red") or (RedoxState == "R") or (RedoxState == "r"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HECHHfnl}: 

 #Define atom groups
   set CYSbb [atomselect top "resid {CYSb} {CYSc} and name N CA C O"]
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set HISd  [atomselect top "resname HIS HIE HID HIP and resid {HISd}"]
   set HEM   [atomselect top "(resname HEC HEM and resid {HECHHshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D) or (resname CYS and resid {CYSb} {CYSc} and name CB SG)"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HECHHshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HECHHshft} and name CAD CBD CGD O1D O2D"]

 #Change the name of the Cys sidechain atoms bonded to the heme
   set CBB2 [atomselect top "resname CYS and resid {CYSb} and name CB"]
   $CBB2 set name CBB2
   set SGB2 [atomselect top "resname CYS and resid {CYSb} and name SG"]
   $SGB2 set name SGB2
   set CBC1 [atomselect top "resname CYS and resid {CYSc} and name CB"]
   $CBC1 set name CBC1
   set SGC1 [atomselect top "resname CYS and resid {CYSc} and name SG"]
   $SGC1 set name SGC1

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CAA) or (resname HEC HEM and resid {HECHHshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CBA) or (resname HEC HEM and resid {HECHHshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name CGA) or (resname HEC HEM and resid {HECHHshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name O1A) or (resname HEC HEM and resid {HECHHshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HECHHshft} and name O2A) or (resname HEC HEM and resid {HECHHshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The backbone of the bonded Cys residues
     $CYSbb set resname CYO
   #The proximal His residue 
     $HISp  set resname PHR; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHR; #Distal   His for oxidized His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCR; #c-type His-His ligated oxidized heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 
 
 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HECHHfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HECHHfnl)+1}
     $PRND set resid {int(HECHHfnl)+2} 

 #Write PDBs
   $HEM  writepdb HCR{HECHHfnl}.pdb
   $PRNA writepdb PRN{int(HECHHfnl)+1}.pdb
   $PRND writepdb PRN{int(HECHHfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break
                                else:
                                    print(" Sorry, I didn't understand your response.")

                        elif ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HM" ):
                            CYSb = line.strip().split(" ")[0]
                            CYSc = line.strip().split(" ")[1]
                            HISp = line.strip().split(" ")[2]
                            METd = line.strip().split(" ")[3]
                            HECHMshft = line.strip().split(" ")[4]
                            HECHMfnl = line.strip().split(" ")[5]

                            print(f""" 
 Heme-{HECHMfnl} is a His-Met-ligated c-type heme.
   Residue ID of Thioether-linked Cys to heme b-ring: {CYSb}
   Residue ID of Thioether-linked Cys to heme c-ring: {CYSc}
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal Met to   Fe-center: {METd}""")

                            while True:
                                if ("RedoxState" in InputDict):
                                    RedoxState = InputDict["RedoxState"][idx-1]
                                else:
                                    RedoxState = input(" Would you like to model this heme as oxidized or reduced (O/R)? ")

                                print(f"{RedoxState}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                                if (RedoxState == "OX") or (RedoxState == "Ox") or (RedoxState == "ox") or (RedoxState == "O") or (RedoxState == "o"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HECHMfnl}: 
 
 #Define atom groups
   set CYSbb [atomselect top "resid {CYSb} {CYSc} and name N CA C O"]
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set METd  [atomselect top "resname MET and resid {METd}"]
   set HEM   [atomselect top "(resname HEC HEM and resid {HECHMshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D) or (resname CYS and resid {CYSb} {CYSc} and name CB SG)"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HECHMshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HECHMshft} and name CAD CBD CGD O1D O2D"]

 #Change the name of the Cys sidechain atoms bonded to the heme
   set CBB2 [atomselect top "resname CYS and resid {CYSb} and name CB"]
   $CBB2 set name CBB2
   set SGB2 [atomselect top "resname CYS and resid {CYSb} and name SG"]
   $SGB2 set name SGB2
   set CBC1 [atomselect top "resname CYS and resid {CYSc} and name CB"]
   $CBC1 set name CBC1
   set SGC1 [atomselect top "resname CYS and resid {CYSc} and name SG"]
   $SGC1 set name SGC1

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CAA) or (resname HEC HEM and resid {HECHMshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CBA) or (resname HEC HEM and resid {HECHMshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CGA) or (resname HEC HEM and resid {HECHMshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name O1A) or (resname HEC HEM and resid {HECHMshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name O2A) or (resname HEC HEM and resid {HECHMshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The backbone of the bonded Cys residues
     $CYSbb set resname CYO
   #The proximal His residue 
     $HISp  set resname PMO; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $METd  set resname DMO; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCO; #c-type His-Met ligated oxidized heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 
 
 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HECHMfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HECHMfnl)+1}
     $PRND set resid {int(HECHMfnl)+2} 

 #Write PDBs
   $HEM  writepdb MCO{HECHMfnl}.pdb
   $PRNA writepdb PRN{int(HECHMfnl)+1}.pdb
   $PRND writepdb PRN{int(HECHMfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break

                                if (RedoxState == "RED") or (RedoxState == "Red") or (RedoxState == "red") or (RedoxState == "R") or (RedoxState == "r"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HECHMfnl}: 

 #Define atom groups
   set CYSbb [atomselect top "resid {CYSb} {CYSc} and name N CA C O"]
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set METd  [atomselect top "resname MET and resid {METd}"]
   set HEM   [atomselect top "(resname HEC HEM and resid {HECHMshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D) or (resname CYS and resid {CYSb} {CYSc} and name CB SG)"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HECHMshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HECHMshft} and name CAD CBD CGD O1D O2D"]

 #Change the name of the Cys sidechain atoms bonded to the heme
   set CBB2 [atomselect top "resname CYS and resid {CYSb} and name CB"]
   $CBB2 set name CBB2
   set SGB2 [atomselect top "resname CYS and resid {CYSb} and name SG"]
   $SGB2 set name SGB2
   set CBC1 [atomselect top "resname CYS and resid {CYSc} and name CB"]
   $CBC1 set name CBC1
   set SGC1 [atomselect top "resname CYS and resid {CYSc} and name SG"]
   $SGC1 set name SGC1

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CAA) or (resname HEC HEM and resid {HECHMshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CBA) or (resname HEC HEM and resid {HECHMshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name CGA) or (resname HEC HEM and resid {HECHMshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name O1A) or (resname HEC HEM and resid {HECHMshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HECHMshft} and name O2A) or (resname HEC HEM and resid {HECHMshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The backbone of the bonded Cys residues
     $CYSbb set resname CYO
   #The proximal His residue 
     $HISp  set resname PMR; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $METd  set resname DMR; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCR; #c-type His-Met ligated oxidized heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 
 
 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HECHMfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HECHMfnl)+1}
     $PRND set resid {int(HECHMfnl)+2} 

 #Write PDBs
   $HEM  writepdb MCR{HECHMfnl}.pdb
   $PRNA writepdb PRN{int(HECHMfnl)+1}.pdb
   $PRND writepdb PRN{int(HECHMfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break
                                else:
                                    print(" Sorry, I didn't understand your response.")


                        elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HH" ):
                            HISp = line.strip().split(" ")[0]
                            HISd = line.strip().split(" ")[1]
                            HEBHHshft = line.strip().split(" ")[2]
                            HEBHHfnl = line.strip().split(" ")[3]
                            print(f""" 
 Heme-{HEBHHfnl} is a His-His-ligated b-type heme.
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal His to   Fe-center: {HISd}""")

                            while True:
                                if ("RedoxState" in InputDict):
                                    RedoxState = InputDict["RedoxState"][idx-1]
                                else:
                                    RedoxState = input(" Would you like to model this heme as oxidized or reduced (O/R)? ")

                                print(f"{RedoxState}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                                if (RedoxState == "OX") or (RedoxState == "Ox") or (RedoxState == "ox") or (RedoxState == "O") or (RedoxState == "o"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HEBHHfnl}: 
 
 #Define atom groups
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set HISd  [atomselect top "resname HIS HIE HID HIP and resid {HISd}"]
   set HEM   [atomselect top "resname HEC HEM and resid {HEBHHshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HEBHHshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HEBHHshft} and name CAD CBD CGD O1D O2D"]

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CAA) or (resname HEC HEM and resid {HEBHHshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CBA) or (resname HEC HEM and resid {HEBHHshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CGA) or (resname HEC HEM and resid {HEBHHshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name O1A) or (resname HEC HEM and resid {HEBHHshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name O2A) or (resname HEC HEM and resid {HEBHHshft} and name O2D)"]
    $O2 set name O2


 #Change residue names for:
   #The proximal His residue 
     $HISp set resname FHO; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHO; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBO; #b-type His-His ligated heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 

 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HEBHHfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HEBHHfnl)+1}
     $PRND set resid {int(HEBHHfnl)+2}

 #Write PDBs
   $HEM  writepdb HBO{HEBHHfnl}.pdb
   $PRNA writepdb PRN{int(HEBHHfnl)+1}.pdb
   $PRND writepdb PRN{int(HEBHHfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break

                                if (RedoxState == "RED") or (RedoxState == "Red") or (RedoxState == "red") or (RedoxState == "R") or (RedoxState == "r"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HEBHHfnl}: 

 #Define atom groups
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set HISd  [atomselect top "resname HIS HIE HID HIP and resid {HISd}"]
   set HEM   [atomselect top "resname HEC HEM and resid {HEBHHshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HEBHHshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HEBHHshft} and name CAD CBD CGD O1D O2D"]

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CAA) or (resname HEC HEM and resid {HEBHHshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CBA) or (resname HEC HEM and resid {HEBHHshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name CGA) or (resname HEC HEM and resid {HEBHHshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name O1A) or (resname HEC HEM and resid {HEBHHshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HEBHHshft} and name O2A) or (resname HEC HEM and resid {HEBHHshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The proximal His residue 
     $HISp set resname FHR; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHR; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBR; #b-type His-His ligated heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 

 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HEBHHfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HEBHHfnl)+1}
     $PRND set resid {int(HEBHHfnl)+2}

 #Write PDBs
   $HEM  writepdb HBR{HEBHHfnl}.pdb
   $PRNA writepdb PRN{int(HEBHHfnl)+1}.pdb
   $PRND writepdb PRN{int(HEBHHfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break
                                else:
                                    print(" Sorry, I didn't understand your response.")

                        elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HM" ):
                            HISp = line.strip().split(" ")[0]
                            METd = line.strip().split(" ")[1]
                            HEBHMshft = line.strip().split(" ")[2]
                            HEBHMfnl = line.strip().split(" ")[3]
                            print(f""" 
 Heme-{HEBHMfnl} is a His-Met-ligated b-type heme.
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal Met to   Fe-center: {METd}""")

                            while True:
                                if ("RedoxState" in InputDict):
                                    RedoxState = InputDict["RedoxState"][idx-1]
                                else:
                                    RedoxState = input(" Would you like to model this heme as oxidized or reduced (O/R)? ")

                                print(f"{RedoxState}", end=" ", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                                if (RedoxState == "OX") or (RedoxState == "Ox") or (RedoxState == "ox") or (RedoxState == "O") or (RedoxState == "o"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HEBHMfnl}: 

 #Define atom groups
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set METd  [atomselect top "resname MET and resid {METd}"]
   set HEM   [atomselect top "resname HEC HEM and resid {HEBHMshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HEBHMshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HEBHMshft} and name CAD CBD CGD O1D O2D"]

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CAA) or (resname HEC HEM and resid {HEBHMshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CBA) or (resname HEC HEM and resid {HEBHMshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CGA) or (resname HEC HEM and resid {HEBHMshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name O1A) or (resname HEC HEM and resid {HEBHMshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name O2A) or (resname HEC HEM and resid {HEBHMshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The proximal His residue 
     $HISp set resname FMO; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $METd set resname RMO; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname MBO; #b-type His-His ligated heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 

 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HEBHMfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HEBHMfnl)+1}
     $PRND set resid {int(HEBHMfnl)+2}

 #Write PDBs
   $HEM  writepdb MBO{HEBHMfnl}.pdb
   $PRNA writepdb PRN{int(HEBHMfnl)+1}.pdb
   $PRND writepdb PRN{int(HEBHMfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break

                                if (RedoxState == "RED") or (RedoxState == "Red") or (RedoxState == "red") or (RedoxState == "R") or (RedoxState == "r"):
                                    print(f"""
 #-------------------------------------------------------------------------
 #Heme-{HEBHMfnl}: 

 #Define atom groups
   set HISp  [atomselect top "resname HIS HIE HID HIP and resid {HISp}"]
   set METd  [atomselect top "resname MET and resid {METd}"]
   set HEM   [atomselect top "resname HEC HEM and resid {HEBHMshft} and not name CAA CAD CBA CBD CGA CGD O1A O1D O2A O2D"]
   set PRNA  [atomselect top "resname HEC HEM and resid {HEBHMshft} and name CAA CBA CGA O1A O2A"]
   set PRND  [atomselect top "resname HEC HEM and resid {HEBHMshft} and name CAD CBD CGD O1D O2D"]

 #Change the names of the propionic acid atoms connected to the heme
    set CA    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CAA) or (resname HEC HEM and resid {HEBHMshft} and name CAD)"]
    $CA set name CA
    set CB    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CBA) or (resname HEC HEM and resid {HEBHMshft} and name CBD)"]
    $CB set name CB
    set CG    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name CGA) or (resname HEC HEM and resid {HEBHMshft} and name CGD)"]
    $CG set name CG
    set O1    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name O1A) or (resname HEC HEM and resid {HEBHMshft} and name O1D)"]
    $O1 set name O1
    set O2    [atomselect top "(resname HEC HEM and resid {HEBHMshft} and name O2A) or (resname HEC HEM and resid {HEBHMshft} and name O2D)"]
    $O2 set name O2

 #Change residue names for:
   #The proximal His residue 
     $HISp set resname FMR; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $METd set resname RMR; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname MBR; #b-type His-His ligated heme
   #The propionic acid groups that eere split-off from the heme into separate residues
     $PRNA set resname PRN 
     $PRND set resname PRN 

 #Change the residue IDs for:
   #The heme group
     $HEM  set resid {HEBHMfnl} 
   #The two propionic acid groups bonded to the heme
     $PRNA set resid {int(HEBHMfnl)+1}
     $PRND set resid {int(HEBHMfnl)+2}

 #Write PDBs
   $HEM  writepdb MBR{HEBHMfnl}.pdb
   $PRNA writepdb PRN{int(HEBHMfnl)+1}.pdb
   $PRND writepdb PRN{int(HEBHMfnl)+2}.pdb
 #-------------------------------------------------------------------------""", file=open('SetupStructure.tcl', 'a')) 
                                    break
                                else:
                                    print(" Sorry, I didn't understand your response.")
                        else:
                            LineNumErr = str(LineNumErr)+str(f" {idx}")

                        idx+=1

                print(f"\n", end="", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                print(f"""
 #-------------------------------------------------------------------------
 #Write PDB of protein
   set prot [atomselect top "all and (noh) and (not resname HCO HCR MCO MCR HBO HBR MBO MBR PRN) and (not water) and (not ions)"]
   $prot writepdb prot.pdb
 #-------------------------------------------------------------------------

 exit """, file=open('SetupStructure.tcl', 'a'))

                if ( LineNumErr != '' ):
                    NumErr = len(LineNumErr.strip().split(" "))
                    if ( NumErr == 1 ):
                        sys.exit(f""" 
 Something is wrong in ResIndexing.txt on line: {LineNumErr}. 
 Please check and correct the file before re-running BioDC.
 Please save the corrected file as CorrectedResIndexing.txt.
 BioDC will automatically detected this file and use it when
 you re-run module #1.""")
                    elif ( NumErr > 1 ):
                        sys.exit(f""" 
 Something is wrong in ResIndexing.txt on lines: {LineNumErr}. 
 Please check and correct the file before re-running BioDC.
 Please save the corrected file as CorrectedResIndexing.txt.
 BioDC will automatically detected this file and use it when
 you re-run module #1.""")
                elif ( LineNumErr == '' ):
                    subprocess.run("vmd -e SetupStructure.tcl > SetupStructure.log", shell=True)
                    print("""
 VMD finished. Please check SetupStructure.log for any errors. You may 
 also want to inspect the generated PDB.
            """, end=" ")
                    break

            elif (ProcessPDBChoice == "NO") or (ProcessPDBChoice == "No") or (ProcessPDBChoice == "no") or (ProcessPDBChoice == "N") or (ProcessPDBChoice == "n"):
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
