################################################################################################################################################
# Generic Modules
import os
import sys
import itertools
import subprocess
from subprocess import Popen
################################################################################################################################################

def insert_ter_records(pdb_file, output_file):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    modified_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Check for N followed by H1, H2, H3
        if line.startswith("ATOM") and line[12:16].strip() == "N":
            # Check if the previous line starts with CRYST1
            if i > 0 and not lines[i-1].startswith("CRYST1"):
                # Check for TER before N
                if not (modified_lines and modified_lines[-1].startswith("TER")):
                    if (i + 3 < len(lines) and
                        lines[i+1].startswith("ATOM") and lines[i+1][12:16].strip() == "H1" and
                        lines[i+2].startswith("ATOM") and lines[i+2][12:16].strip() == "H2" and
                        lines[i+3].startswith("ATOM") and lines[i+3][12:16].strip() == "H3"):
                        # Insert TER if the sequence N, H1, H2, H3 is found
                        modified_lines.append("TER\n")

        # Check for N followed by H2, H3
        if line.startswith("ATOM") and line[12:16].strip() == "N":
            # Check if the previous line starts with CRYST1
            if i > 0 and not lines[i-1].startswith("CRYST1"):
                # Check for TER before N
                if not (modified_lines and modified_lines[-1].startswith("TER")):
                    if (i + 2 < len(lines) and
                        lines[i+1].startswith("ATOM") and lines[i+1][12:16].strip() == "H2" and
                        lines[i+2].startswith("ATOM") and lines[i+2][12:16].strip() == "H3"):
                        # Insert TER if the sequence N, H2, H3 is found
                        modified_lines.append("TER\n")

        modified_lines.append(line)

        # Ensure there is a TER record after OXT
        if line.startswith("ATOM") and line[12:16].strip() == "OXT":
            if (i + 1 >= len(lines) or not lines[i + 1].startswith("TER")):
                modified_lines.append("TER\n")

        i += 1

    with open(output_file, 'w') as file:
        file.writelines(modified_lines)

def PairedChargeAssignment(ForceFieldDir, FFchoice, SelRefRedoxState, InputDict):

#------------------------------------------------------------------------------
# Count how many of each type of heme are present.
#------------------------------------------------------------------------------
    Count_c_HH = 0
    Count_c_HM = 0
    Count_b_HH = 0
    Count_b_HM = 0

    with open("ResIndexing.txt") as fp:
        NumHEC = len(fp.readlines())
        HemID = [0]*NumHEC
        fp.seek(0)

        Lines = fp.readlines()
        for line in Lines:
            EntryLength = len(line.strip().split(" "))
            HemeType = line.strip().split(" ")[-2]
            AxLigType = line.strip().split(" ")[-1]

            if ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HH" ):
                Count_c_HH+=1
            elif ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HM" ):
                Count_c_HM+=1
            elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HH" ):
                Count_b_HH+=1
            elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HM" ):
                Count_b_HM+=1

    if ( len(HemID) != (Count_c_HH + Count_c_HM + Count_b_HH +Count_b_HM)):
        sys.exit(f"""
 The total number of hemes ({len(HemID)}) in ResIndexing.txt 
 does NOT equl the sum of c-type His-His ({Count_c_HH}), 
 b-type His-His ({Count_b_HH}), c-type His-Met ({Count_c_HM}),
 and b-type His-Met ({Count_b_HM}) hemes. These are the only
 types of hemes that can currently be analyzed with BioDC.
 Please revise ResIndexing.txt before re-running BioDC.""")

#------------------------------------------------------------------------------
# Determine the unique pairwise heme combinations among the selected hemes.
# The selected hemes are those specified as the linear sequence to 
# LinearizeHemeSequence.py. 
#------------------------------------------------------------------------------
    idx=0
    if (os.path.isfile("SelResIndexing.txt") == True):
        with open("SelResIndexing.txt") as fp:
            NumHEC = len(fp.readlines())
            print(f"   There are {NumHEC} hemes")
            HEM = [0]*NumHEC
            ActiveIDs = [0]*NumHEC
            fp.seek(0)

            Lines = fp.readlines()
            for line in Lines:
                HEM[idx] = int(line.strip().split(" ")[-3])
                ActiveIDs[idx] = line
                idx+=1
        PairCount = list(itertools.combinations(HEM, r=2))
        print(f"   There is/are {len(PairCount)} unique heme-heme pairwise interactions")

    elif (os.path.isfile("SelResIndexing.txt") == False):
        sys.exit("""
SelResIndexing.txt is missing.
Something went wrong when you defined
the linear sequence of hemes.

This problem must be resolved before 
proceeding.""")

#------------------------------------------------------------------------------
# Loop over the unique heme pairs, writing a TCL script for VMD and an input
# file for TLEaP for each pair. The VMD script changes residue names to be 
# consistent with the redox microstate for the heme pair (OO, OR, RO, and RR, 
# where O = oxidized heme, R = reduced heme.) The generated four PDBs are read
# into TLEaP, which builds the corresponding topologies (prmtop) and coordinate
# (rst7) files. 
#------------------------------------------------------------------------------
    for i in range(len(PairCount)):
        Hi = PairCount[i][0] 
        Hj = PairCount[i][1]
        CYSb = [0]*2
        CYSc = [0]*2
        Ligp = [0]*2
        Ligd = [0]*2

        VMDinput=f"PairInt_{Hi}-{Hj}.tcl"
        TLEaPinput=f"GeneratePairIntTopologiesForHems{Hi}-{Hj}.in" 

        print(f"""
 ------------------------------------------------
 Topologies for Heme Pair {Hi} - {Hj} 
 ------------------------------------------------""")

        print("""
 mol new RefState.pdb""", file=open(VMDinput, 'w'))

        print("""
# Load parameters
 source leaprc.constph
 source leaprc.conste
 source leaprc.gaff
 source leaprc.water.tip3p

 addAtomTypes {""", end=' ', file=open(TLEaPinput, 'w'))

#------------------------------------------------------------------------------
# Identify each type of heme in the pair.
#------------------------------------------------------------------------------
        for m in range(len(ActiveIDs)):
            if (Hi == int(ActiveIDs[m].strip().split(" ")[-3])):
                HiLen = len(ActiveIDs[m].strip().split(" "))
                HiType = ActiveIDs[m].strip().split(" ")[-2]
                HiAxLigType = ActiveIDs[m].strip().split(" ")[-1]
                print(f"****** HiAxLigType = {HiAxLigType}")

                if ( HiLen == 8 ) and ( HiType == "c") and ( HiAxLigType == "HH" ):
                    SelHemIType = "cHH"
                    CYSb[0] = int(ActiveIDs[m].strip().split(" ")[0])
                    CYSc[0] = int(ActiveIDs[m].strip().split(" ")[1])
                    Ligp[0] = int(ActiveIDs[m].strip().split(" ")[2])
                    Ligd[0] = int(ActiveIDs[m].strip().split(" ")[3])
                elif ( HiLen == 8 ) and ( HiType == "c") and ( HiAxLigType == "HM" ):
                    SelHemIType = "cHM"
                    CYSb[0] = int(ActiveIDs[m].strip().split(" ")[0])
                    CYSc[0] = int(ActiveIDs[m].strip().split(" ")[1])
                    Ligp[0] = int(ActiveIDs[m].strip().split(" ")[2])
                    Ligd[0] = int(ActiveIDs[m].strip().split(" ")[3])
                elif ( HiLen == 6 ) and ( HiType == "b") and ( HiAxLigType == "HH" ):
                    SelHemIType = "bHH"
                    Ligp[0] = int(ActiveIDs[m].strip().split(" ")[0])
                    Ligd[0] = int(ActiveIDs[m].strip().split(" ")[1])
                elif ( HiLen == 6 ) and ( HiType == "b") and ( HiAxLigType == "HM" ):
                    SelHemIType = "bHM"
                    Ligp[0] = int(ActiveIDs[m].strip().split(" ")[0])
                    Ligd[0] = int(ActiveIDs[m].strip().split(" ")[1])
                else:
                    sys.exit(f""" *** Missing entries on line number {m} of SelResIndexing.txt!""")

            if (Hj == int(ActiveIDs[m].strip().split(" ")[-3])):
                HjLen = len(ActiveIDs[m].strip().split(" "))
                HjType = ActiveIDs[m].strip().split(" ")[-2]
                HjAxLigType = ActiveIDs[m].strip().split(" ")[-1]
                print(f"****** HjAxLigType = {HjAxLigType}")

                if ( HjLen == 8 ) and ( HjType == "c") and ( HjAxLigType == "HH" ):
                    SelHemJType = "cHH"
                    CYSb[1] = int(ActiveIDs[m].strip().split(" ")[0])
                    CYSc[1] = int(ActiveIDs[m].strip().split(" ")[1])
                    Ligp[1] = int(ActiveIDs[m].strip().split(" ")[2])
                    Ligd[1] = int(ActiveIDs[m].strip().split(" ")[3])
                elif ( HjLen == 8 ) and ( HjType == "c") and ( HjAxLigType == "HM" ):
                    SelHemJType = "cHM"
                    CYSb[1] = int(ActiveIDs[m].strip().split(" ")[0])
                    CYSc[1] = int(ActiveIDs[m].strip().split(" ")[1])
                    Ligp[1] = int(ActiveIDs[m].strip().split(" ")[2])
                    Ligd[1] = int(ActiveIDs[m].strip().split(" ")[3])
                elif ( HjLen == 6 ) and ( HjType == "b") and ( HjAxLigType == "HH" ):
                    SelHemJType = "bHH"
                    Ligp[1] = int(ActiveIDs[m].strip().split(" ")[0])
                    Ligd[1] = int(ActiveIDs[m].strip().split(" ")[1])
                elif ( HjLen == 6 ) and ( HjType == "b") and ( HjAxLigType == "HM" ):
                    SelHemJType = "bHM"
                    Ligp[1] = int(ActiveIDs[m].strip().split(" ")[0])
                    Ligd[1] = int(ActiveIDs[m].strip().split(" ")[1])
                else:
                    print(f" *** Missing entries on line number {m} of SelResIndexing.txt!")

#------------------------------------------------------------------------------
# For the selected heme pair, cycle through all the possible redox microstates
# while all other hemes are left in the reference redox state, which was 
# specified in DefineRefState.py
#------------------------------------------------------------------------------
        for k in ("o", "r"):
            for l in ("o", "r"):

#------------------------------------------------------------------------------
# Each heme in the pair can be oxidized or reduced. 
# These are either possibilities for each heme, and the residue names need to 
# be changed accordingly.
#   - oxidized c-type His-His ligated
#   - reduced  c-type His-His ligated
#   - oxidized c-type His-Met ligated
#   - reduced  c-type His-Met ligated
#   - oxidized b-type His-His ligated
#   - reduced  b-type His-His ligated
#   - oxidized b-type His-Met ligated
#   - reduced  b-type His-Met ligated
#------------------------------------------------------------------------------
                if (k == "o") and (HiType == "c") and (HiAxLigType == "HH"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Oxidized Heme-{Hi}: His-His ligated c-type heme 

 #Define atom groups
   set HISp  [atomselect top "resname PHO PHR and resid {Ligp[0]}"]
   set HISd  [atomselect top "resname DHO DHR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname HCO HCR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PHO; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHO; #Distal   His for oxidized His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCO; #c-type His-His ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (k == "r") and (HiType == "c") and (HiAxLigType == "HH"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Reduced Heme-{Hi}: His-His ligated c-type heme  

 #Define atom groups
   set HISp  [atomselect top "resname PHO PHR and resid {Ligp[0]}"]
   set HISd  [atomselect top "resname DHO DHR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname HCO HCR and resid {Hi}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp  set resname PHR; #Proximal His for reduced  His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHR; #Distal   His for reduced  His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCR; #c-type His-His ligated reduced  heme""", file=open(VMDinput, 'a'))
      
                if (k == "o") and (HiType == "c") and (HiAxLigType == "HM"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Oxidized Heme-{Hi}: His-Met ligated c-type heme  

 #Define atom groups
   set HISp  [atomselect top "resname PMO PMR and resid {Ligp[0]}"]
   set METd  [atomselect top "resname DMO DMR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname MCO MCR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PMO; #Proximal His for oxidized His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMO; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCO; #c-type His-Met ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (k == "r") and (HiType == "c") and (HiAxLigType == "HM"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Reduced Heme-{Hi}: His-Met ligated c-type heme   

 #Define atom groups
   set HISp  [atomselect top "resname PMO PMR and resid {Ligp[0]}"]
   set METd  [atomselect top "resname DMO DMR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname MCO MCR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PMR; #Proximal His for oxidized His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMR; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCR; #c-type His-Met ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (k == "o") and (HiType == "b") and (HiAxLigType == "HH"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Oxidized Heme-{Hi}: His-His ligated b-type heme   
 
 #Define atom groups
   set HISp  [atomselect top "resname FHO FHR and resid {Ligp[0]}"]
   set HISd  [atomselect top "resname RHO RHR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname HBO HBR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FHO; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHO; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBO; #b-type His-His ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (k == "r") and (HiType == "b") and (HiAxLigType == "HH"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Reduced Heme-{Hi}: His-His ligated b-type heme    
 
 #Define atom groups
   set HISp  [atomselect top "resname FHO FHR and resid {Ligp[0]}"]
   set HISd  [atomselect top "resname RHO RHR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname HBO HBR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FHR; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHR; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBR; #b-type His-His ligated oxidized heme""", file=open(VMDinput, 'a'))
    
                if (k == "o") and (HiType == "b") and (HiAxLigType == "HM"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Oxidized Heme-{Hi}: His-Met ligated b-type heme    

 #Define atom groups
   set HISp  [atomselect top "resname FMO FMR and resid {Ligp[0]}"]
   set METd  [atomselect top "resname RMO RMR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname MBO MBR and resid {Hi}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FMO; #Proximal His for oxidized His-Met ligated b-type heme.
   #The distal His residue 
     $METd set resname RMO; #Distal   Met for oxidized His-Met ligated b-type heme.
   #The heme group
     $HEM  set resname MBO; #b-type His-Met ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (k == "r") and (HiType == "b") and (HiAxLigType == "HM"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Reduced Heme-{Hi}: His-Met ligated b-type heme    

 #Define atom groups
   set HISp  [atomselect top "resname FMO FMR and resid {Ligp[0]}"]
   set METd  [atomselect top "resname RMO RMR and resid {Ligd[0]}"]
   set HEM   [atomselect top "resname MBO MBR and resid {Hi}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp set resname FMR; #Proximal His for reduced His-Met ligated b-type heme.
   #The distal His residue 
     $METd set resname RMR; #Distal   Met for reduced His-Met ligated b-type heme.
   #The heme group
     $HEM  set resname MBR; #b-type His-Met ligated reduced heme""", file=open(VMDinput, 'a'))

                if (l == "o") and (HjType == "c") and (HjAxLigType == "HH"):
                    print(f"""
 #Oxidized Heme-{Hj}: His-His ligated c-type heme    

 #Define atom groups
   set HISp  [atomselect top "resname PHO PHR and resid {Ligp[1]}"]
   set HISd  [atomselect top "resname DHO DHR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname HCO HCR and resid {Hj}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PHO; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHO; #Distal   His for oxidized His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCO; #c-type His-His ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (l == "r") and (HjType == "c") and (HjAxLigType == "HH"):
                    print(f"""
 #Reduced Heme-{Hj}: His-His ligated c-type heme     

 #Define atom groups
   set HISp  [atomselect top "resname PHO PHR and resid {Ligp[1]}"]
   set HISd  [atomselect top "resname DHO DHR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname HCO HCR and resid {Hj}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp  set resname PHR; #Proximal His for reduced His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHR; #Distal   His for reduced His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCR; #c-type His-His ligated reduced heme""", file=open(VMDinput, 'a'))
      
                if (l == "o") and (HjType == "c") and (HjAxLigType == "HM"):
                    print(f"""
 #Oxidized Heme-{Hj}: His-Met ligated c-type heme     

 #Define atom groups
   set HISp  [atomselect top "resname PMO PMR and resid {Ligp[1]}"]
   set METd  [atomselect top "resname DMO DMR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname MCO MCR and resid {Hj}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PMO; #Proximal His for oxidized His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMO; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCO; #c-type His-Met ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (l == "r") and (HjType == "c") and (HjAxLigType == "HM"):
                    print(f"""
 #Reduced Heme-{Hj}: His-Met ligated c-type heme     

 #Define atom groups
   set HISp  [atomselect top "resname PMO PMR and resid {Ligp[1]}"]
   set METd  [atomselect top "resname DMO DMR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname MCO MCR and resid {Hj}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp  set resname PMR; #Proximal His for reduced His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMR; #Distal   Met for reduced His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCR; #c-type His-Met ligated reduced heme""", file=open(VMDinput, 'a'))

                if (l == "o") and (HjType == "b") and (HjAxLigType == "HH"):
                    print(f"""
 #Oxidized Heme-{Hj}: His-His ligated b-type heme      
 
 #Define atom groups
   set HISp  [atomselect top "resname FHO FHR and resid {Ligp[1]}"]
   set HISd  [atomselect top "resname RHO RHR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname HBO HBR and resid {Hj}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FHO; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHO; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBO; #b-type His-His ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (l == "r") and (HjType == "b") and (HjAxLigType == "HH"):
                    print(f"""
 #Reduced Heme-{Hj}: His-His ligated b-type heme      
 
 #Define atom groups
   set HISp  [atomselect top "resname FHO FHR and resid {Ligp[1]}"]
   set HISd  [atomselect top "resname RHO RHR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname HBO HBR and resid {Hj}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp set resname FHR; #Proximal His for reduced His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHR; #Distal   HIS for reduced His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBR; #b-type His-His ligated reduced heme""", file=open(VMDinput, 'a'))
    
                if (l == "o") and (HjType == "b") and (HjAxLigType == "HM"):
                    print(f"""
 #Oxidized Heme-{Hj}: His-Met ligated b-type heme      

 #Define atom groups
   set HISp  [atomselect top "resname FMO FMR and resid {Ligp[1]}"]
   set METd  [atomselect top "resname RMO RMR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname MBO MBR and resid {Hj}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FMO; #Proximal His for oxidized His-Met ligated b-type heme.
   #The distal His residue 
     $METd set resname RMO; #Distal   Met for oxidized His-Met ligated b-type heme.
   #The heme group
     $HEM  set resname MBO; #b-type His-Met ligated oxidized heme""", file=open(VMDinput, 'a'))

                if (l == "r") and (HjType == "b") and (HjAxLigType == "HM"):
                    print(f"""
 #-------------------------------------------------------------------------
 #Reduced Heme-{Hj}: His-Met ligated b-type heme       

 #Define atom groups
   set HISp  [atomselect top "resname FMO FMR and resid {Ligp[1]}"]
   set METd  [atomselect top "resname RMO RMR and resid {Ligd[1]}"]
   set HEM   [atomselect top "resname MBO MBR and resid {Hj}"]

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp set resname FMR; #Proximal His for reduced His-Met ligated b-type heme.
   #The distal His residue 
     $METd set resname RMR; #Distal   Met for reduced His-Met ligated b-type heme.
   #The heme group
     $HEM  set resname MBR; #b-type His-Met ligated reduced heme""", file=open(VMDinput, 'a'))

                print(f"""
 set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
 $sel writepdb temp_{k}{Hi}-{l}{Hj}.pdb
 #-------------------------------------------------------------------------""", file=open(VMDinput, 'a'))

        print(f""" \n exit """, file=open(VMDinput, 'a'))

        print(f"""
 Using VMD to generate redox microstate PDBs for Heme-{Hi} and Heme-{Hj}... """)
        subprocess.run(f"vmd -e PairInt_{Hi}-{Hj}.tcl > PairInt_{Hi}-{Hj}.log", shell=True)

        chk=0
        for k in ("o", "r"):
            for l in ("o", "r"):
                if (os.path.isfile(f"temp_{k}{Hi}-{l}{Hj}.pdb") == True):
                    print(f"  VMD successfully generated temp_{k}{Hi}-{l}{Hj}.pdb")
                    insert_ter_records(f"temp_{k}{Hi}-{l}{Hj}.pdb", f"{k}{Hi}-{l}{Hj}.pdb")
                    print(f"  Added TER records to generated PDB: {k}{Hi}-{l}{Hj}.pdb")
                    if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.pdb") == True):
                        os.remove(f"temp_{k}{Hi}-{l}{Hj}.pdb")
                        print(f"  Deleted temporary PDB: temp_{k}{Hi}-{l}{Hj}.pdb")
                        print(f"""  √ {k}{Hi} -- {l}{Hj}: VMD finished successfully!""")
                if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.pdb") == False):
                    chk+=1
                    print(f"""
   X temp_{k}{Hi} -- {l}{Hj}: VMD failed. 
   Please check PairInt_{Hi}-{Hj}.log.""")
    
        if (chk != 0):
            sys.exit(f"""
 VMD failed to generate the PDB for one or more redox microstates. 
 Please inspect PairInt_{Hi}-{Hj}.log before re-running this module \n""")

#------------------------------------------------------------------------------
# Now create the rest of the TLEaP input file, making sure only to specify 
# atom types and forcefield files if one or more hemes in the structure is 
# of that type, but not to specify multiple types of those entires.
#------------------------------------------------------------------------------
        if (HiType == "c" and HiAxLigType == "HH") or (HjType == "c" and HiAxLigType == "HH"):
            print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }
        { "M8"  "Fe" "sp3" } #M8&T1-T6:
        { "T1"  "N" "sp3" }  #Reduced
        { "T2"  "N" "sp3" }  #His-His
        { "T3"  "N" "sp3" }  #Ligated
        { "T4"  "N" "sp3" }  #c-Heme
        { "T5"  "N" "sp3" }
        { "T6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (HiType == "c" and HiAxLigType == "HM") or (HjType == "c" and HiAxLigType == "HM"):
            print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (HiType == "b" and HiAxLigType == "HH") or (HjType == "b" and HiAxLigType == "HH"):
            print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (HiType == "b" and HiAxLigType == "HM") or (HjType == "b" and HjAxLigType == "HM"):
            print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "S" "sp3" }  #Oxidized
        { "W2"  "N" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "S" "sp3" }  #Reduced
        { "X2"  "N" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_b_HH != 0 and SelHemIType != "bHH") or (SelRefRedoxState == "O" and Count_b_HH != 0 and SelHemJType != "bHH"):
            print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_b_HH != 0 and SelHemIType != "bHH") or (SelRefRedoxState == "R" and Count_b_HH != 0 and SelHemJType != "bHH"):
            print("""
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_b_HM != 0 and SelHemIType != "bHM") or (SelRefRedoxState == "O" and Count_b_HM != 0 and SelHemJType != "bHM"):
            print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "S" "sp3" }  #Oxidized
        { "W2"  "N" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_b_HM != 0 and SelHemIType != "bHM") or (SelRefRedoxState == "R" and Count_b_HM != 0 and SelHemJType != "bHM"):
            print("""
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "S" "sp3" }  #Reduced
        { "X2"  "N" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_c_HM != 0 and SelHemIType != "cHM") or (SelRefRedoxState == "O" and Count_c_HM != 0 and SelHemJType != "cHM"):
            print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_c_HM != 0 and SelHemIType != "cHM") or (SelRefRedoxState == "R" and Count_c_HM != 0 and SelHemJType != "cHM"):
            print("""
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_c_HH != 0 and SelHemIType != "cHH") or (SelRefRedoxState == "O" and Count_c_HH != 0 and SelHemJType != "cHH"):
            print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_c_HH != 0 and SelHemIType != "cHH") or (SelRefRedoxState == "R" and Count_c_HH != 0 and SelHemJType != "cHH"):
            print("""
        { "M8"  "Fe" "sp3" } #M8&T1-T6:
        { "T1"  "N" "sp3" }  #Reduced
        { "T2"  "N" "sp3" }  #His-His
        { "T3"  "N" "sp3" }  #Ligated
        { "T4"  "N" "sp3" }  #c-Heme
        { "T5"  "N" "sp3" }
        { "T6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

        print("""\n } """, file=open(TLEaPinput, 'a'))

        if ( Count_b_HH != 0 ) or ( Count_b_HM != 0 ):
            print("""
# References for b-type heme forcefield parameters:
#    Bonded parameters for the macrocycle come from:
#      Yang, Longhua, Åge A. Skjevik, Wen-Ge Han Du, Louis Noodleman, Ross C. Walker, and Andreas W. Götz.
#      Data for molecular dynamics simulations of B-type cytochrome c oxidase with the Amber force field.
#      Data in brief 8 (2016): 1209-1214.
#
#    Bonded parameters for the Fe center and atomic partial charges were derived by Guberman-Pfeffer 
#    using the Metal Center Parameter Builder. The B3LYP approximate density functional was used with 
#    the z mixed basis set (LANL2TZ(f) for Fe and 6-31G(d) for 2nd row elements. 
#
#    A different set of charges is available in the literature (below reference), but only for the 
#    oxidized redox state. Also, in the developmenet of BioDC, Guberman-Pfeffer liked the idea of
#    having a consistently-derived set of parameters for b- and c-type hemes with His-His and 
#    His-Met ligation.
#
#    Alternative set of charges are available at:
#      L.Noodleman et al. Inorg. Chem., 53 (2014)
#      6458;
#      J.A.Fee et al. J.Am.Chem.Soc., 130 (2008) 15002. 
""", end=" ", file=open(TLEaPinput, 'a'))

        if (HiType == "b" and HiAxLigType == "HH") or (HjType == "b" and HiAxLigType == "HH"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", file=open(TLEaPinput, 'a'))
        if (HiType == "b" and HiAxLigType == "HM") or (HjType == "b" and HiAxLigType == "HM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib""", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_b_HH != 0 and SelHemIType != "bHH") or (SelRefRedoxState == "O" and Count_b_HH != 0 and SelHemJType != "bHH"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_b_HH != 0 and SelHemIType != "bHH") or (SelRefRedoxState == "R" and Count_b_HH != 0 and SelHemJType != "bHH"):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "O" and Count_b_HM != 0 and SelHemIType != "bHM") or (SelRefRedoxState == "O" and Count_b_HM != 0 and SelHemJType != "bHM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_b_HM != 0 and SelHemIType != "bHM") or (SelRefRedoxState == "R" and Count_b_HM != 0 and SelHemJType != "bHM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

        if ( Count_c_HH != 0 ) or ( Count_c_HM != 0 ):
            print("""
# References for c-type heme forcefield parameters:
#    Bonded parameters for the macrocycle come from:
#      Crespo, A.; Martí, M. A.; Kalko, S. G.; Morreale, A.; Orozco, M.; Gelpi, J. L.; Luque, F. J.; 
#      Estrin, D. A. Theoretical Study of the Truncated Hemoglobin HbN: Exploring the Molecular Basis 
#      of the NO Detoxification Mechanism. J. Am. Chem. Soc. 2005, 127 (12), 4433–4444.
#
#    Bonded parameters for the Fe center and atomic partial charges were derived by Guberman-Pfeffer 
#    using the Metal Center Parameter Builder. The B3LYP approximate density functional was used with 
#    the z mixed basis set (LANL2TZ(f) for Fe and 6-31G(d) for 2nd row elements. 
#
#    A different set of charges is available in the literature (below reference), but in the 
#    developmenet of BioDC, Guberman-Pfeffer liked the idea of having a consistently-derived 
#    set of parameters for b- and c-type hemes with His-His and His-Met ligation.
#
#    Alternative set of charges are available at:
#      Henriques, J.; Costa, P. J.; Calhorda, M. J.; Machuqueiro, M. Charge Parametrization 
#      of the DvH-c3 Heme Group: Validation Using Constant-(pH,E) Molecular Dynamics 
#      Simulations. J. Phys. Chem. B 2013, 117 (1), 70–82.
""", end=" ", file=open(TLEaPinput, 'a'))

        if (HiType == "c" and HiAxLigType == "HH") or (HjType == "c" and HiAxLigType == "HH"):
            if (FFchoice.lower() in ["henriques", "h"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", file=open(TLEaPinput, 'a'))
            elif (FFchoice.lower() in ["guberman-pfeffer", "gp"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", file=open(TLEaPinput, 'a'))
        if (HiType == "c" and HiAxLigType == "HM") or (HjType == "c" and HiAxLigType == "HM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_c_HM != 0 and SelHemIType != "cHM") or (SelRefRedoxState == "O" and Count_c_HM != 0 and SelHemJType != "cHM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_c_HM != 0 and SelHemIType != "cHM") or (SelRefRedoxState == "R" and Count_c_HM != 0 and SelHemJType != "cHM"):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

        if (SelRefRedoxState == "O" and Count_c_HH != 0 and SelHemIType != "cHH") or (SelRefRedoxState == "O" and Count_c_HH != 0 and SelHemJType != "cHH"):
            if (FFchoice.lower() in ["henriques", "h"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            elif (FFchoice.lower() in ["guberman-pfeffer", "gp"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
        if (SelRefRedoxState == "R" and Count_c_HH != 0 and SelHemIType != "cHH") or (SelRefRedoxState == "R" and Count_c_HH != 0 and SelHemJType != "cHH"):
            if (FFchoice.lower() in ["henriques", "h"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            elif (FFchoice.lower() in ["guberman-pfeffer", "gp"]):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

        print(f"""
# Load PDB
  oo = loadpdb o{Hi}-o{Hj}.pdb
  or = loadpdb o{Hi}-r{Hj}.pdb
  ro = loadpdb r{Hi}-o{Hj}.pdb
  rr = loadpdb r{Hi}-r{Hj}.pdb""", file=open(TLEaPinput, 'a'))

#------------------------------------------------------------------------------
# Specify the disulfide linkages. The residue IDs for each pair are read from 
# DisulfideDefinitions.txt, which was created by SelectDisulfides.py 
#------------------------------------------------------------------------------
        if (os.path.isfile("DisulfideDefinitions.txt") == True):
            with open("DisulfideDefinitions.txt") as dsl:
                NumDisulfide = int(len(dsl.readlines()))
                DisulfPairID = [0]*NumDisulfide
                dsl.seek(0)

                idx=0
                Lines_dsl = dsl.readlines()
                for line in Lines_dsl:
                    SelPairIDs = line
                    DisulfPairID[idx] = list(map(int,SelPairIDs.split()))
                    idx+=1

            if (len(DisulfPairID) != 0):
                print("# Define Disulfide linkages ")
                for sbi in range(len(DisulfPairID)):
                    print(f" bond  oo.{DisulfPairID[sbi][0]}.SG oo.{DisulfPairID[sbi][1]}.SG", file=open(TLEaPinput, 'a'))
                    print(f" bond  or.{DisulfPairID[sbi][0]}.SG or.{DisulfPairID[sbi][1]}.SG", file=open(TLEaPinput, 'a'))
                    print(f" bond  ro.{DisulfPairID[sbi][0]}.SG ro.{DisulfPairID[sbi][1]}.SG", file=open(TLEaPinput, 'a'))
                    print(f" bond  rr.{DisulfPairID[sbi][0]}.SG rr.{DisulfPairID[sbi][1]}.SG", file=open(TLEaPinput, 'a'))

#------------------------------------------------------------------------------
# Regardless of which hemes constitute the selected pair, bond definitions 
# need to be provided for each heme in the structure. Thus, we read in 
# ResIndexing.txt now instead of SelResIndexing.txt. The bond definitions are
# quadrupled because each heme is present in the structure for the four 
# microstates of the selected heme pair, namely, OO, OR, RO, and RR.
#------------------------------------------------------------------------------

        # Add this set to keep track of defined bonds
        defined_bonds = set()

        idx=0
        with open("ResIndexing.txt") as fp:
            Lines = fp.readlines()
            for line in Lines:
                EntryLength = len(line.strip().split(" "))
                HemeType = line.strip().split(" ")[-2]
                AxLigType = line.strip().split(" ")[-1]

                if ( EntryLength == 8 ) and ( HemeType == "c"):
                    idx+=1
                    CYSb = int(line.strip().split(" ")[0])
                    CYSc = int(line.strip().split(" ")[1])
                    Ligp = int(line.strip().split(" ")[2])
                    Ligd = int(line.strip().split(" ")[3])
                    HEM = int(line.strip().split(" ")[5])
                if ( EntryLength == 6 ) and ( HemeType == "b"):
                    idx+=1
                    Ligp = int(line.strip().split(" ")[0])
                    Ligd = int(line.strip().split(" ")[1])
                    HEM = int(line.strip().split(" ")[3])

                print(f"""
#------------------------------------------------------------
#For hemes {HEM}

#Bond ligating atoms to Fe center""", file=open(TLEaPinput, 'a'))

                print(f""" bond  oo.{Ligp}.NE2   oo.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                if (AxLigType == "HH"):
                    print(f""" bond  oo.{Ligd}.NE2   oo.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                elif (AxLigType == "HM"):
                    print(f""" bond  oo.{Ligd}.SD   oo.{HEM}.FE""", file=open(TLEaPinput, 'a'))

                print(f"""\n bond  or.{Ligp}.NE2   or.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                if (AxLigType == "HH"):
                    print(f""" bond  or.{Ligd}.NE2   or.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                elif (AxLigType == "HM"):
                    print(f""" bond  or.{Ligd}.SD   or.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                    
                print(f"""\n bond  ro.{Ligp}.NE2   ro.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                if (AxLigType == "HH"):
                    print(f""" bond  ro.{Ligd}.NE2   ro.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                elif (AxLigType == "HM"):
                    print(f""" bond  ro.{Ligd}.SD   ro.{HEM}.FE""", file=open(TLEaPinput, 'a'))

                print(f"""\n bond  rr.{Ligp}.NE2   rr.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                if (AxLigType == "HH"):
                    print(f""" bond  rr.{Ligd}.NE2   rr.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                elif (AxLigType == "HM"):
                    print(f""" bond  rr.{Ligd}.SD   rr.{HEM}.FE""", file=open(TLEaPinput, 'a'))
  
                print(f"""
#Bond axially coordinated residues to preceeding and proceeding residues """, file=open(TLEaPinput, 'a'))

                # Check and write bonds only if they haven't been defined yet
                for bond in [
                    (Ligp-1, Ligp),
                    (Ligp, Ligp+1),
                    (Ligd-1, Ligd),
                    (Ligd, Ligd+1)
                ]:
                    if bond not in defined_bonds:
                        print(f" bond oo.{bond[0]}.C   oo.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                        print(f" bond or.{bond[0]}.C   or.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                        print(f" bond ro.{bond[0]}.C   ro.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                        print(f" bond rr.{bond[0]}.C   rr.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                        defined_bonds.add(bond)

                if (HemeType == "c"):
                    print(f"""
#Bond heme thioethers to protein backbone""", end=" ", file=open(TLEaPinput, 'a'))
                if (HemeType == "c"):
                    print(f"""
 bond  oo.{CYSb}.CA   oo.{HEM}.CBB2
 bond  oo.{CYSc}.CA   oo.{HEM}.CBC1

 bond  or.{CYSb}.CA   or.{HEM}.CBB2
 bond  or.{CYSc}.CA   or.{HEM}.CBC1

 bond  ro.{CYSb}.CA   ro.{HEM}.CBB2
 bond  ro.{CYSc}.CA   ro.{HEM}.CBC1

 bond  rr.{CYSb}.CA   rr.{HEM}.CBB2
 bond  rr.{CYSc}.CA   rr.{HEM}.CBC1
 """, end=" ", file=open(TLEaPinput, 'a'))

                print(f"""
#Bond propionic acids to heme
 bond  oo.{HEM}.C2A   oo.{HEM+1}.CA
 bond  oo.{HEM}.C3D   oo.{HEM+2}.CA

 bond  or.{HEM}.C2A   or.{HEM+1}.CA
 bond  or.{HEM}.C3D   or.{HEM+2}.CA

 bond  ro.{HEM}.C2A   ro.{HEM+1}.CA
 bond  ro.{HEM}.C3D   ro.{HEM+2}.CA

 bond  rr.{HEM}.C2A   rr.{HEM+1}.CA
 bond  rr.{HEM}.C3D   rr.{HEM+2}.CA
 """, file=open(TLEaPinput, 'a'))

        print(f"""
# Save topology and coordinate files
 saveamberparm  oo o{Hi}-o{Hj}.prmtop o{Hi}-o{Hj}.rst7
 saveamberparm  or o{Hi}-r{Hj}.prmtop o{Hi}-r{Hj}.rst7
 saveamberparm  ro r{Hi}-o{Hj}.prmtop r{Hi}-o{Hj}.rst7
 saveamberparm  rr r{Hi}-r{Hj}.prmtop r{Hi}-r{Hj}.rst7

quit""", file=open(TLEaPinput, 'a'))

        print(f"""
 Using TLEaP to build the redox microstate topologies for Heme-{Hi} and Heme-{Hj}... """)
        subprocess.run(f"tleap -s -f GeneratePairIntTopologiesForHems{Hi}-{Hj}.in > GeneratePairIntTopologiesForHems{Hi}-{Hj}.log", shell=True)

        chk=0
        for k in ("o", "r"):
            for l in ("o", "r"):
                if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == True) and (os.path.isfile(f"{k}{Hi}-{l}{Hj}.rst7") == True):
                    print(f"""  √ {k}{Hi} -- {l}{Hj}: TLEaP finished successfully!""")
                if (os.path.isfile(f"{k}{Hi}-{l}{Hj}.prmtop") == False) or (os.path.isfile(f"{k}{Hi}-{l}{Hj}.rst7") == False):
                    chk+=1
                    print(f"""
   X {k}{Hi} -- {l}{Hj}: TLEaP failed. 
   Please check GeneratePairIntTopologiesForHems{Hi}-{Hj}.log""")
    
        if (chk != 0):
            sys.exit(f"""
 TLEaP failed to build the topologies for one or more redox microstates. 
 Please inspect PairInt_{Hi}-{Hj}.log and GeneratePairIntTopologiesForHems{Hi}-{Hj}.log 
 before re-running this module \n""")
            
