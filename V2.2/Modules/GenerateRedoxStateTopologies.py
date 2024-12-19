################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

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

def GenerateRedoxStateTopologies(ForceFieldDir, FFchoice, SelRefRedoxState):

    with open("SelResIndexing.txt") as sri:
        Len_sri = len(sri.readlines())
        SelHEM = [0]*Len_sri
        sri.seek(0)

        idxs = 0
        Lines_sri = sri.readlines()
        for line in Lines_sri:
            SelHEM[idxs] = int(line.strip().split(" ")[-3])

            EntryLength = len(line.strip().split(" "))
            HemeType = line.strip().split(" ")[-2]
            AxLigType = line.strip().split(" ")[-1]

            if ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HH" ):
                SelHemIDType = "cHH"
                CYSb = int(line.strip().split(" ")[0])
                CYSc = int(line.strip().split(" ")[1])
                HISp = int(line.strip().split(" ")[2])
                HISd = int(line.strip().split(" ")[3])
                HEM  = int(line.strip().split(" ")[5])

                print(f""" 
 Reading Info. for Heme-{HEM} as a bis-His-ligated c-type heme.
   Residue ID of Thioether-linked Cys to heme b-ring: {CYSb}
   Residue ID of Thioether-linked Cys to heme c-ring: {CYSc}
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal His to   Fe-center: {HISd}\n""")

                print(f"""
 #-------------------------------------------------------------------------
 #Input Structure
  mol new RefState.pdb

 #-------------------------------------------------------------------------
 #Oxidized Heme-{HEM}: 

 #Define atom groups
   set HISp  [atomselect top "resname PHO PHR and resid {HISp}"]
   set HISd  [atomselect top "resname DHO DHR and resid {HISd}"]
   set HEM   [atomselect top "resname HCO HCR and resid {HEM}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PHO; #Proximal His for oxidized His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHO; #Distal   His for oxidized His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCO; #c-type His-His ligated oxidized heme

 #Write PDB for oxidized state 
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_o{HEM}.pdb

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp  set resname PHR; #Proximal His for reduced  His-His ligated heme.
   #The distal His residue
     $HISd  set resname DHR; #Distal   His for reduced  His-His ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname HCR; #c-type His-His ligated reduced heme

 #Write PDB for reduced state 
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_r{HEM}.pdb
 #-------------------------------------------------------------------------

 exit""", file=open(f"SetRedoxStatesForHem{HEM}.tcl", 'w')) 

            elif ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HM" ):
                SelHemIDType = "cHM"
                CYSb = int(line.strip().split(" ")[0])
                CYSc = int(line.strip().split(" ")[1])
                HISp = int(line.strip().split(" ")[2])
                METd = int(line.strip().split(" ")[3])
                HEM = int(line.strip().split(" ")[5])

                print(f""" 
 Reading Info. for Heme-{HEM} as a His-Met-ligated c-type heme.
   Residue ID of Thioether-linked Cys to heme b-ring: {CYSb}
   Residue ID of Thioether-linked Cys to heme c-ring: {CYSc}
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal Met to   Fe-center: {METd} \n""")

                print(f"""
 #-------------------------------------------------------------------------
 #Input Structure
  mol new RefState.pdb

 #-------------------------------------------------------------------------
 #Oxidized Heme-{HEM}: 

 #Define atom groups
   set HISp  [atomselect top "resname PMO PMR and resid {HISp}"]
   set METd  [atomselect top "resname DMO DMR and resid {METd}"]
   set HEM   [atomselect top "resname MCO MCR and resid {HEM}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp  set resname PMO; #Proximal His for oxidized His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMO; #Distal   Met for oxidized His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCO; #c-type His-Met ligated oxidized heme

 #Write PDB for oxidized state 
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_o{HEM}.pdb

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp  set resname PMR; #Proximal His for reduced  His-Met ligated heme.
   #The distal His residue
     $METd  set resname DMR; #Distal   Met for reduced  His-Met ligated heme.
   #The heme group including the thioether linkages from Cys sidechains
     $HEM   set resname MCR; #c-type His-Met ligated redued heme

 #Write PDB for reduced state: 
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_r{HEM}.pdb
 #-------------------------------------------------------------------------

 exit""", file=open(f"SetRedoxStatesForHem{HEM}.tcl", 'w')) 

            elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HH" ):
                SelHemIDType = "bHH"
                HISp = int(line.strip().split(" ")[0])
                HISd = int(line.strip().split(" ")[1])
                HEM = int(line.strip().split(" ")[3])

                print(f""" 
 Reading Info. for Heme-{HEM} as a bis-His-ligated b-type heme.
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal His to   Fe-center: {HISd} \n""")

                print(f"""
 #-------------------------------------------------------------------------
 #Input Structure
  mol new RefState.pdb

 #-------------------------------------------------------------------------
 #Oxidized Heme-{HEM}: 
 
 #Define atom groups
   set HISp  [atomselect top "resname FHO FHR and resid {HISp}"]
   set HISd  [atomselect top "resname RHO RHR and resid {HISd}"]
   set HEM   [atomselect top "resname HBO HBR and resid {HEM}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FHO; #Proximal His for oxidized His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHO; #Distal   HIS for oxidized His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBO; #b-type His-His ligated oxidized heme

 #Write PDB for oxidizeds
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_o{HEM}.pdb

 #Change residue names for:
   #The proximal His residue 
     $HISp set resname FHR; #Proximal His for reduced  His-His ligated b-type heme.
   #The distal His residue 
     $HISd set resname RHR; #Distal   HIS for reduced  His-His ligated b-type heme.
   #The heme group
     $HEM  set resname HBR; #b-type His-His ligated reduced heme

 #Write PDBs for reduced state:
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_r{HEM}.pdb
 #-------------------------------------------------------------------------

 exit""", file=open(f"SetRedoxStatesForHem{HEM}.tcl", 'w')) 

            elif (EntryLength == 6) and (HemeType == "b") and (AxLigType == "HM"):
                SelHemIDType = "bHM"
                HISp = int(line.strip().split(" ")[0])
                METd = int(line.strip().split(" ")[1])
                HEM = int(line.strip().split(" ")[3])
                
                print(f""" 
 Reading Info. for Heme-{HEM} as a His-Met-ligated b-type heme.
   Residue ID of         proximal His to   Fe-center: {HISp}
   Residue ID of           distal Met to   Fe-center: {METd}\n""")

                print(f"""
 #-------------------------------------------------------------------------
 #Input Structure
  mol new RefState.pdb

 #-------------------------------------------------------------------------
 #Oxidized Heme-{HEM}: 

 #Define atom groups
   set HISp  [atomselect top "resname FMO FMR and resid {HISp}"]
   set METd  [atomselect top "resname RMO RMR and resid {METd}"]
   set HEM   [atomselect top "resname MBO MBR and resid {HEM}"]

 #Change residue names for oxidized state:
   #The proximal His residue 
     $HISp set resname FMO; #Proximal His for oxidized His-Met ligated b-type heme.
   #The distal His residue 
     $METd set resname RMO; #Distal   Met for oxidized His-Met ligated b-type heme.
   #The heme group
     $HEM  set resname MBO; #b-type His-Met ligated oxidized heme

 #Write PDBs
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_o{HEM}.pdb

 #Change residue names for reduced state:
   #The proximal His residue 
     $HISp set resname FMR; #Proximal His for reduced  His-His ligated b-type heme.
   #The distal His residue 
     $METd set resname RMR; #Distal   Met for reduced  His-His ligated b-type heme.
   #The heme group
     $HEM  set resname MBR; #b-type His-Met ligated reduced heme

 #Write PDBs for reduced state
   set sel [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
   $sel writepdb temp_r{HEM}.pdb
 #-------------------------------------------------------------------------

 exit""", file=open(f"SetRedoxStatesForHem{HEM}.tcl", 'w')) 

            print(f"   Running VMD to generate o{HEM}.pdb and r{HEM}.pdb ...")
            subprocess.run(f"vmd -dispdev text -e SetRedoxStatesForHem{HEM}.tcl > SetRedoxStatesForHem{HEM}.log", shell=True)

            if (os.path.isfile(f"temp_o{HEM}.pdb") == True):
                print(f"  VMD successfully generated o{HEM}.PDB")
                insert_ter_records(f"temp_o{HEM}.pdb", f"o{HEM}.pdb")
                print(f"   Added TER records to generated PDB: o{HEM}.pdb")
                if (os.path.isfile(f"o{HEM}.pdb") == True):
                    os.remove(f"temp_o{HEM}.pdb")
                    print(f"  Deleted temporary PDB: temp_o{HEM}.pdb")
            elif (os.path.isfile(f"temp_o{HEM}.pdb") == False):
                sys.exit(f"""     
 VMD failed to generated o{HEM}.pdb. 
 Please check SetRedoxStatesForHem{HEM}.log 
 to diagnose the problem. \n""")

            if (os.path.isfile(f"temp_r{HEM}.pdb") == True):
                print(f"  VMD successfully generated r{HEM}.PDB")
                insert_ter_records(f"temp_r{HEM}.pdb", f"r{HEM}.pdb")
                print(f"   Added TER records to generated PDB: r{HEM}.pdb")
                if (os.path.isfile(f"o{HEM}.pdb") == True):
                    os.remove(f"temp_r{HEM}.pdb")
                    print(f"  Deleted temporary PDB: temp_r{HEM}.pdb")
            elif (os.path.isfile(f"temp_r{HEM}.pdb") == False):
                sys.exit(f"""     
 VMD failed to generated temp_r{HEM}.pdb. 
 Please check SetRedoxStatesForHem{HEM}.log 
 to diagnose the problem. \n""")

            Count_c_HH = 0  
            Count_c_HM = 0  
            Count_b_HH = 0  
            Count_b_HM = 0  
            TLEaPinput=f"GenerateRedoxStateTopologiesForHem{SelHEM[idxs]}.in"

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
                sys.exit("""
 The total number of hemes in ResIndexing.txt does NOT equl
 the sum of His-His and His-Met ligated b- and c-type hemes.
 Other types of hemes cannot yet be analyzed with BiodC.
 Please revise ResIndexing.txt before re-running BioDC.""")

            print("""
# Load parameters
 source leaprc.constph
 source leaprc.conste
 source leaprc.gaff
 source leaprc.water.tip3p

 addAtomTypes {""", end=' ', file=open(TLEaPinput, 'w'))
    
            if (SelHemIDType == "bHH"):
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

            if (SelHemIDType == "bHM"):
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

            if (SelHemIDType == "cHM"):
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

            if (SelHemIDType == "cHH"):
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

            if (SelRefRedoxState == "O") and (Count_b_HH != 0) and (SelHemIDType != "bHH"):
                print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelRefRedoxState == "R") and (Count_b_HH != 0) and (SelHemIDType != "bHH"):
                print("""
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

            if (SelRefRedoxState == "O") and (Count_b_HM != 0) and (SelHemIDType != "bHM"):
                print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "S" "sp3" }  #Oxidized
        { "W2"  "N" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelRefRedoxState == "R") and (Count_b_HM != 0) and (SelHemIDType != "bHM"):
                print("""
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "S" "sp3" }  #Reduced
        { "X2"  "N" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

            if (SelRefRedoxState == "O") and (Count_c_HM != 0) and (SelHemIDType != "cHM"):
                print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelRefRedoxState == "R") and (Count_c_HM != 0) and (SelHemIDType != "cHM"):
                print("""
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))

            if (SelRefRedoxState == "O") and (Count_c_HH != 0) and (SelHemIDType != "cHH"):
                print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelRefRedoxState == "R") and (Count_c_HH != 0) and (SelHemIDType != "cHH"):
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
#      L.Noodleman et al. Inorg. Chem., 53 (2014) 6458;
#      J.A.Fee et al. J.Am.Chem.Soc., 130 (2008) 15002. 
""", end=" ", file=open(TLEaPinput, 'a'))

            if (SelHemIDType == "bHH"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelHemIDType == "bHM"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

            if ( Count_b_HH != 0 ) and (SelRefRedoxState == "O") and (SelHemIDType != "bHH"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if ( Count_b_HH != 0 ) and (SelRefRedoxState == "R") and (SelHemIDType != "bHH"):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if ( Count_b_HM != 0 ) and (SelRefRedoxState == "O") and (SelHemIDType != "bHM"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if ( Count_b_HM != 0 ) and (SelRefRedoxState == "R") and (SelHemIDType != "bHM"):
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

            if (SelHemIDType == "cHM"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if (SelHemIDType == "cHH"):
                if (FFchoice == "HENRIQUES") or (FFchoice == "Henriques") or (FFchoice == "henriques") or (FFchoice == "H") or (FFchoice == "h"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
                elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

            if ( Count_c_HH != 0 ) and (SelRefRedoxState == "O") and (SelHemIDType != "cHH"):
                if (FFchoice == "HENRIQUES") or (FFchoice == "Henriques") or (FFchoice == "henriques") or (FFchoice == "H") or (FFchoice == "h"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
                elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if ( Count_c_HH != 0 ) and (SelRefRedoxState == "R") and (SelHemIDType != "cHH"):
                if (FFchoice == "HENRIQUES") or (FFchoice == "Henriques") or (FFchoice == "henriques") or (FFchoice == "H") or (FFchoice == "h"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
                elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

            if ( Count_c_HM != 0 ) and (SelRefRedoxState == "O") and (SelHemIDType != "cHM"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))
            if ( Count_c_HM != 0 ) and (SelRefRedoxState == "R") and (SelHemIDType != "cHM"):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open(TLEaPinput, 'a'))

            print(f"""

# Load PDB
 ox = loadpdb o{SelHEM[idxs]}.pdb
 rd = loadpdb r{SelHEM[idxs]}.pdb""", file=open(TLEaPinput, 'a'))
      
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
                        print(f" bond mol.{DisulfPairID[sbi][0]}.SG mol.{DisulfPairID[sbi][1]}.SG", file=open(TLEaPinput, 'a'))

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

                    print(f""" bond  ox.{Ligp}.NE2   ox.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                    if (AxLigType == "HH"):
                        print(f""" bond  ox.{Ligd}.NE2   ox.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                    elif (AxLigType == "HM"):
                        print(f""" bond  ox.{Ligd}.SD   ox.{HEM}.FE""", file=open(TLEaPinput, 'a'))

                    print(f"""\n bond  rd.{Ligp}.NE2   rd.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                    if (AxLigType == "HH"):
                        print(f""" bond  rd.{Ligd}.NE2   rd.{HEM}.FE""", file=open(TLEaPinput, 'a'))
                    elif (AxLigType == "HM"):
                        print(f""" bond  rd.{Ligd}.SD   rd.{HEM}.FE""", file=open(TLEaPinput, 'a'))
  
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
                            print(f" bond ox.{bond[0]}.C   ox.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                            print(f" bond rd.{bond[0]}.C   rd.{bond[1]}.N", file=open(TLEaPinput, 'a'))
                            defined_bonds.add(bond)

                    if (HemeType == "c"):
                        print(f"""
#Bond heme thioethers to protein backbone""", end=" ", file=open(TLEaPinput, 'a'))
                    if (HemeType == "c"):
                        print(f"""
 bond  ox.{CYSb}.CA   ox.{HEM}.CBB2
 bond  ox.{CYSc}.CA   ox.{HEM}.CBC1

 bond  rd.{CYSb}.CA   rd.{HEM}.CBB2
 bond  rd.{CYSc}.CA   rd.{HEM}.CBC1""", file=open(TLEaPinput, 'a'))

                    print(f"""
#Bond propionic acids to heme
 bond  ox.{HEM}.C2A   ox.{HEM+1}.CA
 bond  ox.{HEM}.C3D   ox.{HEM+2}.CA

 bond  rd.{HEM}.C2A   rd.{HEM+1}.CA
 bond  rd.{HEM}.C3D   rd.{HEM+2}.CA""", file=open(TLEaPinput, 'a'))

            print(f"""
# Save topology and coordinate files
 saveamberparm  ox o{SelHEM[idxs]}.prmtop o{SelHEM[idxs]}.rst7
 saveamberparm  rd r{SelHEM[idxs]}.prmtop r{SelHEM[idxs]}.rst7

quit""", file=open(TLEaPinput, 'a'))


            print(f"   Running TLEaP to generate o{SelHEM[idxs]}.prmtop, o{SelHEM[idxs]}.rst7, r{SelHEM[idxs]}.prmtop, and r{SelHEM[idxs]}.rst7 ...")
            subprocess.run(f"tleap -s -f {TLEaPinput} > GenerateRedoxStateTopologiesForHem{SelHEM[idxs]}.log", shell=True)

            if (os.path.isfile(f"o{SelHEM[idxs]}.prmtop") == True) and (os.path.isfile(f"o{SelHEM[idxs]}.rst7") == True):
               print(f"    TLEaP successfully generated o{SelHEM[idxs]}.prmtop and o{SelHEM[idxs]}.rst7")
            elif (os.path.isfile(f"o{SelHEM[idxs]}.prmtop") == False) or (os.path.isfile(f"o{SelHEM[idxs]}.rst7") == False):
                sys.exit(f"""     
 TLEaP failed to generated o{SelHEM[idxs]}.prmtop and/or o{SelHEM[idxs]}.rst7
 Please check GenerateRedoxStateTopologiesForHem{SelHEM[idxs]}.log 
 to diagnose the problem. \n""")

            if (os.path.isfile(f"r{SelHEM[idxs]}.prmtop") == True) and (os.path.isfile(f"r{SelHEM[idxs]}.rst7") == True):
                print(f"     TLEaP successfully generated r{SelHEM[idxs]}.prmtop and r{SelHEM[idxs]}.rst7")
            elif (os.path.isfile(f"r{SelHEM[idxs]}.prmtop") == False) or (os.path.isfile(f"r{SelHEM[idxs]}.rst7") == False):
                sys.exit(f"""     
 TLEaP failed to generated r{SelHEM[idxs]}.prmtop and/or r{SelHEM[idxs]}.rst7
 Please check GenerateRedoxStateTopologiesForHem{SelHEM[idxs]}.log 
 to diagnose the problem. \n""")

            idxs+=1
################################################################################################################################################
