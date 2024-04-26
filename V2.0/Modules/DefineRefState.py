################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

def DefineRefState(ForceFieldDir, OutPrefix):

    while True:
        ChooseRef = input(" Should the reference state be all-hemes oxidized (ox) or all-hemes reduced (red)? ")

        if (ChooseRef == "OX") or (ChooseRef == "Ox") or (ChooseRef == "ox") or (ChooseRef == "O") or (ChooseRef == "o"):
            SelRefRedoxState = "O"

            print("""
 mol new min.pdb

 #Change resname of proximal His
 set PHR [atomselect top "resname PHR"]
 $PHR set resname PHO
 set PMR [atomselect top "resname PMR"]
 $PMR set resname PMO
 set FHR [atomselect top "resname FHR"]
 $FHR set resname FHO
 set FMR [atomselect top "resname FMR"]
 $FMR set resname FMO

 #Change resname of distal His/Met
 set DHR [atomselect top "resname DHR"]
 $DHR set resname DHO
 set DMR [atomselect top "resname DMR"]
 $DMR set resname DMO
 set RHR [atomselect top "resname RHR"]
 $RHR set resname RHO
 set RMR [atomselect top "resname RMR"]
 $RMR set resname RMO
  
 #Change resname of the heme group
 set HCR [atomselect top "resname HCR"]
 $HCR set resname HCO
 set MCR [atomselect top "resname MCR"]
 $MCR set resname MCO
 set HBR [atomselect top "resname HBR"]
 $HBR set resname HBO
 set MBR [atomselect top "resname MBR"]
 $MBR set resname MBO

 set all [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
 $all writepdb RefState.pdb

 exit""", file=open('SetRefRedoxState.tcl', 'w'))

            subprocess.run("vmd -e SetRefRedoxState.tcl > SetRefRedoxState.log", shell=True)
            break

        elif (ChooseRef == "RED") or (ChooseRef == "Red") or (ChooseRef == "red") or (ChooseRef == "R") or (ChooseRef == "r"):
            SelRefRedoxState = "R"

            print("""
 mol new min.pdb

 #Change resname of proximal His
 set PHO [atomselect top "resname PHO"]
 $PHO set resname PHR
 set PMO [atomselect top "resname PMO"]
 $PMO set resname PMR
 set FHO [atomselect top "resname FHO"]
 $FHO set resname FHR
 set FMO [atomselect top "resname FMO"]
 $FMO set resname FMR

 #Change resname of distal His/Met
 set DHO [atomselect top "resname DHO"]
 $DHO set resname DHR
 set DMO [atomselect top "resname DMO"]
 $DMO set resname DMR
 set RHO [atomselect top "resname RHO"]
 $RHO set resname RHR
 set RMO [atomselect top "resname RMO"]
 $RMO set resname RMR
  
 #Change resname of the heme group
 set HCO [atomselect top "resname HCO"]
 $HCO set resname HCR
 set MCO [atomselect top "resname MCO"]
 $MCO set resname MCR
 set HBO [atomselect top "resname HBO"]
 $HBO set resname HBR
 set MBO [atomselect top "resname MBO"]
 $MBO set resname MBR

 set all [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]
 $all writepdb RefState.pdb

 exit""", file=open('SetRefRedoxState.tcl', 'w'))

            subprocess.run("vmd -e SetRefRedoxState.tcl > SetRefRedoxState.log", shell=True)
            break
        else:
            print(" Sorry, I didn't understand your response")

    if (os.path.isfile(f"RefState.pdb") == False):
        sys.exit("""
 Something went wrong running VMD to generate the 
 reference state PDB. Please check SetRefRedoxState.tcl
 and SetRefRedoxState.log for errors. \n""")
    elif (os.path.isfile(f"RefState.pdb") == True):
        Format = f"sed -i '/OXT/a TER' RefState.pdb"
        subprocess.run(Format, shell=True)

    if (os.path.isfile(f"ResIndexing.txt") == False):
        sys.exit("""
 ResIndexing.txt is missing!
 I, unfortunately, do not know how to
 proceed without this file. 

 Please check CreateResIndexing.log and 
 and SetupStructure.log for problems that 
 may have occured in the previous steps.\n""") 
    elif (os.path.isfile(f"ResIndexing.txt") == True):

        idx=0
        Count_c_HH = 0 
        Count_c_HM = 0
        Count_b_HH = 0
        Count_b_HM = 0
         
        with open("ResIndexing.txt") as fp:
            LineNumErr = []
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
                else:
                    LineNumErr.append(idx+1)
                    #str(LineNumErr)+str(f", {idx}")
                idx+=1

        if ( LineNumErr != [] ):
            sys.exit(f"""
 Something is wrong in ResIndexing.txt on line(s): {LineNumErr}.
 Please check and correct the file before re-running BioDC.
 Please save the corrected file as CorrectedResIndexing.txt.
 BioDC will automatically detected this file and use it when
 you re-run module #1.""")

        elif ( LineNumErr == [] ):
            print("""
# Load parameters
 source leaprc.constph
 source leaprc.conste
 source leaprc.gaff
 source leaprc.water.tip3p

 addAtomTypes {""", end=' ', file=open('SetupRefStateTopology.in', 'w'))

            if ( Count_b_HH != 0) and (SelRefRedoxState == "O"):
                print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_b_HH != 0) and (SelRefRedoxState == "R"):
                print("""
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_b_HM != 0) and (SelRefRedoxState == "O"):
                print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "N" "sp3" }  #Oxidized
        { "W2"  "S" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_b_HM != 0) and (SelRefRedoxState == "R"):
                print("""
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "N" "sp3" }  #Reduced
        { "X2"  "S" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_c_HM != 0) and (SelRefRedoxState == "O"):
                print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_c_HM != 0) and (SelRefRedoxState == "R"):
                print("""
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_c_HH != 0) and (SelRefRedoxState == "O"):
                print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_c_HH != 0) and (SelRefRedoxState == "R"):
                print("""
        { "M8"  "Fe" "sp3" } #M8&T1-T6:
        { "T1"  "N" "sp3" }  #Reduced
        { "T2"  "N" "sp3" }  #His-His
        { "T3"  "N" "sp3" }  #Ligated
        { "T4"  "N" "sp3" }  #c-Heme
        { "T5"  "N" "sp3" }
        { "T6"  "N" "sp3" }""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            print("""\n } """, file=open('SetupRefStateTopology.in', 'a'))

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
""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_b_HH != 0 ) and (SelRefRedoxState == "O"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
            if ( Count_b_HH != 0 ) and (SelRefRedoxState == "R"):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
            if ( Count_b_HM != 0 ) and (SelRefRedoxState == "O"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
            if ( Count_b_HM != 0 ) and (SelRefRedoxState == "R"):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

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
""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if ( Count_c_HH != 0 ) and (SelRefRedoxState == "O"):
                while True:
                    FFchoice = input("""
 Do you wish to use the previously published set of atomic 
 partial charges for oxidized bis-histidine-ligated c-type 
 hemes developed by Henriques et al.[1] or the set developed
 for BioDC by Guberman-Pfeffer[2]?

    Reference:
      [1] Henriques, J.; Costa, P. J.; Calhorda, M. J.; 
          Machuqueiro, M. Charge Parametrization of the
          DvH-c3 Heme Group: Validation Using Constant-(pH,E) 
          Molecular Dynamics Simulations. J. Phys. Chem. B 
          2013, 117 (1), 70–82. 

      [2] In preparation. RESP charges were commputed using the 
          MK scheme with the B3LYP approximate density functional 
          and a mixed basis set (LANL2TZ(f) for Fe, and 6-31G(d) for
          second row elements.)
 
 Please Henriques (H) or Guberman-Pfeffer (GP) charge sets: """)
                    if (FFchoice == "HENRIQUES") or (FFchoice == "Henriques") or (FFchoice == "henriques") or (FFchoice == "H") or (FFchoice == "h"):
                        print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
                        break
                    elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                        print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
                        break
                    else:
                        print("Sorry, I didn't understand your response.")
            if ( Count_c_HH != 0 ) and (SelRefRedoxState == "R"):
                while True:
                    FFchoice = input("""
 Do you wish to use the previously published set of atomic 
 partial charges for reduced bis-histidine-ligated c-type 
 hemes developed by Henriques et al.[1] or the set developed
 for BioDC by Guberman-Pfeffer[2]?

    Reference:
      [1] Henriques, J.; Costa, P. J.; Calhorda, M. J.; 
          Machuqueiro, M. Charge Parametrization of the
          DvH-c3 Heme Group: Validation Using Constant-(pH,E) 
          Molecular Dynamics Simulations. J. Phys. Chem. B 
          2013, 117 (1), 70–82. 

      [2] In preparation. RESP charges were commputed using the 
          MK scheme with the B3LYP approximate density functional 
          and a mixed basis set (LANL2TZ(f) for Fe, and 6-31G(d) for
          second row elements.)
 
 Please Henriques (H) or Guberman-Pfeffer (GP) charge sets: """)
                    if (FFchoice == "HENRIQUES") or (FFchoice == "Henriques") or (FFchoice == "henriques") or (FFchoice == "H") or (FFchoice == "h"):
                        print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
                        break
                    elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                        print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
                        break
                    else:
                        print("Sorry, I didn't understand your response.")

            if ( Count_c_HM != 0 ) and (SelRefRedoxState == "O"):
                print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))
            if ( Count_c_HM != 0 ) and (SelRefRedoxState == "R"):
                print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            print(f"""

# Load PDB
 mol = loadpdb RefState.pdb""", file=open('SetupRefStateTopology.in', 'a'))

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
                        print(f" bond mol.{DisulfPairID[sbi][0]}.SG mol.{DisulfPairID[sbi][1]}.SG", file=open('SetupRefTopology.in', 'a'))

            idx=0
            with open("ResIndexing.txt") as fp:
                NumHEC = len(fp.readlines())
                HemID = [0]*NumHEC
                fp.seek(0)

                Lines = fp.readlines()
                for line in Lines:
                    EntryLength = len(line.strip().split(" "))
                    HemID[idx] = line.strip().split(" ")[-3]
                    HemeType = line.strip().split(" ")[-2]
                    AxLigType = line.strip().split(" ")[-1]

                    if ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HH" ):
                        CYSb = int(line.strip().split(" ")[0])
                        CYSc = int(line.strip().split(" ")[1])
                        Ligp = int(line.strip().split(" ")[2])
                        Ligd = int(line.strip().split(" ")[3])
                        HEM = int(line.strip().split(" ")[5])

                        print(f"""
#------------------------------------------------------------
#For heme {HEM}:

#Bond ligating atoms to Fe center
 bond mol.{Ligp}.NE2   mol.{HEM}.FE
 bond mol.{Ligd}.NE2   mol.{HEM}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond mol.{Ligp-1}.C   mol.{Ligp}.N
 bond mol.{Ligp}.C   mol.{Ligp+1}.N
 bond mol.{Ligd-1}.C   mol.{Ligd}.N
 bond mol.{Ligd}.C   mol.{Ligd+1}.N

#Bond heme thioethers to protein backbone
 bond mol.{CYSb}.CA   mol.{HEM}.CBB2
 bond mol.{CYSc}.CA   mol.{HEM}.CBC1

#Bond propionic acids to heme
 bond mol.{HEM}.C2A   mol.{HEM+1}.CA
 bond mol.{HEM}.C3D   mol.{HEM+2}.CA
#------------------------------------------------------------""", file=open('SetupRefStateTopology.in', 'a'))

                    elif ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HM" ):
                        CYSb = int(line.strip().split(" ")[0])
                        CYSc = int(line.strip().split(" ")[1])
                        Ligp = int(line.strip().split(" ")[2])
                        Ligd = int(line.strip().split(" ")[3])
                        HEM = int(line.strip().split(" ")[5])

                        print(f"""
#------------------------------------------------------------
#For heme {HEM}:

#Bond ligating atoms to Fe center
 bond mol.{Ligp}.NE2   mol.{HEM}.FE
 bond mol.{Ligd}.SD   mol.{HEM}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond mol.{Ligp-1}.C   mol.{Ligp}.N
 bond mol.{Ligp}.C   mol.{Ligp+1}.N
 bond mol.{Ligd-1}.C   mol.{Ligd}.N
 bond mol.{Ligd}.C   mol.{Ligd+1}.N

#Bond heme thioethers to protein backbone
 bond mol.{CYSb}.CA   mol.{HEM}.CBB2
 bond mol.{CYSc}.CA   mol.{HEM}.CBC1

#Bond propionic acids to heme
 bond mol.{HEM}.C2A   mol.{HEM+1}.CA
 bond mol.{HEM}.C3D   mol.{HEM+2}.CA
#------------------------------------------------------------""", file=open('SetupRefStateTopology.in', 'a'))

                    elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HH" ):
                        Ligp = int(line.strip().split(" ")[0])
                        Ligd = int(line.strip().split(" ")[1])
                        HEM = int(line.strip().split(" ")[3])

                        print(f"""
#------------------------------------------------------------
#For heme {HEM}:

#Bond ligating atoms to Fe center
 bond mol.{Ligp}.NE2   mol.{HEM}.FE
 bond mol.{Ligd}.NE2   mol.{HEM}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond mol.{Ligp-1}.C   mol.{Ligp}.N
 bond mol.{Ligp}.C   mol.{Ligp+1}.N
 bond mol.{Ligd-1}.C   mol.{Ligd}.N
 bond mol.{Ligd}.C   mol.{Ligd+1}.N

#Bond propionic acids to heme
 bond mol.{HEM}.C2A   mol.{HEM+1}.CA
 bond mol.{HEM}.C3D   mol.{HEM+2}.CA
#------------------------------------------------------------""", file=open('SetupRefStateTopology.in', 'a'))

                    elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HM" ):
                        Ligp = int(line.strip().split(" ")[0])
                        Ligd = int(line.strip().split(" ")[1])
                        HEM = int(line.strip().split(" ")[3])

                        print(f"""
#------------------------------------------------------------
#For heme {HEM}:

#Bond ligating atoms to Fe center
 bond mol.{Ligp}.NE2   mol.{HEM}.FE
 bond mol.{Ligd}.SD   mol.{HEM}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond mol.{Ligp-1}.C   mol.{Ligp}.N
 bond mol.{Ligp}.C   mol.{Ligp+1}.N
 bond mol.{Ligd-1}.C   mol.{Ligd}.N
 bond mol.{Ligd}.C   mol.{Ligd+1}.N

#Bond propionic acids to heme
 bond mol.{HEM}.C2A   mol.{HEM+1}.CA
 bond mol.{HEM}.C3D   mol.{HEM+2}.CA
#------------------------------------------------------------""", file=open('SetupRefStateTopology.in', 'a'))

                    else:
                        print(f""" 
 #Problem defining bond definitions for heme {HemID[idx]}.
 #Please inspect ResIndexing.txt for missing or incomplete entries.""")

                    idx+=1

            print(f"""
# Save topology and coordinate files
saveamberparm mol RefState.prmtop RefState.rst7

quit""", end=" ", file=open('SetupRefStateTopology.in', 'a'))

            if (os.path.isfile(f"SetupRefStateTopology.in") == True):
                print("""
 Using TLEaP to build the reference state topology... 
    """, end=" ")
                subprocess.run("tleap -s -f SetupRefStateTopology.in > SetupRefStateTopology.log", shell=True)

            if (os.path.isfile(f"RefState.prmtop") == True) and (os.path.isfile(f"RefState.rst7") == True):
                print("""
 TLEaP finished successfully! """)
            if (os.path.isfile(f"RefState.prmtop") == False) or (os.path.isfile(f"RefState.rst7") == False):
                sys.exit("""
 TLEaP failed!
 Please inspect SetupRefStateTopology.in and SetupRefStateTopology.log for problems. \n""")

    try:
        return SelRefRedoxState, FFchoice
    except UnboundLocalError:
        FFchoice=''
        return SelRefRedoxState, FFchoice
################################################################################################################################################
