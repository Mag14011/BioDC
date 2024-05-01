################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import SelectpHActiveSites
################################################################################################################################################

def ReBuildStructure(ForceFieldDir, PDB, SelCpH, InputDict, LaunchDir):

    if (os.path.isfile(f"ResIndexing.txt") == False):
        sys.exit("""
 ResIndexing.txt is missing!
 I, unfortunately, do not know how to
 proceed without this file. 

 Please check CreateResIndexing.log and 
 and SetupStructure.log for problems that 
 may have occured in the previous steps.\n""") 

    if (os.path.isfile(f"ResIndexing.txt") == True):
        print("""
 Now, we need to stitch the edited PDBs of the protein, hemes, 
 and heme propionic acid groups into a single PDB. Then, this 
 re-constructed PDB of the cytochrome will be processed with 
 TLEaP of the AmberTools package to generate topology and 
 coordinate files.
    """, end=" ")
   
        if ("OutPrefix" in InputDict):
            OutPrefix = InputDict["OutPrefix"]
        else:
            OutPrefix = input(" \n Prefix for output parm/rst7: ")

        print(f"OutPrefix = {OutPrefix}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

        subprocess.run(f"cat prot.pdb > temp.pdb", shell=True)

        idx=0
        LineNumErr = ''
        CountHCO = 0; CountHCR = 0
        CountHBO = 0; CountHBR = 0
        CountMCO = 0; CountMCR = 0
        CountMBO = 0; CountMBR = 0

        with open("ResIndexing.txt") as fp:
            NumHEC = len(fp.readlines())
            HemID = [0]*NumHEC
            fp.seek(0)
 
            Lines = fp.readlines()
            for line in Lines:
                HemID[idx] = line.strip().split(" ")[-3]

                if (os.path.isfile(f"HCO{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountHCO+=1
                    print(f"   Adding HCO{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat HCO{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for an oxidized His-His ligated c-type heme! \n")
                    else:
                        print("   ! Missing the forcefield files for an oxidized His-His ligated c-type heme! \n")
                elif (os.path.isfile(f"HCR{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountHCR+=1
                    print(f"   Adding HCR{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat HCR{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for a  reduced  His-His ligated c-type heme! \n")
                    else:
                        print("   ! Missing the forcefield files for a  reduced  His-His ligated c-type heme! \n")
                elif (os.path.isfile(f"MCO{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountMCO+=1
                    print(f"   Adding MCO{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat MCO{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for an oxidized His-Met ligated c-type heme!\n")
                    else:
                        print("   ! Missing the forcefield files for an oxidized His-Met ligated c-type heme!\n")
                elif (os.path.isfile(f"MCR{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountMCR+=1
                    print(f"   Adding MCR{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat MCR{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for a  reduced  His-Met ligated c-type heme! \n")
                    else:
                        print("   ! Missing the forcefield files for a  reduced  His-Met ligated c-type heme! \n")
                elif (os.path.isfile(f"HBO{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountHBO+=1
                    print(f"   Adding HBO{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat HBO{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for an oxidized His-His ligated b-type heme!\n")
                    else:
                        print("   ! Missing the forcefield files for an oxidized His-His ligated b-type heme!\n")
                elif (os.path.isfile(f"HBR{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountHBR+=1
                    print(f"   Adding HBR{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat HBR{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for a  reduced  His-His ligated b-type heme!\n")
                    else:
                        print("   ! Missing the forcefield files for a  reduced  His-His ligated b-type heme!\n")
                elif (os.path.isfile(f"MBO{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountMBO+=1
                    print(f"   Adding MBO{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat MBO{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for an oxidized His-Met ligated b-type heme!\n")
                    else:
                        print("   ! Missing the forcefield files for an oxidized His-Met ligated b-type heme!\n")
                elif (os.path.isfile(f"MBR{HemID[idx]}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == True) and (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == True):
                    CountMBR+=1
                    print(f"   Adding MBR{HemID[idx]} and its propionic acid groups.""")
                    subprocess.run(f"cat MBR{HemID[idx]}.pdb PRN{int(HemID[idx])+1}.pdb PRN{int(HemID[idx])+2}.pdb >> temp.pdb", shell=True)

                    if (os.path.isfile(f"{ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib") == True) and (os.path.isfile(f"{ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod") == True):
                        print("   √ Found   the forcefield files for a  reduced  His-Met ligated b-type heme!\n")
                    else:
                        print("   ! Missing the forcefield files for a  reduced  His-Met ligated b-type heme!\n")
                elif (os.path.isfile(f"HCO{HemID[idx]}.pdb") == False) or (os.path.isfile(f"HCR{HemID[idx]}.pdb") == False) or (os.path.isfile(f"MCO{HemID[idx]}.pdb") == False) or (os.path.isfile(f"MCR{HemID[idx]}.pdb") == False) or (os.path.isfile(f"HBO{HemID[idx]}.pdb") == False) or (os.path.isfile(f"HBR{HemID[idx]}.pdb") == False) or (os.path.isfile(f"MBO{HemID[idx]}.pdb") == False) or (os.path.isfile(f"MBR{HemID[idx]}.pdb") == False):
                    sys.exit(f""" 
 The PDB for heme-{HemID[idx]} is missing!
 Something went wrong in the previous step.

 Please check the log file and correct the 
 error before re-running this module.

 I apologize for the inconvenience!""")
                elif (os.path.isfile(f"PRN{int(HemID[idx])+1}.pdb") == False) or (os.path.isfile(f"PRN{int(HemID[idx])+2}.pdb") == False):
                    sys.exit(f""" 
 The PDB for one or both of the propionic acid groups from heme-{HemID[idx]} is missing!
 Something went wrong in the previous step.

 Please check the log file and correct the 
 error before re-running this module.

 I apologize for the inconvenience!""")

                idx+=1

        subprocess.run(f"grep -v CRYST1 temp.pdb | grep -v END > {OutPrefix}-{PDB}.pdb", shell=True)
        Format = f"sed -i '/OXT/a TER' {OutPrefix}-{PDB}.pdb"
        subprocess.run(Format, shell=True)

        print(f" Heme Type \t Ligation state \t Redox state \t Count")
        print(f"     c     \t      His-His   \t   oxidized  \t   {CountHCO}")
        print(f"     c     \t      His-His   \t    reduced  \t   {CountHCR}")
        print(f"     c     \t      His-Met   \t   oxidized  \t   {CountMCO}")
        print(f"     c     \t      His-Met   \t    reduced  \t   {CountMCR}")
        print(f"     b     \t      His-His   \t   oxidized  \t   {CountHBO}")
        print(f"     b     \t      His-His   \t    reduced  \t   {CountHBR}")
        print(f"     b     \t      His-Met   \t   oxidized  \t   {CountMBO}")
        print(f"     b     \t      His-Met   \t    reduced  \t   {CountMBR}")

        print("""
# Load parameters
 source leaprc.constph
 source leaprc.conste
 source leaprc.gaff
 source leaprc.water.tip3p

 addAtomTypes {""", end=' ', file=open('tleap.in', 'w'))

        if ( CountHBO != 0 ):
            print("""
        { "M1"  "Fe" "sp3" } #M1&Y1-Y6:
        { "Y1"  "N" "sp3" }  #Oxidized
        { "Y2"  "N" "sp3" }  #His-His
        { "Y3"  "N" "sp3" }  #Ligated
        { "Y4"  "N" "sp3" }  #b-Heme
        { "Y5"  "N" "sp3" }
        { "Y6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHBR != 0 ):
            print("""
        { "M2"  "Fe" "sp3" } #M2&Z1-Z6:
        { "Z1"  "N" "sp3" }  #Reduced
        { "Z2"  "N" "sp3" }  #His-His
        { "Z3"  "N" "sp3" }  #Ligated
        { "Z4"  "N" "sp3" }  #b-Heme
        { "Z5"  "N" "sp3" }
        { "Z6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountMBO != 0 ):
            print("""
        { "M3"  "Fe" "sp3" } #M3&W1-W6:
        { "W1"  "N" "sp3" }  #Oxidized
        { "W2"  "S" "sp3" }  #His-Met
        { "W3"  "N" "sp3" }  #Ligated
        { "W4"  "N" "sp3" }  #b-Heme
        { "W5"  "N" "sp3" }
        { "W6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountMBR != 0 ):
            print("""
        { "M4"  "Fe" "sp3" } #M4&X1-X6:
        { "X1"  "N" "sp3" }  #Reduced
        { "X2"  "S" "sp3" }  #His-Met
        { "X3"  "N" "sp3" }  #Ligated
        { "X4"  "N" "sp3" }  #b-Heme
        { "X5"  "N" "sp3" }
        { "X6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountMCO != 0 ):
            print("""
        { "M5"  "Fe" "sp3" } #M5&U1-U6:
        { "U1"  "N" "sp3" }  #Oxidized
        { "U2"  "S" "sp3" }  #His-Met
        { "U3"  "N" "sp3" }  #Ligated
        { "U4"  "N" "sp3" }  #c-Heme
        { "U5"  "N" "sp3" }
        { "U6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountMCR != 0 ):
            print("""
        { "M6"  "Fe" "sp3" } #M6&V1-V6:
        { "V1"  "N" "sp3" }  #Reduced
        { "V2"  "S" "sp3" }  #His-Met
        { "V3"  "N" "sp3" }  #Ligated
        { "V4"  "N" "sp3" }  #c-Heme
        { "V5"  "N" "sp3" }
        { "V6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHCO != 0 ):
            print("""
        { "M7"  "Fe" "sp3" } #M7&S1-S6:
        { "S1"  "N" "sp3" }  #Oxidized
        { "S2"  "N" "sp3" }  #His-His
        { "S3"  "N" "sp3" }  #Ligated
        { "S4"  "N" "sp3" }  #c-Heme
        { "S5"  "N" "sp3" }
        { "S6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHCR != 0 ):
            print("""
        { "M8"  "Fe" "sp3" } #M8&T1-T6:
        { "T1"  "N" "sp3" }  #Reduced
        { "T2"  "N" "sp3" }  #His-His
        { "T3"  "N" "sp3" }  #Ligated
        { "T4"  "N" "sp3" }  #c-Heme
        { "T5"  "N" "sp3" }
        { "T6"  "N" "sp3" }""", end=" ", file=open('tleap.in', 'a'))

        print("""\n } """, file=open('tleap.in', 'a'))

        if ( CountHBO != 0 ) or ( CountHBR != 0 ) or ( CountMBO != 0 ) or ( CountMBR != 0 ):
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
""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHBO != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
        if ( CountHBR != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_b-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
        if ( CountMBO != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
        if ( CountMBR != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_b-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_b-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHCO != 0 ) or ( CountHCR != 0 ) or ( CountMCO != 0 ) or ( CountMCR != 0 ):
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
""", end=" ", file=open('tleap.in', 'a'))

        if ( CountHCO != 0 ):
            while True:
                if ("FFchoice" in InputDict):
                    FFchoice = InputDict["FFchoice"]
                else:
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
 loadoff {ForceFieldDir}/Henriques_Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
                    break
                elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
                    break
                else:
                    print("Sorry, I didn't understand your response.")

        if ( CountHCR != 0 ):
            while True:
                if ("FFchoice" in InputDict):
                    FFchoice = InputDict["FFchoice"]
                else:
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
 loadoff {ForceFieldDir}/Henriques_Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
                    break
                elif (FFchoice == "GUBERMAN-PFEFFER") or (FFchoice == "Guberman-Pfeffer") or (FFchoice == "guberman-pfeffer") or (FFchoice == "GP") or (FFchoice == "gp") or (FFchoice == "Gp") or (FFchoice == "gP"):
                    print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisHisLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisHisLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
                    break
                else:
                    print("Sorry, I didn't understand your response.")

        if ( CountMCO != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Oxidized_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Oxidized_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))
        if ( CountMCR != 0 ):
            print(f"""
 loadamberparams {ForceFieldDir}/Reduced_HisMetLigated_c-heme.frcmod
 loadoff {ForceFieldDir}/Reduced_HisMetLigated_c-heme_RESP.lib""", end=" ", file=open('tleap.in', 'a'))

        print(f"""

# Load PDB
 {OutPrefix} = loadpdb {OutPrefix}-{PDB}.pdb""", file=open('tleap.in', 'a'))
      
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
            Lines = fp.readlines()
            for line in Lines:
                EntryLength = len(line.strip().split(" "))
                HemeType = line.strip().split(" ")[-2]
                AxLigType = line.strip().split(" ")[-1]

                if ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HH" ):
                    idx+=1
                    CYSb = int(line.strip().split(" ")[0])
                    CYSc = int(line.strip().split(" ")[1])
                    Ligp = int(line.strip().split(" ")[2])
                    Ligd = int(line.strip().split(" ")[3])
                    HECHHfnl = int(line.strip().split(" ")[5])

                    print(f"""
#------------------------------------------------------------
#For heme {HECHHfnl}:

#Bond ligating atoms to Fe center
 bond {OutPrefix}.{Ligp}.NE2   {OutPrefix}.{HECHHfnl}.FE
 bond {OutPrefix}.{Ligd}.NE2   {OutPrefix}.{HECHHfnl}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond {OutPrefix}.{Ligp-1}.C   {OutPrefix}.{Ligp}.N
 bond {OutPrefix}.{Ligp}.C   {OutPrefix}.{Ligp+1}.N
 bond {OutPrefix}.{Ligd-1}.C   {OutPrefix}.{Ligd}.N
 bond {OutPrefix}.{Ligd}.C   {OutPrefix}.{Ligd+1}.N

#Bond heme thioethers to protein backbone
 bond {OutPrefix}.{CYSb}.CA   {OutPrefix}.{HECHHfnl}.CBB2
 bond {OutPrefix}.{CYSc}.CA   {OutPrefix}.{HECHHfnl}.CBC1

#Bond propionic acids to heme
 bond {OutPrefix}.{HECHHfnl}.C2A   {OutPrefix}.{HECHHfnl+1}.CA
 bond {OutPrefix}.{HECHHfnl}.C3D   {OutPrefix}.{HECHHfnl+2}.CA
#------------------------------------------------------------""", file=open('tleap.in', 'a'))

                elif ( EntryLength == 8 ) and ( HemeType == "c") and ( AxLigType == "HM" ):
                    idx+=1
                    CYSb = int(line.strip().split(" ")[0])
                    CYSc = int(line.strip().split(" ")[1])
                    Ligp = int(line.strip().split(" ")[2])
                    Ligd = int(line.strip().split(" ")[3])
                    HECHMfnl = int(line.strip().split(" ")[5])

                    print(f"""
#------------------------------------------------------------
#For heme {HECHMfnl}:

#Bond ligating atoms to Fe center
 bond {OutPrefix}.{Ligp}.NE2   {OutPrefix}.{HECHMfnl}.FE
 bond {OutPrefix}.{Ligd}.SD   {OutPrefix}.{HECHMfnl}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond {OutPrefix}.{Ligp-1}.C   {OutPrefix}.{Ligp}.N
 bond {OutPrefix}.{Ligp}.C   {OutPrefix}.{Ligp+1}.N
 bond {OutPrefix}.{Ligd-1}.C   {OutPrefix}.{Ligd}.N
 bond {OutPrefix}.{Ligd}.C   {OutPrefix}.{Ligd+1}.N

#Bond heme thioethers to protein backbone
 bond {OutPrefix}.{CYSb}.CA   {OutPrefix}.{HECHMfnl}.CBB2
 bond {OutPrefix}.{CYSc}.CA   {OutPrefix}.{HECHMfnl}.CBC1

#Bond propionic acids to heme
 bond {OutPrefix}.{HECHMfnl}.C2A   {OutPrefix}.{HECHMfnl+1}.CA
 bond {OutPrefix}.{HECHMfnl}.C3D   {OutPrefix}.{HECHMfnl+2}.CA
#------------------------------------------------------------""", file=open('tleap.in', 'a'))

                elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HH" ):
                    idx+=1
                    Ligp = int(line.strip().split(" ")[0])
                    Ligd = int(line.strip().split(" ")[1])
                    HEBHHfnl = int(line.strip().split(" ")[3])

                    print(f"""
#------------------------------------------------------------
#For heme {HEBHHfnl}:

#Bond ligating atoms to Fe center
 bond {OutPrefix}.{Ligp}.NE2   {OutPrefix}.{HEBHHfnl}.FE
 bond {OutPrefix}.{Ligd}.NE2   {OutPrefix}.{HEBHHfnl}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond {OutPrefix}.{Ligp-1}.C   {OutPrefix}.{Ligp}.N
 bond {OutPrefix}.{Ligp}.C   {OutPrefix}.{Ligp+1}.N
 bond {OutPrefix}.{Ligd-1}.C   {OutPrefix}.{Ligd}.N
 bond {OutPrefix}.{Ligd}.C   {OutPrefix}.{Ligd+1}.N

#Bond propionic acids to heme
 bond {OutPrefix}.{HEBHHfnl}.C2A   {OutPrefix}.{HEBHHfnl+1}.CA
 bond {OutPrefix}.{HEBHHfnl}.C3D   {OutPrefix}.{HEBHHfnl+2}.CA
#------------------------------------------------------------""", file=open('tleap.in', 'a'))

                elif ( EntryLength == 6 ) and ( HemeType == "b") and ( AxLigType == "HM" ):
                    idx+=1
                    Ligp = int(line.strip().split(" ")[0])
                    Ligd = int(line.strip().split(" ")[1])
                    HEBHMfnl = int(line.strip().split(" ")[3])

                    print(f"""
#------------------------------------------------------------
#For heme {HEBHMfnl}:

#Bond ligating atoms to Fe center
 bond {OutPrefix}.{Ligp}.NE2   {OutPrefix}.{HEBHMfnl}.FE
 bond {OutPrefix}.{Ligd}.SD   {OutPrefix}.{HEBHMfnl}.FE

#Bond axially coordinated residues to preceeding and proceeding residues
 bond {OutPrefix}.{Ligp-1}.C   {OutPrefix}.{Ligp}.N
 bond {OutPrefix}.{Ligp}.C   {OutPrefix}.{Ligp+1}.N
 bond {OutPrefix}.{Ligd-1}.C   {OutPrefix}.{Ligd}.N
 bond {OutPrefix}.{Ligd}.C   {OutPrefix}.{Ligd+1}.N

#Bond propionic acids to heme
 bond {OutPrefix}.{HEBHMfnl}.C2A   {OutPrefix}.{HEBHMfnl+1}.CA
 bond {OutPrefix}.{HEBHMfnl}.C3D   {OutPrefix}.{HEBHMfnl+2}.CA
#------------------------------------------------------------""", file=open('tleap.in', 'a'))

                else:
                    print(f""" 
 #Problem defining bond definitions for heme {HemID[idx]}.
 #Please inspect ResIndexing.txt for missing or incomplete entries.""")

        try:
            print(f"FFchoice = {FFchoice}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
        except UnboundLocalError:
            pass

        while True:
            if ("SolvEnv" in InputDict):
                SolvEnv = InputDict["SolvEnv"]
            else:
                SolvEnv = input(""" 
 Should the structure be prepared with 
  an explicit or implicit solvent (explicit/implicit)? """)

            print(f"SolvEnv = {SolvEnv}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if (SolvEnv == "EXPLICIT") or (SolvEnv == "Explicit") or (SolvEnv == "explicit") or (SolvEnv == "E") or (SolvEnv == "e"):
                while True:
                    if ("BoxShape" in InputDict):
                        BoxShape = InputDict["BoxShape"]
                    else:
                        BoxShape = input(" Using a rectangular or an octahedral box (rec/oct)? ")

                    print(f"BoxShape = {BoxShape}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                    while True:
                        try:
                            if ("BufferSize" in InputDict):
                                BufferSize = int(InputDict["BufferSize"])
                                print(f"BufferSize = {BufferSize}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                            else:
                                BufferSize = int(input(" With how much of a solvent buffer (in angstroms)? "))
                                print(f"BufferSize = {BufferSize}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        except ValueError:
                            print(" Your entry needs to be an integer. \n")
                        else:
                            break

            
                    if (BoxShape == "RECTANGULAR") or (BoxShape == "Rectangular") or (BoxShape == "rectangular") or (BoxShape == "REC") or (BoxShape == "Rec") or (BoxShape == "rec") or (BoxShape == "R") or (BoxShape == "r"):
                        print("""
#Solvate
solvateBox %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                        break
                    elif (BoxShape == "OCTAHEDRAL") or (BoxShape == "Octahedral") or (BoxShape == "octahedral") or (BoxShape == "OCT") or (BoxShape == "Oct") or (BoxShape == "oct") or (BoxShape == "O") or (BoxShape == "o"): 
                        print("""

#Solvate
solvateOct %s TIP3PBOX %0d""" %(OutPrefix, BufferSize), end=" ", file=open('tleap.in', 'a'))
                        break
                    else:
                        print(" Sorry, I didn't understand your response.")

                while True:
                    try:
                        if ("NaCount" in InputDict):
                            NaCount = int(InputDict["NaCount"])
                            print(f"NaCount = {NaCount}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        else:
                            NaCount = int(input(" And how many Na+ ions; 0 = enough for charge neutrality? "))
                            print(f"NaCount = {NaCount}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    except ValueError:
                        print(" Your entry needs to be an integer. \n")
                    else:
                        break


                while True:
                    try:
                        if ("ClCount" in InputDict):
                            ClCount = int(InputDict["ClCount"])
                            print(f"ClCount = {ClCount}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                        else:
                            ClCount = int(input(" And how many Cl- ions; 0 = enough for charge neutrality? "))
                            print(f"ClCount = {ClCount}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))
                    except ValueError:
                        print(" Your entry needs to be an integer. \n")
                    else:
                        break


                print("""

#Add ions
addions %s Na+ %0d""" %(OutPrefix, NaCount), end=" ", file=open('tleap.in', 'a'))
                print("""
addions %s Cl- %0d""" %(OutPrefix, ClCount), end=" ", file=open('tleap.in', 'a'))

                break
            elif (SolvEnv == "IMPLICIT") or (SolvEnv == "Implicit") or (SolvEnv == "implicit") or (SolvEnv == "I") or (SolvEnv == "i"):
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
    
        if (SelCpH == "YES") or (SelCpH == "Yes") or (SelCpH == "yes") or (SelCpH == "Y") or (SelCpH == "y"):
            print(""" 
 Now Generating the cpin file because you indicated that you want to run molecular dynamics
 with titratable residues. In order to take this step, the residues may need to be renumbered. 
 
 The foregoing structure preparation steps required you to place all 
 the heme residues from all the chains at the end of the PDB with 
 sequential numbering. The cpinutil.py program used to generate the
 cpin file needed for constant pH dynamics wants instead the residues 
 to come in the order of the connectivity; that is, the hemes of 
 chain A should come before any residue in chain B. 

 To oblige this different numbering convention, we'll use 
 CPPTRAJ of the AmberTools package to re-order the residues. 
 This process will write a new topology and coordinate file,

 Unfortunately, because the residue numbering may change, you will
 need to re-select which residues you want to be titratable.""")

            if (os.path.isfile(f"{OutPrefix}.prmtop") == True):
                print(f"""
parm {OutPrefix}.prmtop
trajin {OutPrefix}.rst7
fixatomorder parmout {OutPrefix}_reord.prmtop
trajout {OutPrefix}_reord.rst7 topresnum
trajout {OutPrefix}_reord.pdb  topresnum
run
quit
                """, file=open("ReorderRes.in", "w"))
                subprocess.run("cpptraj -i ReorderRes.in > ReorderRes.log 2> /dev/null", shell=True)
                SelASPIDs, SelGLUIDs, SelHISIDs, SelLYSIDs, SelTYRIDs, SelPRNIDs = SelectpHActiveSites.SelectpHActiveSites(f"{OutPrefix}_reord", "SecondPass", InputDict, LaunchDir)

            RESNAMES = " "
            if (len(SelASPIDs) != 0):
                RESNAMES+=str(" AS4 ")
            if (len(SelGLUIDs) != 0):
                RESNAMES+=str(" GL4 ")
            if (len(SelHISIDs) != 0):
                RESNAMES+=str(" HIP ")
            if (len(SelLYSIDs) != 0):
                RESNAMES+=str(" LYS ")
            if (len(SelTYRIDs) != 0):
                RESNAMES+=str(" TYR ")
            #if (len(PRNA) != 0) and (len(PRND) != 0):
            #    RESNAMES+=str(" PRN ")
            if (len(SelPRNIDs) != 0):
                RESNAMES+=str(" PRN ")

            RESID = " "
            word1 = f"CA  HIO"
            if (os.path.isfile(f"{OutPrefix}_reord.pdb") == True):

                with open(f"{OutPrefix}_reord.pdb", 'r') as fp:
                    lines = fp.readlines()
                    for line in lines:
                        if (line.find(word1) != -1):
                            idx1 = 0
                            for x in line:
                                if (idx1 >= 22) and (idx1 <= 25):
                                    RESID+=x
                                idx1+=1
                            RESID+=" "

                SelHISIDLIST=list(SelHISIDs.split())
                RESIDLIST=list(RESID.split())

                for i in SelHISIDLIST[:]:
                    if i in RESIDLIST:
                        RESIDLIST.remove(i)
                        SelHISIDLIST.remove(i)
            SelHisIDs = " ".join(SelHISIDLIST) 

            print(f"\n Command: \n  cpinutil.py -resnames {RESNAMES} -resnums {SelASPIDs} {SelGLUIDs} {SelHisIDs} {SelLYSIDs} {SelTYRIDs} {SelPRNIDs} -p {OutPrefix}_reord.prmtop -igb 2 -op {OutPrefix}_new.prmtop -o {OutPrefix}.cpin \n")
            subprocess.run(f"cpinutil.py -resnames {RESNAMES} -resnums {SelASPIDs} {SelGLUIDs} {SelHisIDs} {SelLYSIDs} {SelTYRIDs} {SelPRNIDs} -p {OutPrefix}_reord.prmtop -igb 2 -op {OutPrefix}_new.prmtop -o {OutPrefix}.cpin", shell=True)

    return OutPrefix, SolvEnv 

################################################################################################################################################
