################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################
def CreateResIndexing(PDB, InputDict, LaunchDir):

    PDBbasename = os.path.basename(PDB)

    if (os.path.isfile("CorrectedResIndexing.txt") == True):
        print(""" 
 Found CorrectedResIndexing.txt! 
 Copying CorrectedResIndexing.txt to ResIndexing.txt.""")
        subprocess.run("cp CorrectedResIndexing.txt ResIndexing.txt", shell=True)
    else:
        while True:
            if ("CreateResIndexingMethod" in InputDict):
                CreateResIndexingMethod = InputDict["CreateResIndexingMethod"]
            else:
                CreateResIndexingMethod = input("""
 We need to create a file (ResIndexing.txt) that identifies 
 the IDs of the Cys, His or Met residues bounded to the heme group.     
 Would you like to create it automatically or manually (auto/man)? """)

            print(f"CreateResIndexingMethod = {CreateResIndexingMethod}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

            if (CreateResIndexingMethod == "AUTO") or (CreateResIndexingMethod == "Auto") or (CreateResIndexingMethod == "auto") or (CreateResIndexingMethod == "A") or (CreateResIndexingMethod == "a"):
                print ("""
 The automated creation of ResIndexing.txt has two parts:
   1) Writing CreateResIndexing.tcl 
   2) Submitting the TCL script to Visual Molecular Dynamics (VMD)

 The TCL script identifies the Cys, His and Met residues within a 
 distance range of each heme group and assumes that these residues 
 are bonded to that particular heme.

 The distance from the heme is scanned in 0.1 Å increments 
 from the minimum until the maximum distance you specify.

 If only one axial ligand of type His is found, the script 
 then searches for a Met ligand within the same distance range.

 Recommended minimum and maximum distances are 2.0 and 4.0 Å, respectively """)

                while True:
                    try:
                        if ("DistMin" in InputDict):
                            DistMin = InputDict["DistMin"]
                        else:
                            DistMin = float(input("""
 What minimum distance threshold would you like to use? """))
                    except ValueError: 
                        print(" Your entry needs to be an integer. \n")
                    else:
                        break
                print(f"DistMin = {DistMin}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

                while True:
                    try:
                        if ("DistMax" in InputDict):
                            DistMax = InputDict["DistMax"]
                        else:
                            DistMax = float(input(""" What maximum distance threshold would you like to use? """))
                    except ValueError: 
                        print(" Your entry needs to be an integer. \n")
                    else:
                        break
                print(f"DistMax = {DistMax}", file=open(f"{LaunchDir}/InteractiveInput.txt", 'a'))

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

                print(f"""
 #--------------------------------------------------
 set DistMin {DistMin}
 set DistMax {DistMax}
 mol new {PDB}.pdb
 #--------------------------------------------------""", file=open('CreateResIndexing.tcl', 'w'))

                print("""
 #--------------------------------------------------
 #Create output files
 set out1 [open "ResIndexing.txt" w]
 set out2 [open "ResIndexing_HisHisHeme.txt" w]
 set out3 [open "ResIndexing_HisMetHeme.txt" w]
 #--------------------------------------------------

 #--------------------------------------------------
 #Identify the heme residue IDs
 set HEC [[atomselect top "resname HEC HEM HEH and name FE"] get resid]
 #--------------------------------------------------

 #--------------------------------------------------
 #Determine whether each heme is of the b- or c-type,
 #and if it is His-His or His-Met ligated.

 set i {0}
 foreach HResID $HEC {

   #--------------------------------------------------
   #Shift the residue ID of each heme
   set ShiftHResID [expr {$HResID + 1000}]
   set heme [atomselect top "resname HEC HEM and resid $HResID"]
   $heme set resid $ShiftHResID
   set NewHResID [expr {$HResID + ($i * 2)}]
   #--------------------------------------------------

   #--------------------------------------------------
   #Identify His residues nearby the selected heme
   set NumHIS {0}  
   set DistMinAx $DistMin
   set DistMaxAx $DistMax
   while {$NumHIS < 2 && $DistMinAx <= $DistMax} {
     set HIS [lsort -integer [[atomselect top "(resname HIS HIE HID HIP and name NE2 and within $DistMinAx of resname HEC HEM and resid $ShiftHResID) or (resname HIS HIE HID HIP and name N and within $DistMinAx of resname HEC HEM and resid $ShiftHResID)"] get resid]]
     set NumHIS [llength $HIS]
     set DistMinAx [expr {$DistMinAx+0.1}]
   }
   #--------------------------------------------------

   #--------------------------------------------------
   #If only one His residue was found, check for a nearby Met
   set NumMET {0}  
   set DistMinAx $DistMin
   set DistMaxAx $DistMax
   while {$NumMET < 1 && $DistMinAx <= $DistMax && $NumHIS == 1} {
     set MET [lsort -integer [[atomselect top "(resname MET and name SD and within $DistMinAx of resname HEC HEM and resid $ShiftHResID) or (resname MET and name N and within $DistMinAx of resname HEC HEM and resid $ShiftHResID)"] get resid]]
     set NumMET [llength $MET]
     set DistMinAx [expr {$DistMinAx+0.1}]
   }
   #--------------------------------------------------

   #--------------------------------------------------
   #Make sure the numbers of His and Met residues make sense
   if { $NumHIS == 0 && $NumMET == 0 } {
     puts " Problem: No His or Met ligands were found between $DistMin and $DistMax angstroms of heme-$NewHResID. "

   } elseif { $NumHIS == 1 && $NumMET == 0 } {
     puts " Problem: Only $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$NewHResID."
     puts "   There must be either two His or one His and one Met for the strucutre to be analyzed with BioDC. "

   } elseif { $NumHIS == 2 && $NumMET == 0 } {
     set HISresult [format " HIS: %5s %.2f %2d %10s " $ShiftHResID $DistMinAx $NumHIS $HIS]
     puts $HISresult
     puts " Heme-$NewHResID appears to be bis-histidine-ligated."

   } elseif { $NumHIS >= 2 && $NumMET == 0 } {
     puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$NewHResID."
     puts "   This may be a bis-histidine-ligated heme, but there are too many His residues within the selected distance range."
     puts "   It is unclear which of these His residues are actually coordinated. Correct this issue before continuing with BioDC. "

   } elseif { $NumHIS == 1 && $NumMET == 1 } {
     set HISresult [format " His: %5s %.2f %2d %10s " $ShiftHResID $DistMinAx $NumHIS $HIS]
     set METresult [format " MET: %5s %.2f %2d %10s " $ShiftHResID $DistMinAx $NumMET $MET]
     puts $HISresult
     puts $METresult
     puts " Heme-$NewHResID appears to be His-Met ligated."

   } elseif { $NumHIS == 2 && $NumMET == 1 } {
     puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$NewHResID."
     puts "   It is unclear if this heme is His-His or His-Met ligated. Setup with BioDc is not possible with this ambiguity."

   } elseif { $NumHIS >= 2 && $NumMET == 1 } {
     puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$NewHResID."
     puts "   It is unclear if this heme is His-His or His-Met ligated. Setup with BioDc is not possible with this ambiguity."
   }
   #--------------------------------------------------

   #--------------------------------------------------
   #Identify Cys residues nearby the selected heme
   set NumCYS {0}  
   set DistMinCYS $DistMin
   set DistMaxCYS $DistMax
   while {$NumCYS < 2 && $DistMinCYS <= $DistMaxCYS} {
    set CYS [lsort -integer [[atomselect top "resname CYS and name SG and within $DistMinCYS of resname HEC HEM and resid $ShiftHResID"] get resid]]
    set NumCYS [llength $CYS]
    set DistMinCYS [expr {$DistMinCYS+0.1}]
   }

   #--------------------------------------------------
   if { $NumCYS == 0 } {
     puts " No Cys residues found within $DistMinCYS - $DistMaxCYS angstroms of heme-$NewHResID"
     puts " This heme will be treated as a b-type heme. "

     #--------------------------------------------------
     if { $NumHIS == 2 && $NumMET == 0 } {
       if {[lindex $HIS 0] > [lindex $HIS 1]} {
         set HISp [lindex $HIS 1]
         set HISd [lindex $HIS 0]
       } elseif {[lindex $HIS 0] < [lindex $HIS 1]} {
         set HISp [lindex $HIS 0]
         set HISd [lindex $HIS 1]    
       }

       puts       "$HISp $HISd $ShiftHResID $NewHResID b HH"
       puts $out1 "$HISp $HISd $ShiftHResID $NewHResID b HH"
       puts $out2 "$HISp $HISd $ShiftHResID $NewHResID b HH"
     #--------------------------------------------------

     #--------------------------------------------------
     } elseif { $NumHIS == 1 && $NumMET == 1 } {
       set HISp [lindex $HIS 0]
       puts       "$HISp $MET $ShiftHResID $NewHResID b HM"
       puts $out1 "$HISp $MET $ShiftHResID $NewHResID b HM"
       puts $out3 "$HISp $MET $ShiftHResID $NewHResID b HM"
     #--------------------------------------------------

     #--------------------------------------------------
     } else {
       puts " Problem: I do not know what type of axial ligation is present."
       puts       "UNKNOWN LIGATION: $ShiftHResID $NewHResID "
       puts $out1 "UNKNOWN LIGATION! $ShiftHResID $NewHResID "
     }
     #--------------------------------------------------

   #--------------------------------------------------
   } elseif { $NumCYS == 1 } {
     puts " Problem: Only one Cys residue found within $DistMinCYS - $DistMaxCYS angstroms of heme-$NewHResID"
     puts " Cannot identify the heme as b- or c-type. This problem will cause the setup to fail. \n"
   #--------------------------------------------------

   #--------------------------------------------------
   } elseif { $NumCYS == 2 } {

     #--------------------------------------------------
     set CYSresult [format " CYS: %5s %.2f %2d %10s " $ShiftHResID $DistMinAx $NumCYS $CYS]
     puts $CYSresult
     puts " Two Cys residues found within $DistMinCYS - $DistMaxCYS angstroms of heme-$NewHResID"
     puts " This heme will therefore be treated as a c-type heme."

     set CYSb [lindex $CYS 0]
     set CYSc [lindex $CYS 1]
     set HISp [expr {$CYSc + 1}]

     #--------------------------------------------------
     if { $NumHIS == 2 && $NumMET == 0 } {
       if {[lindex $HIS 0] != $HISp} {
         set HISd [lindex $HIS 0]
       } else {
         set HISd [lindex $HIS 1]    
       }

       puts       "$CYSb $CYSc $HISp $HISd $ShiftHResID $NewHResID c HH"
       puts $out1 "$CYSb $CYSc $HISp $HISd $ShiftHResID $NewHResID c HH"
       puts $out2 "$CYSb $CYSc $HISp $HISd $ShiftHResID $NewHResID c HH"
     }
     #--------------------------------------------------

     #--------------------------------------------------
     if { $NumHIS == 1 && $NumMET == 1 } {
       puts       "$CYSb $CYSc $HISp $MET $ShiftHResID $NewHResID c HM"
       puts $out1 "$CYSb $CYSc $HISp $MET $ShiftHResID $NewHResID c HM"
       puts $out3 "$CYSb $CYSc $HISp $MET $ShiftHResID $NewHResID c HM"
     }
   }
 puts "=============================================================="
 set i [expr {$i+1}]
 }      

 close $out1
 close $out2
 close $out3

 set all [atomselect top all] """, file=open('CreateResIndexing.tcl', 'a'))

                print(f"""
 $all writepdb {PDBbasename}_renumd.pdb

 exit """, file=open('CreateResIndexing.tcl', 'a'))
        
                subprocess.run("vmd -e CreateResIndexing.tcl > CreateResIndexing.log", shell=True)
                return PDBbasename
                break
            elif (CreateResIndexingMethod == "MANUAL") or (CreateResIndexingMethod == "Manual") or (CreateResIndexingMethod == "manual") or (CreateResIndexingMethod == "MAN") or (CreateResIndexingMethod == "Man") or (CreateResIndexingMethod == "man") or (CreateResIndexingMethod == "M") or (CreateResIndexingMethod == "m"): 
                print("""
 To create ResIndexing.txt by hand:
    > Create a txt file with an editor of your choosing (e.g. 
      vi ResIndexing.txt). 

    > In this file, there must be one line for each heme cofactor 
      in your structure. 

    > An entry for a c-type heme has eight space-separated fields from left-to-right:
        1) Residue ID of Cys attached to the B-ring of the heme macrocycle 
        2) Residue ID of Cys attached to the C-ring of the heme macrocycle 
        3) Residue ID of proximal His ligated to the Fe center of the heme 
        4) Residue ID of distal His or Met ligand to the Fe center of the heme 
        5) Residue ID of the heme, shifted by +1000 (see below explaination)
        6) Residue ID that the heme will have in the fully prepared structure
        7) A "c" or "b" to indicate "c-type" or "b-type" heme
        8) A "HH" or "HM" to indicate His-His or His-Met ligation

    > An entry for a b-type heme has six space-separated fields from left-to-right:
        1) Residue ID of proximal His ligated to the Fe center of the heme 
        2) Residue ID of distal His or Met ligand to the Fe center of the heme 
        3) Residue ID of the heme, shifted by +1000 (see below explaination)
        4) Residue ID that the heme will have in the fully prepared structure
        5) A "c" or "b" to indicate "c-type" or "b-type" heme
        6) A "HH" or "HM" to indicate His-His or His-Met ligation

      A shift in the residue ID of the heme is needed to (potentially) 
      avoid assigning the propionates, whcih are made into separate 
      residues, to whatever came immediately after the heme in the 
      original PDB. (See the below example.)

      The residue ID the will be assigned in to the heme int he 
      properly formatted PDB for use with TLEaP is the original 
      residue ID for the heme + (i * 2), where i is a zero-based 
      index that counts the number of hemes in your system.

      To given an example:
        Let's say you have a di-heme system with bis-histidine-ligated
        c-type hemes. The residue IDs of the hemes are 114 and 115. 
        The entries in ResIndexing.txt should be:

        CysB1 CysC1 HisP1 HisD1 1114 114 c HH
        CysB2 CysC2 HisP2 HisD2 1115 117 c HH

        where CysB, CysC, HisP, and HisD are respectively defined by 
        points 1-4 above, and the 1 and 2 are used to distinguish the 
        Cys/His for the first and second hemes, respectively.

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
