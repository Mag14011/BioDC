import subprocess

def generate_tcl_script(DistMin, DistMax, pdb_filename):
    tcl_script = f"""
#--------------------------------------------------
set DistMin {DistMin}
set DistMax {DistMax}
mol new {pdb_filename}
#--------------------------------------------------

#--------------------------------------------------
#Create output files
set out1 [open "ResIndexing.txt" w]
set out2 [open "ResIndexing_HisHisHeme.txt" w]
set out3 [open "ResIndexing_HisMetHeme.txt" w]
#--------------------------------------------------

#--------------------------------------------------
#Identify the heme residue IDs
set HEM [[atomselect top "name FE M1 M2 M3 M4 M5 M6 M7 M8"] get resid]
#--------------------------------------------------

#--------------------------------------------------
#Determine whether each heme is of the b- or c-type,
#and if it is His-His or His-Met ligated.

set i 0
foreach HResID $HEM {{

    #--------------------------------------------------
    #Identify His residues nearby the selected heme
    set NumHIS 0  
    set DistMinAx $DistMin
    set DistMaxAx $DistMax

    while {{$NumHIS < 2 && $DistMinAx <= $DistMax}} {{

        set HIS [lsort -integer [[atomselect top "
            (resname PHO and name NE2 S1 and within $DistMinAx of resid $HResID and name FE M7) or 
            (resname DHO and name NE2 S2 and within $DistMinAx of resid $HResID and name FE M7) or 
            (resname PHR and name NE2 T1 and within $DistMinAx of resid $HResID and name FE M8) or 
            (resname DHR and name NE2 T2 and within $DistMinAx of resid $HResID and name FE M8) or 
            (resname FHO and name NE2 Y1 and within $DistMinAx of resid $HResID and name FE M1) or 
            (resname RHO and name NE2 Y2 and within $DistMinAx of resid $HResID and name FE M1) or 
            (resname FHR and name NE2 Z1 and within $DistMinAx of resid $HResID and name FE M2) or 
            (resname RHR and name NE2 Z2 and within $DistMinAx of resid $HResID and name FE M2) or 
            (resname PMO and name NE2 U1 and within $DistMinAx of resid $HResID and name FE M5) or 
            (resname PMR and name NE2 V1 and within $DistMinAx of resid $HResID and name FE M6) or 
            (resname FMO and name NE2 W2 and within $DistMinAx of resid $HResID and name FE M3) or 
            (resname FMR and name NE2 X2 and within $DistMinAx of resid $HResID and name FE M4)"] get resid]]

        set NumHIS [llength $HIS]
        set DistMinAx [expr {{$DistMinAx+0.1}}]
    }}
    #--------------------------------------------------

    #--------------------------------------------------
    #If only one His residue was found, check for a nearby Met
    set NumMET 0  
    set DistMinAx $DistMin
    set DistMaxAx $DistMax

    while {{$NumMET < 1 && $DistMinAx <= $DistMax && $NumHIS == 1}} {{

        set MET [lsort -integer [[atomselect top "
            (resname DMO and name SD U2 and within $DistMinAx of resid $HResID and name FE M5) or 
            (resname DMR and name SD V2 and within $DistMinAx of resid $HResID and name FE M6) or 
            (resname RMO and name SD W1 and within $DistMinAx of resid $HResID and name FE M3) or 
            (resname RMR and name SD X1 and within $DistMinAx of resid $HResID and name FE M4)"] get resid]]

        set NumMET [llength $MET]
        set DistMinAx [expr {{$DistMinAx+0.1}}]
    }}
    #--------------------------------------------------

    #--------------------------------------------------
    #Make sure the numbers of His and Met residues make sense
    if {{ $NumHIS == 0 && $NumMET == 0 }} {{
        puts " Problem: No His or Met ligands were found between $DistMin and $DistMax angstroms of heme-$HResID. "

    }} elseif {{ $NumHIS == 1 && $NumMET == 0 }} {{
        puts " Problem: Only $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$HResID."
        puts "   There must be either two His or one His and one Met for the structure to be analyzed with BioDC. "

    }} elseif {{ $NumHIS == 2 && $NumMET == 0 }} {{
        set HISresult [format " HIS: %5s %.2f %2d %10s " $HResID $DistMinAx $NumHIS $HIS]
        puts $HISresult
        puts " Heme-$HResID appears to be bis-histidine-ligated."

    }} elseif {{ $NumHIS >= 2 && $NumMET == 0 }} {{
        puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$HResID."
        puts "   This may be a bis-histidine-ligated heme, but there are too many His residues within the selected distance range."
        puts "   It is unclear which of these His residues are actually coordinated. Correct this issue before continuing with BioDC. "

    }} elseif {{ $NumHIS == 1 && $NumMET == 1 }} {{
        set HISresult [format " His: %5s %.2f %2d %10s " $HResID $DistMinAx $NumHIS $HIS]
        set METresult [format " MET: %5s %.2f %2d %10s " $HResID $DistMinAx $NumMET $MET]
        puts $HISresult
        puts $METresult
        puts " Heme-$HResID appears to be His-Met ligated."

    }} elseif {{ $NumHIS == 2 && $NumMET == 1 }} {{
        puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$HResID."
        puts "   It is unclear if this heme is His-His or His-Met ligated. Setup with BioDC is not possible with this ambiguity."

    }} elseif {{ $NumHIS >= 2 && $NumMET == 1 }} {{
        puts " Problem: $NumHIS His and $NumMET Met residues were found between $DistMin and $DistMax angstroms of heme-$HResID."
        puts "   It is unclear if this heme is His-His or His-Met ligated. Setup with BioDC is not possible with this ambiguity."
    }}
    #--------------------------------------------------

    #--------------------------------------------------
    #Identify Cys residues nearby the selected heme
    set NumCYS 0  
    set DistMinCYS $DistMin
    set DistMaxCYS $DistMax

    while {{$NumCYS < 2 && $DistMinCYS <= $DistMaxCYS}} {{
        set CYS [lsort -integer [[atomselect top "
        (resname CYO and name CA and within $DistMinCYS of resid $HResID and name CBB2) or 
        (resname CYO and name CA and within $DistMinCYS of resid $HResID and name CBC1)"] get resid]]

        set NumCYS [llength $CYS]
        set DistMinCYS [expr {{$DistMinCYS+0.1}}]
    }}

    #--------------------------------------------------
    if {{ $NumCYS == 0 }} {{
        puts " No Cys residues found within $DistMinCYS - $DistMaxCYS angstroms of heme-$HResID"
        puts " This heme will be treated as a b-type heme. "

        #--------------------------------------------------
        if {{ $NumHIS == 2 && $NumMET == 0 }} {{
            if {{[lindex $HIS 0] > [lindex $HIS 1]}} {{
                set HISp [lindex $HIS 1]
                set HISd [lindex $HIS 0]
            }} elseif {{[lindex $HIS 0] < [lindex $HIS 1]}} {{
                set HISp [lindex $HIS 0]
                set HISd [lindex $HIS 1]    
            }}

            puts       "$HISp $HISd $HResID $HResID b HH"
            puts $out1 "$HISp $HISd $HResID $HResID b HH"
            puts $out2 "$HISp $HISd $HResID $HResID b HH"
            #--------------------------------------------------

            #--------------------------------------------------
        }} elseif {{ $NumHIS == 1 && $NumMET == 1 }} {{
            set HISp [lindex $HIS 0]
            puts       "$HISp $MET $HResID $HResID b HM"
            puts $out1 "$HISp $MET $HResID $HResID b HM"
            puts $out3 "$HISp $MET $HResID $HResID b HM"
            #--------------------------------------------------

            #--------------------------------------------------
        }} else {{
            puts " Problem: I do not know what type of axial ligation is present."
            puts       "UNKNOWN LIGATION: $HResID $HResID "
            puts $out1 "UNKNOWN LIGATION! $HResID $HResID "
        }}
        #--------------------------------------------------

        #--------------------------------------------------
    }} elseif {{ $NumCYS == 1 }} {{
        puts " Problem: Only one Cys residue found within $DistMinCYS - $DistMaxCYS angstroms of heme-$HResID"
        puts " Cannot identify the heme as b- or c-type. This problem will cause the setup to fail."
    #--------------------------------------------------

    #--------------------------------------------------
    }} elseif {{ $NumCYS == 2 }} {{

        #--------------------------------------------------
        set CYSresult [format " CYS: %5s %.2f %2d %10s " $HResID $DistMinAx $NumCYS $CYS]
        puts $CYSresult
        puts " Two Cys residues found within $DistMinCYS - $DistMaxCYS angstroms of heme-$HResID"
        puts " This heme will therefore be treated as a c-type heme."

        set CYSb [lindex $CYS 0]
        set CYSc [lindex $CYS 1]
        set HISp [expr {{$CYSc + 1}}]

        #--------------------------------------------------
        if {{ $NumHIS == 2 && $NumMET == 0 }} {{
            if {{[lindex $HIS 0] != $HISp}} {{
                set HISd [lindex $HIS 0]
            }} else {{
                set HISd [lindex $HIS 1]    
            }}

            puts       "$CYSb $CYSc $HISp $HISd $HResID $HResID c HH"
            puts $out1 "$CYSb $CYSc $HISp $HISd $HResID $HResID c HH"
            puts $out2 "$CYSb $CYSc $HISp $HISd $HResID $HResID c HH"
        }}
        #--------------------------------------------------

        #--------------------------------------------------
        if {{ $NumHIS == 1 && $NumMET == 1 }} {{
            puts       "$CYSb $CYSc $HISp $MET $HResID $HResID c HM"
            puts $out1 "$CYSb $CYSc $HISp $MET $HResID $HResID c HM"
            puts $out3 "$CYSb $CYSc $HISp $MET $HResID $HResID c HM"
        }}
    }}
    puts "=============================================================="
    set i [expr {{$i+1}}]
}}      

close $out1
close $out2
close $out3

exit 
"""

    with open('RecreateResIndexing.tcl', 'w') as file:
        file.write(tcl_script)

    return 'RecreateResIndexing.tcl'

def submit_to_vmd(DistMin, DistMax, pdb_filename, log_filename="vmd_output.log"):
    tcl_script_path = generate_tcl_script(DistMin, DistMax, pdb_filename)

    # Run VMD with the generated TCL script and redirect output to a log file
    try:
        with open(log_filename, 'w') as log_file:
            process = subprocess.Popen(
                ["vmd", "-dispdev", "text", "-e", tcl_script_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True
            )

            # Read the output line by line and write it to the log file
            for line in process.stdout:
                log_file.write(line)
                log_file.flush()  # Ensure the output is written immediately

            # Wait for the process to complete
            return_code = process.wait()

            if return_code == 0:
                print(f"VMD script executed successfully. Output saved to {log_filename}")
            else:
                print(f"VMD script encountered an error. Check {log_filename} for details.")

    except Exception as e:
        print(f"Error running VMD: {e}")
        # Handle the error appropriately

# Example usage
if __name__ == "__main__":
    DistMin = 2.0
    DistMax = 4.0
    pdb_filename = "min.pdb"
    log_filename = "RecreateResIndexing.log"
    submit_to_vmd(DistMin, DistMax, pdb_filename, log_filename)
