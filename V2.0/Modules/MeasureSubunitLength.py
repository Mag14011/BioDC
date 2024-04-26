################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def MeasureSubunitLength():

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEM = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEM[idx] = int(line.strip().split(" ")[1])
            idx += 1

    print("""
 mol new min.pdb
 set output [open  "SubunitLength.txt" w]

 set FirstHeme [[atomselect top "resname HEM HEC HEH HCO HCR HBO HBR MCO MCR MBO MBR and resid %0d and name FE"] get index]
 set LastHeme  [[atomselect top "resname HEM HEC HEH HCO HCR HBO HBR MCO MCR MBO MBR and resid %0d and name FE"] get index]

 set dista  [measure bond [list $FirstHeme $LastHeme]]
 set distcm [expr {$dista * 1E-8}]

 puts $output "Proposed Subunit Length (cm) = $distcm"
 exit
    """ %(HEM[0], HEM[x-1]), file=open("MeasureSubunitLength.tcl", "w"))

    subprocess.run("vmd -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)
    #subprocess.run("/Applications/VMD\ 1.9.4a51-x86_64-Rev9.app/Contents/vmd/vmd_MACOSXX86_64 -e MeasureSubunitLength.tcl > MeasureSubunitLength.log", shell=True)

    with open("SubunitLength.txt") as fp:
        SubunitLength = float(fp.readline().strip().split()[5])

    return SubunitLength

################################################################################################################################################
