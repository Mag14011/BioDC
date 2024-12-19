import os
import sys
import subprocess
from subprocess import Popen

def LambdaFromSASA(OutPrefix, InputDict):

    if not os.path.isfile("min.pdb"):
        sys.exit("""
 The minimized structure (min.pdb) is missing.
 Something went wrong in a prior step and 
 we cannot proceed. I apologize for the 
 inconvenience!""")

    if os.path.isfile("min.pdb"):
        idx = 0
        if os.path.isfile("SelResIndexing.txt"):
            with open("SelResIndexing.txt") as fp:
                NumHEC = len(fp.readlines())
                HEM = [0] * NumHEC
                ActiveIDs = [0] * NumHEC
                fp.seek(0)
 
                Lines = fp.readlines()
                for line in Lines:
                    EntryLength = len(line.strip().split(" "))
                    HemeType = line.strip().split(" ")[-2]
                    AxLigType = line.strip().split(" ")[-1]

                    if EntryLength == 8 and HemeType == "c" and AxLigType == "HH":
                        CYSb = int(line.strip().split(" ")[0])
                        CYSc = int(line.strip().split(" ")[1])
                        HISp = int(line.strip().split(" ")[2])
                        HISd = int(line.strip().split(" ")[3])
                        HEM[idx] = int(line.strip().split(" ")[5])
                        #Includes propionic acids
                        #ActiveIDs[idx] = f"{HISp} {HISd} {HEM[idx]} {HEM[idx]+1} {HEM[idx]+2}"
                        #Excludes propionic acids
                        ActiveIDs[idx] = f"{HISp} {HISd} {HEM[idx]}"
                    elif EntryLength == 8 and HemeType == "c" and AxLigType == "HM":
                        CYSb = int(line.strip().split(" ")[0])
                        CYSc = int(line.strip().split(" ")[1])
                        HISp = int(line.strip().split(" ")[2])
                        METd = int(line.strip().split(" ")[3])
                        HEM[idx] = int(line.strip().split(" ")[5])
                        #Includes propionic acids
                        #ActiveIDs[idx] = f"{HISp} {METd} {HEM[idx]} {HEM[idx]+1} {HEM[idx]+2}"
                        #Excludes propionic acids
                        ActiveIDs[idx] = f"{HISp} {METd} {HEM[idx]}"
                    elif EntryLength == 6 and HemeType == "b" and AxLigType == "HH":
                        HISp = int(line.strip().split(" ")[0])
                        HISd = int(line.strip().split(" ")[1])
                        HEM[idx] = int(line.strip().split(" ")[3])
                        #Includes propionic acids
                        #ActiveIDs[idx] = f"{HISp} {HISd} {HEM[idx]} {HEM[idx]+1} {HEM[idx]+2}"
                        #Excludes propionic acids
                        ActiveIDs[idx] = f"{HISp} {HISd} {HEM[idx]}"
                    elif EntryLength == 6 and HemeType == "b" and AxLigType == "HM":
                        HISp = int(line.strip().split(" ")[0])
                        METd = int(line.strip().split(" ")[1])
                        HEM[idx] = int(line.strip().split(" ")[3])
                        #Includes propionic acids
                        #ActiveIDs[idx] = f"{HISp} {METd} {HEM[idx]} {HEM[idx]+1} {HEM[idx]+2}"
                        #Excludes propionic acids
                        ActiveIDs[idx] = f"{HISp} {METd} {HEM[idx]}"
                    else:
                        print(f" *** Missing entries on line number {idx+1} of SelResIndexing.txt!")

                    idx += 1
        else:
            sys.exit("""
 SelResIndexing.txt is missing.
 Something went wrong when you defined
 the linear sequence of hemes.

 This problem must be resolved before 
 proceeding.""")

        for idx in range(len(HEM) - 1):
            print(f"""
 mol new min.pdb
 set HEM1 {HEM[idx]}; set HEM2 {HEM[idx+1]}""", file=open('SASACalc.tcl', 'w'))  

            print("""
 set output [open "${HEM1},${HEM2}_SASAanalysis.dat" a]""", file=open('SASACalc.tcl', 'a'))

            print(f"""
 set allsel     [atomselect top "all and not resname WAT 'Na+' 'Cl-'"]

 set donor      [atomselect top "not water and not ions and resid {ActiveIDs[idx]} and not name N H CA HA C O"]
 set Dfe        [atomselect top "resid {HEM[idx]} and name FE"]
 set DfeIDX     [$Dfe get index]
 set DfeID      [$Dfe get resid]

 set acceptor   [atomselect top "not water and not ions and resid {ActiveIDs[idx+1]} and not name N H CA HA C O"]
 set Afe        [atomselect top "resid {HEM[idx+1]} and name FE"]
 set AfeIDX     [$Afe get index]
 set AfeID      [$Afe get resid]
 set DA         [list $DfeIDX $AfeIDX]""", file=open('SASACalc.tcl', 'a'))

            print("""
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

            if os.path.isfile(f"{HEM[idx]},{HEM[idx+1]}_SASAanalysis.dat"):
                print(f""" Found {HEM[idx]},{HEM[idx+1]}_SASAanalysis.dat from a prior execution.
 This prior output will be used for the analysis.
  """) 
            else:
                print(f" Now using VMD to compute SASA Donor = {HEM[idx]} & Acceptor = {HEM[idx+1]}...")
                subprocess.run("vmd -e SASACalc.tcl > SASACalc.log", shell=True)

        print(" Computing Reorganization Energy from Solvent Accessibility...")
        alpha = 5.18
        beta = 0.016
        Rd = 4.6
        Ra = 4.6
        Eopt = 1.84 

        Dsasa = [0] * (len(HEM) - 1)
        Asasa = [0] * (len(HEM) - 1)
        TotalSASA = [0] * (len(HEM) - 1)
        Es = [0] * (len(HEM) - 1)
        M = [0] * (len(HEM) - 1)
        Rda = [0] * (len(HEM) - 1)
        R = [0] * (len(HEM) - 1)
        Lambda = [0] * (len(HEM) - 1)
        for idx in range(len(HEM) - 1):
            with open(f"{HEM[idx]},{HEM[idx+1]}_SASAanalysis.dat") as fp:
                Lines = fp.readlines()
                for line in Lines:
                    Dsasa[idx] = float(line.strip().split(" ")[1])
                    Asasa[idx] = float(line.strip().split(" ")[2])
                    Rda[idx] = float(line.strip().split(" ")[3])

            TotalSASA[idx] = Dsasa[idx] + Asasa[idx]
            Es[idx] = alpha + (beta * TotalSASA[idx])
            M[idx] = (1 / Eopt) - (1 / Es[idx])
            R[idx] = (1 / ((2 * Rd) / 0.53)) + (1 / ((2 * Ra) / 0.53)) - (1 / (Rda[idx] / 0.53))
            Lambda[idx] = ((-1) ** 2) * (M[idx]) * (R[idx]) * (27.2114)

        for idx in range(len(HEM) - 1):
            if idx == 0:
                print(""" 
 HEH-%0d -> HEM-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" % (HEM[idx], HEM[idx + 1], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'w'))

            if idx != 0:
                print(""" 
 HEH-%0d -> HEM-%0d --------- 
 Dsasa     = %.3f
 Asasa     = %.3f
 Rda       = %.3f
 TotalSASA = %.3f
 Es        = %.3f
 ----------------------------
 Reorg. Eng. = %.3f""" % (HEM[idx], HEM[idx + 1], Dsasa[idx], Asasa[idx], Rda[idx], TotalSASA[idx], Es[idx], Lambda[idx]), file=open('Lambda.txt', 'a'))

    print(" Done!")

    return Lambda
