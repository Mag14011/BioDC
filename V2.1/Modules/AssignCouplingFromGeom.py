################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def AssignCouplingFromGeom(LaunchDir, OutPrefix, InputDict):

    if (os.path.exists(f"{LaunchDir}/SPR")):
        StrucDir = f"{LaunchDir}/SPR"
    else:
        StrucDir = f"{LaunchDir}"

    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEM = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEM[idx] = int(line.strip().split(" ")[1])
            idx += 1

    ang = [0]*(len(HEM)-1)
    Hda = [0]*(len(HEM)-1)
    for idx in range(len(HEM)-1):
        print(f"""
 parm    {StrucDir}/{OutPrefix}.prmtop
 trajin  {StrucDir}/min.rst7
 
 vector Hem{HEM[idx]} corrplane :{HEM[idx]}&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 vector Hem{HEM[idx+1]} corrplane :{HEM[idx+1]}&@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA
 run
 
 vectormath vec1 Hem{HEM[idx]} vec2 Hem{HEM[idx+1]} dotangle out hemplaneorient_{HEM[idx]}-{HEM[idx+1]}.dat
 
 run
 quit""", file=open("CalcOrientation.in", "w"))

        if (os.path.isfile(f"hemplaneorient_{HEM[idx]}-{HEM[idx+1]}.dat")):
            print(f""" 
 Found hemplaneorient_{HEM[idx]}-{HEM[idx+1]}.dat from a prior execution.
 This prior output will be used for the analysis.""")

        if (not os.path.isfile(f"hemplaneorient_{HEM[idx]}-{HEM[idx+1]}.dat")):
            subprocess.run("cpptraj -i CalcOrientation.in > CalcOrientation.log", shell=True)

        lc = 0
        with open(f"hemplaneorient_{HEM[idx]}-{HEM[idx+1]}.dat") as fp:
            Lines = fp.readlines()
            for line in Lines:
                if (lc == 1):
                    ang[idx] = float(line.strip().split()[1])
                    if (ang[idx] > 90):
                        ang[idx] = 180 - ang[idx]
                        if (ang[idx] < 45):
                            Hda[idx] = 8.0
                        elif (ang[idx] >= 45):
                            Hda[idx] = 2.0
                    if (ang[idx] <= 90):
                        if (ang[idx] < 45):
                            Hda[idx] = 8.0
                        elif (ang[idx] >= 45):
                            Hda[idx] = 2.0
                    else:
                        Hda[idx] = "?"
                lc += 1

    print("\n Assigning coupling values based on inter-macrocycle planar angle ... ")
    for idx in range(len(HEM)-1):

        if (idx == 0):
            print("Hda(HEM-%0d <-> HEM-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" % (HEM[idx], HEM[idx+1], ang[idx], Hda[idx]), file=open('Hda.txt', 'w'))
        else:
            print("Hda(HEM-%0d <-> HEM-%0d) ang. = %10.3f deg.; Hda = %6.3f meV" % (HEM[idx], HEM[idx+1], ang[idx], Hda[idx]), file=open('Hda.txt', 'a'))

    print(" Done!")

    return Hda

################################################################################################################################################