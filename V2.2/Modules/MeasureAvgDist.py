################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def write_r_txt(distances, LaunchDir):
    with open(f"{LaunchDir}/EE/r.txt", "w") as f:
        for idx, distance in enumerate(distances):
            f.write(f"r = {distance:.2f} Å\n")
#   print("r.txt file has been created successfully.")

def MeasureAvgDist(LaunchDir):
    idx = 0
    with open("LinearizedHemeSequence.txt") as fp:
        x = len(fp.readlines())
        HEM = [0]*x

        fp.seek(0)
        Lines = fp.readlines()
        for line in Lines:
            HEM[idx] = int(line.strip().split(" ")[1])
            idx += 1

    dist = [0]*(len(HEM)-1)
    for idx in range(len(HEM)-1):
        print(f"""
parm min.pdb
trajin min.pdb

nativecontacts :{HEM[idx]}&(@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA) :{HEM[idx+1]}&(@FE,NA,C1A,C2A,C3A,C4A,CHB,C1B,NB,C2B,C3B,C4B,CHC,C1C,NC,C2C,C3C,C4C,CHD,C1D,ND,C2D,C3D,C4D,CHA) out HEM{HEM[idx]}-{HEM[idx+1]}_dist.dat mindist
run
        """, file=open(f"MeasureDistance_{HEM[idx]}-{HEM[idx+1]}.in", "w"))

        subprocess.run(f"cpptraj -i MeasureDistance_{HEM[idx]}-{HEM[idx+1]}.in > MeasureDistance_{HEM[idx]}-{HEM[idx+1]}.log", shell=True)

        if (os.path.isfile(f"HEM{HEM[idx]}-{HEM[idx+1]}_dist.dat") == True):
            with open(f"HEM{HEM[idx]}-{HEM[idx+1]}_dist.dat",'r') as fp:
                Lines = fp.readlines()[1:]
                for line in Lines:
                    dist[idx] = float(line.strip().split()[3])

        idx += 1

    # Write the r.txt file
    write_r_txt(dist, LaunchDir)

    return dist

################################################################################################################################################
