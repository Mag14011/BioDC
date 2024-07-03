################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
################################################################################################################################################

################################################################################################################################################

def StructRelax(StrucDir, OutPrefix, SolvEnv, InputDict):

    if (os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop") == True and os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7")) or (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop") == True and os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.rst7")) or (os.path.isfile(f"{StrucDir}/{OutPrefix}.prmtop") == True and os.path.isfile(f"{StrucDir}/{OutPrefix}.rst7")):
        if (os.path.isfile(f"{StrucDir}/{OutPrefix}_new.prmtop")):
            pass
            print(f" Found {StrucDir}/{OutPrefix}_new.prmtop and {StrucDir}/{OutPrefix}_reord.rst7")
        elif (os.path.isfile(f"{StrucDir}/{OutPrefix}_reord.prmtop")):
            pass
            print(f" Found {StrucDir}/{OutPrefix}_reord.prmtop and {StrucDir}/{OutPrefix}_reord.rst7") 
        else:
            pass
            print(" Found {StrucDir}/{OutPrefix}.prmtop and {StrucDir}/{OutPrefix}.rst7")

        print(" Preparing to relax the geometry")

        cwd = os.getcwd()
        os.chdir(f"{StrucDir}")

        if (SolvEnv == "explicit") or (SolvEnv == "Explicit") or (SolvEnv == "e") or (SolvEnv == "E"):
            print("""
 Energy Minimization Stage in Implicit Solvent
 &cntrl
  imin=1,            ! Perform an energy minimization
  ntb=1,             ! Constant volumne
  cut=10.0,          ! Non-bonded cutoff in angstroms
  ntmin=1,           ! Steepest descent + conjugate gradient method
  ncyc=1000,         ! Number of steepest descent cycles
  maxcyc=5000,       ! Maximum number of minimization cycles
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
 /
            """, file=open('min.in', 'w'))

            while True:
                if ("CompChoice" in InputDict):
                    CompChoice = InputDict["CompChoice"]
                    print(f"CompChoice = {CompChoice}", file=open("InteractiveInput.txt", 'a'))
                else:
                    CompChoice = input("\n Run the minimization using SANDER (S) or PMEMD (P)? ")
                    print(f"CompChoice = {CompChoice}", file=open("InteractiveInput.txt", 'a'))

                if (CompChoice == "SANDER") or (CompChoice == "Sander") or (CompChoice == "sander") or (CompChoice == "S") or (CompChoice == "s"):
                    print(" Running minimization ...")
                    if (os.path.isfile(OutPrefix+"_new.prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}_new.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+"_reord.prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}_reord.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+".prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}.rst7", shell=True)
                    print(" Minimization finished!")
                    break
                elif (CompChoice == "PMEMD") or (CompChoice == "Pmemd") or (CompChoice == "pmemd") or (CompChoice == "P") or (CompChoice == "p"): 
                    while True:
                        if ("NProc" in InputDict):
                            NProc = InputDict["NProc"]
                            print(NProc)
                            print(f"NProc = {NProc}", file=open("InteractiveInput.txt", 'a'))
                            break
                        else:
                            NProc = input(" Parallelize the minimization over how many CPUs? ")
                            print(f"NProc = {NProc}", file=open("InteractiveInput.txt", 'a'))

                    print(" Running minimization ...")
                    if (os.path.isfile(OutPrefix+"_new.prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}_new.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+"_reord.prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}_reord.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+".prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}.rst7", shell=True)
                    print(" Minimization finished!")
                    break
                else:
                    print(" Sorry, I didn't understand your response. Please try again.")
            os.chdir(cwd)
        else:
            pass

        if (SolvEnv == "IMPLICIT") or (SolvEnv == "Implicit") or (SolvEnv == "implicit") or (SolvEnv == "I") or (SolvEnv == "i"):
            print("""
Energy Minimization in Implicit Solvent
&cntrl
  imin=1,            ! Perform an energy minimization
  ntb=0,             ! Non-periodic
  cut=9999           ! Non-bonded cutoff in Å
  ntmin=1,           ! steepest descent + conjugate gradient method 
  ncyc=200,          ! Number of steepest descent cycles
  maxcyc=500,        ! Maximum number of minimization cycles
  igb=2,             ! Generalized Born implicit solvent model
  saltcon=0.1,       ! salt concentration in M
  ntwr=100,          ! Restart file written every ntwr steps
  ntwx=100,          ! Trajectory file written every ntwx steps
  ntpr=100,          ! The mdout and mdinfo files written every ntpr steps
  ntr=1,             ! Turn on positional restraints
  restraintmask='@CA,C,O,N&!:WAT|@FE,NA,NB,NC,ND,C3D,C2A,C3B,C2C,CA,CB',
  restraint_wt=10.0, ! 10 kcal/mol.A**2 restraint force constant
/
            """, file=open('min.in', 'w'))

            while True:
                CompChoice = input("\n Run the minimization using SANDER (S) or PMEMD (P)? ")

                if (CompChoice == "SANDER") or (CompChoice == "Sander") or (CompChoice == "sander") or (CompChoice == "S") or (CompChoice == "s"):
                    print(" Running minimization ...")
                    if (os.path.isfile(OutPrefix+"_new.prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}_new.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+"_reord.prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}_reord.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+".prmtop")):
                        subprocess.run(f"sander -O -i min.in -o min.out -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}.rst7", shell=True)
                    print(" Minimization finished!")
                    break
                elif (CompChoice == "PMEMD") or (CompChoice == "Pmemd") or (CompChoice == "pmemd") or (CompChoice == "P") or (CompChoice == "p"): 
                    NProc = input(" parallelize the minimization over how many CPUs? ")
                   #NProc = 2

                    print(" Running minimization ...")
                    if (os.path.isfile(OutPrefix+"_new.prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}_new.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+"_reord.prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}_reord.prmtop -c {OutPrefix}_reord.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}_reord.rst7", shell=True)
                    elif (os.path.isfile(OutPrefix+".prmtop")):
                        subprocess.run(f"mpirun -np {NProc} pmemd.MPI -O -i min.in -o min.out -p {OutPrefix}.prmtop -c {OutPrefix}.rst7 -inf min.mdinfo -r min.rst7 -ref {OutPrefix}.rst7", shell=True)
                    print(" Minimization finished!")
                    break
                else:
                    print(" Sorry, I didn't understand your response. Please try again.")
    else:
        print(f" {StrucDir}/{OutPrefix}_new.prmtop/{StrucDir}/{OutPrefix}_reord.rst7 or {StrucDir}/{OutPrefix}.prmtop/{StrucDir}/{OutPrefix}.rst7 were not found", end=" ")
        sys.exit(" Nothing to minimize. Something went wrong in the preceeding steps!")

    os.chdir(cwd)

################################################################################################################################################
