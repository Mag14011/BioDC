################################################################################################################################################
# Generic Modules
import math
from pandas import read_csv
################################################################################################################################################

################################################################################################################################################

def ComputeMarcusRates(Hda, Lambda, DG):

    NumSteps = len(DG)
    prefactor = [0]*(NumSteps)
    Eactf = [0]*(NumSteps)
    ketf = [0]*(NumSteps)
    Eactb = [0]*(NumSteps)
    ketb = [0]*(NumSteps)
    for idx in range(NumSteps):
        T = 300.0
        pi = 3.141592654
        kB = 8.6173304E-5
        hbar = 6.582119514E-16

        prefactor[idx] = (2 * pi * (Hda[idx]/1000)**2) / (hbar * math.sqrt(4 * pi * Lambda[idx] * kB* T))

        Eactf[idx] = ((DG[idx] + Lambda[idx])**2)/(4 * Lambda[idx])
        ketf[idx] = prefactor[idx] * math.exp((-1 * Eactf[idx])/(kB * T))

        Eactb[idx] = (((-1 * DG[idx]) + Lambda[idx])**2)/(4 * Lambda[idx])
        ketb[idx] = prefactor[idx] * math.exp((-1 * Eactb[idx])/(kB * T))

    for idx in range(NumSteps):
        print("""
 Step #%0d:
    Activation Energy:
        Forward: %.3E eV
        Reverse: %.3E eV
    Rates:
        Forward: %.3E 
        Reverse: %.3E\n""" %(idx, Eactf[idx], Eactb[idx], ketf[idx], ketb[idx]), end=" ")

        if (idx == 0):
            print("ketf,ketb", file=open("rates.txt", "w"))
            print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "a"))
        else:
            print("%.3E,%3E" %(ketf[idx], ketb[idx]), file=open("rates.txt", "a"))

################################################################################################################################################
