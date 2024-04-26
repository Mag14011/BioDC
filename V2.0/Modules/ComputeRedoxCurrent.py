################################################################################################################################################
# Generic Modules
import os
import sys
import math
import subprocess
from subprocess import Popen
################################################################################################################################################
import MeasureAvgDist
import ComputeDiffusionCoefficient
import MeasureSubunitLength
################################################################################################################################################

def ComputeRedoxCurrent(Javg):

    r=0.75E-7      #cm
    c=3.00E10      #cm/s
    e=1.602E-19    #c
    kb=1.38E-23    #J/k
    hbar=1.05E-34  #J/s 
    deltabar=3     #cm^-1

    print(" Please provide the following parmaeters: ")

    while True:
        try:
            T = float(input("  Temperature (K)? "))
        except ValueError:
            print(" Your entry must be a floating-point number.")
        else:
            break

    while True:
        try:
            cps = int(input("  Number of charges per subunit? "))
        except ValueError:
            print(" Your entry must be an integer.")
        else:
            break

    ahs  = (round(sum(MeasureAvgDist.MeasureAvgDist())/len(MeasureAvgDist.MeasureAvgDist()), 3))*1E-8
    print(f"  Average heme spacing is {ahs} cm")

    while True:
        try:
            Ncnt = float(input("  Fewest number of contacts at either protein-electrode interface? "))
        except ValueError:
            print(" Your entry must be an integer.")
        else:
            break
    
    DiffusionCoefficient = ComputeDiffusionCoefficient.ComputeDiffusionCoefficient(ahs)
    
    SubunitLength = MeasureSubunitLength.MeasureSubunitLength()

    print(""" 
 The length of a subunit of the cytochrome polymer is needed.
  The subunit length measured between the first and the last heme
  specified in LinearizedHemeSequence.txt is %.2E cm\n""" %(SubunitLength))

    lsub = SubunitLength

    while True:
        try:
            lw = float(input("  Length of wire (cm)? "))
        except ValueError:
            print(" Your entry must be a floating-point number.")
        else:
            break

    cpsul = ((cps)/(lsub))
    csa = (math.pi * (r)**2) 
    crgden = (cpsul)/(csa)

    while True:
        IsGKnown = input(f"""
 To interconvert between electrical conductivity and electron flux, the
 voltage must be specified.

 The multi-particle flux kinetic model of Blumberger and co-workers
 assumes here a hard-coded input/output electorn rate of 1E12 electrons/sec.
 If the expeirmental conductance is known, the correpsonding voltage to 
 supply this current can be calculated. If the experimental conductance is
 not known, the voltage must be chosen to characterize the electrical 
 properties.

 Is the conductance (in S) known (yes/no)? """)

        if (IsGKnown == "YES") or (IsGKnown == "Yes") or (IsGKnown == "yes") or (IsGKnown == "Y") or (IsGKnown == "y"):
            while True:
                try:
                    Gexp = float(input("""  Experimental Conductance (S)? """))
                except ValueError:
                    print(" Your entry must be a floating-point number.")
                else:
                    break

            V = ((1E12)*(e))/Gexp
            s_exp = Gexp*(lw/csa)
            j_exp = ((s_exp*csa)/(e*Ncnt*lw))*(V)
            u_exp = s_exp/((e)*(crgden))
            d_exp = ((Gexp * kb * T * lw) / (csa * crgden * (e)**2))
            k_exp = (d_exp)/((ahs)**2)

            break

        elif (IsGKnown == "NO") or (IsGKnown == "No") or (IsGKnown == "no") or (IsGKnown == "N") or (IsGKnown == "n"):
            while True:
                try:
                    print(""" 
 Because the experimental conductance is not known, 
 it will be set to zero in the output""")

                    Gexp = 0.0 
                    V = float(input("""  Properties should be computed at what voltage (V)? """))
                except ValueError:
                    print(" Your entry must be a floating-point number.")
                else:
                    break
            break

        else:
            print(" Sorry, I didn't understand your response.")

    if (os.path.isfile("D.txt") == True):
        with open("D.txt") as fp:
            Dcalc = float(fp.readline().strip().split()[3])

    k_dif = (Dcalc)/((ahs)**2)
    u_dif = (Dcalc)*(e/(kb*T))
    s_dif = e*(crgden)*(u_dif)
    j_dif = ((s_dif*csa)/(e*Ncnt*lw))*(V)

    s_flx = ((Javg*e*Ncnt*lw)/(V*csa)) 
    u_flx = s_flx/((e)*(crgden))
    d_flx = ((kb*T)/(e))*(u_flx)
    k_flx = (d_flx)/((ahs)**2)

    u_bt = ((e)*(ahs**2))/(2*hbar)
    d_bt = ((kb*T)/(e))*(u_bt)
    k_bt = (d_bt)/((ahs)**2)
    s_bt = (e)*(crgden)*(u_bt)
    j_bt = ((s_bt*csa)/(e*Ncnt*lw))*(V)

    u_hp = (((2)*(math.pi)*(c)*(e)*(ahs**2))/((kb)*(T)))*(deltabar)
    d_hp = ((kb*T)/(e))*(u_hp)
    k_hp = (d_hp)/((ahs)**2)
    s_hp = (e)*(crgden)*(u_hp)
    j_hp = ((s_hp*csa)/(e*Ncnt*lw))*(V)

    print("""
 Entered Quantities 
   Temperature                                       = %.1f K    
   Average Heme Spacing                              = %e cm
   Charge per Subunit Length                         = %e q/cm 
   Length of Filament                                = %e cm 
   Fewest number of protein-electrode contacts       = %e cm 
   Experimental Conductance                          = %s S

 Computed Quantities 
   1) Structure Properties:
        Cross-Sectional Area                         = %e cm^2
        Charge Density                               = %e q/cm^2

   2) Single Particle Diffusion Model
        Diffusion Constant                           = %e cm^2/s
        Homogenious Chain Hopping Rate               = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   3) Multi-Particle Steady-State Model
        Diffusion Constant                           = %e cm^2/s
        Homogenious Chain Hopping Rate               = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   4) Band Theory Minimum Requirements
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s

   5) Hopping Rate Maximum Limits
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at at %.1f V                            = %e electrons/s
    """ %(T, ahs, cps, lw, Ncnt, Gexp, csa, crgden, Dcalc, k_dif, u_dif, s_dif, V, j_dif, d_flx, k_flx, u_flx, s_flx, V, Javg, d_bt, k_bt, u_bt, s_bt, V, j_bt, d_hp, k_hp, u_hp, s_hp, V, j_hp), file=open("RedoxCurrentAnalysis.txt", 'w')) 

    if (IsGKnown == "YES") or (IsGKnown == "Yes") or (IsGKnown == "yes") or (IsGKnown == "Y") or (IsGKnown == "y"):
        print("""
   6) Experiment-Based Quantities
        Diffusion Constant                           = %e cm^2/s
        Hopping Rate                                 = %e s^-1
        Charge Mobility                              = %e cm^2/Vs
        Conductivity                                 = %e S/cm
        Flux at %.1f V                               = %e electrons/s
        """ %(d_exp, k_exp, u_exp, s_exp, V, j_exp), file=open("RedoxCurrentAnalysis.txt", 'a')) 

    print(" %6s %8s %8s" %("Voltage (V)", "Exp. Current (pA)", "Computed Current (pA)"), file=open("RedoxCurrentAnalysis.txt", 'a'))
    for V in [x / 10.0 for x in range(-5, 5, 1)]:
        Iexp = ((Gexp * V) * 1E12)
        Icmp = (((csa * crgden * (e)**2 * Dcalc) / (kb * T * lw)) * V) * 1E12 
        print(" %11.3f %17.3f %21.3E" %(V, Iexp, Icmp), file=open("RedoxCurrentAnalysis.txt", 'a'))

    print(" The following analysis has been saved to RedoxCurrentAnalysis.txt:") 
    with open("RedoxCurrentAnalysis.txt") as rca:
        print(rca.read())

################################################################################################################################################
