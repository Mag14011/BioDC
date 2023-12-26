#!/usr/bin/env python
import hopping as hp
import math
import sys
import numpy as np

def flux(ketf, ketb):

    kji = ketf 
    kij = ketb 
    kji_r = kji[:]
    kji_r.reverse()
    kij_r = kij[:]
    kij_r.reverse()
    #print '###########################################################'

    print('Rate:')
    for i in range(len(kji)):
        print(" | %4.2e | %4.2e |" %(kji[i], kij[i]))

    sys.stdout.flush()
    sol = hp.solve_flux(kji, kij)
    print("Solutions are: ", sol)
    Jf = sol[-1]
    sys.stdout.flush()

    sol = hp.solve_flux(kij_r, kji_r)
    print("Solutions are:", sol)
    Jb = sol[-1]
    sys.stdout.flush()
    ######################################################################################

    return Jf, Jb
