################################################################################################################################################
# Generic Modules
from pandas import read_csv
################################################################################################################################################
# Custom Modules 
import blumberger
################################################################################################################################################

################################################################################################################################################

def ComputeFlux():

    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()

    Jf,Jb = blumberger.flux(ketf, ketb)
    Javg = (Jf+Jb)/2

    print("Forward Flux: %.2E" %(Jf)) 
    print("Reverse Flux: %.2E" %(Jb))
    print("Average forward/backward Flux: %.2E" %(Javg))

    return Jf, Jb, Javg
################################################################################################################################################
