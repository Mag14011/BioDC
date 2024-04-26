################################################################################################################################################
# Generic Modules
import os
import sys
import subprocess
from subprocess import Popen
from pandas import read_csv
################################################################################################################################################
# Custom Modules 
import derrida
################################################################################################################################################

def ComputeDiffusionCoefficient(AvgHemeSpacing):

    data = read_csv("rates.txt")
    ketf = data['ketf'].tolist()
    ketb = data['ketb'].tolist()
    V,D = derrida.VD(ketf, ketb)
    print("  Diffusion constant = %E cm^2/S"  % (D * ((AvgHemeSpacing)**2))) 
    print("Diffusion constant = %E (cm^2/S)" % (D * ((AvgHemeSpacing)**2)), file=open('D.txt', 'w'))
################################################################################################################################################
