################################################################################################################################################
# Generic Modules
import os
import sys
import shutil
import subprocess
from subprocess import Popen
################################################################################################################################################
# Custom Modules 
import ComputeFlux
import ComputeRedoxCurrent
################################################################################################################################################

################################################################################################################################################
def RedoxCurrentPrediction(LaunchDir, InputDict):

    if (os.path.exists(f"{LaunchDir}/EE") == True):
        EngDir = f"{LaunchDir}/EE"
    else:
        EngDir = f"{LaunchDir}"

    print("""
 This division of the BioDC workflow computes, via the analytical
 Derrida formula [Ref #1], the charge diffuction coefficient 
 based on the non-adiabatic Marcus theory rates estimate in the 
 Energetic Evalulation module. 

 The diffusion coefficient is then related to the electrical 
 resistance and used in Ohm's law to compute the current as a 
 function of applied bias [Ref #2]. Note that this approach is 
 only rigorously correct in the limit of zero bias and single 
 (or non-interacting) mobile charges. The implementaiton of the 
 Derrida formula was kindly provided by Fredrik Jansson 
 [Refs #3 & #4].

 To relax the single-particle condition, the multi-particle 
 steady-state flux is computed according to an analytical 
 expression derived by Blumberger and co-workers (Refs. #5-#7).
 The original code was kindly provided by Jochen Blumberger 
 and Xiuyun Jiang.

 References
   [1] B. Derrida, 
       "Velocity and diffusion constant of a periodic one-dimensional hopping model" 
       J. Stat. Phys. 31, 433 (1983).

   [2] M. J. Guberman-Pfeffer
       "Assessing Thermal Response of Redox Conduction for Anti-Arrhenius Kinetics 
       in a Microbial Cytochrome Nanowire"
       J. Phys. Chem. B 2022, 126, 48, 10083–10097

   [3] Thesis, F. Jansson, Charge transport in disordered materials -
       simulations, theory, and numerical modeling of hopping transport and
       electron-hole recombination. Åbo Akademi University, 2011
       https://urn.fi/URN:NBN:fi-fe201311277464

       Implementation details in section 3.6, especially a method to evaluate
       the expressions in linear time. Application in section 6.2.

   [4] Effect of Electric Field on Diffusion in Disordered Materials I.
       One-dimensional Hopping Transport, A. V. Nenashev, F. Jansson,
       S. D. Baranovskii, R. Österbacka, A. V. Dvurechenskii, F. Gebhard,
       Phys. Rev. B 81, 115203 (2010)
       http://dx.doi.org/10.1103/PhysRevB.81.115203

   [5] M. Breuer, K. M. Rosso, and J. Blumberger, 
       “Electron flow in multi-heme bacterial cytochromes is a balancing act 
       between heme electronic interaction and redox potentials,”
       Proc. Nat. Acad. Sci. USA, vol. 111, p. 611, 2014.

   [6] X. Jiang, Z. Futera, M. E. Ali, F. Gajdos, G. F. von Rudorff, A. Carof, M. Breuer, and J. Blumberger, 
       “Cysteine linkages accelerate electron flow through tetra-heme protein STC,” 
       J. Am. Chem. Soc., vol. 139, p. 17237–17240, 2017.

   [7] X. Jiang, J. H. van Wonderen, J. N. Butt, M. J. Edwards, T. A. Clarke, and J. Blumberger, 
       “Which multi-heme protein complex transfers electrons more efficiently? Comparing MtrCAB from Shewanella 
       with OmcS from Geobacter,” 
       J. Phys. Chem. Lett., vol. 11, pp. 9421-9425, 2020.
    """, end=" ")

    if (os.path.isfile(f"{EngDir}/min.pdb") == True) and (os.path.isfile(f"{EngDir}/rates.txt") == True) and (os.path.isfile(f"{EngDir}/LinearizedHemeSequence.txt") == True):
        print("""
 The following files were found and will be copied to the current working directory:
   min.pdb
   LinearizedHemeSequence.txt will also be copied.
   rates.pdb

 All the information needed for the analysis of 
 redox conductivity is provided by these files,
 so let's proceed!

 We will first compute the multi-particle, stead-state flux.
 """)
        shutil.copy(f"{EngDir}/min.pdb", f"{os.getcwd()}/min.pdb")
        shutil.copy(f"{EngDir}/LinearizedHemeSequence.txt", f"{os.getcwd()}/LinearizedHemeSequence.txt")
        shutil.copy(f"{EngDir}/rates.txt", f"{os.getcwd()}/rates.txt")
        Jf,Jb,Javg = ComputeFlux.ComputeFlux()

        print("""
 We will now (at last) compute the redox current. 
 To do this, some system-specific information is needed.
        """)
        ComputeRedoxCurrent.ComputeRedoxCurrent(Javg, LaunchDir)
        print(" Done!\n")
    else:
        sys.exit(f"""
 min.pdb, rates.txt, and/or LinearizedHemeSequence.txt are missing 
 from {EngDir}. 

 Please run modules #1 and/or #2 of the BioDC program to generate 
 these files. Alternatively, both text files can be created by hand.

 In the case of rates.txt, the syntax on each line should be: 
    [forward rate],[revers rate]
 where each line is for a charge transfer step. The lines
 are assumed to be in the linear sequence of charge 
 transfer steps.

 LinearizedHemeSequence.txt should be a two column file 
 where each line has a zero-based index and the heme 
 residue ID, seaprated by a space. \n""")
################################################################################################################################################
