#!/bin/bash

##################################################################################
# This script obtains the relevant PDB, selects what portion of the biological   #
# assembly is desired for modeling, reorders and renumbers the residues so all   #
# the protein chains come before any cofactors with sequential residue numbering #
# and renames the final PDB.                                                     #
#                                                                                #
# The reordering and renumbering of the residues are essential for BioDC to      #
# process the PDB.                                                                #
##################################################################################

# 1) Fetch biological assembly from PDB
 pdb_fetch -biounit 6EF8 > 6EF8.pdb

# 2) Split the PDB by chain 
 pdb_splitchain 6EF8.pdb

# 3) Concatonate the chains in the order of their geometrical linear sequence
 cat 6EF8_D.pdb 6EF8_B.pdb 6EF8_A.pdb 6EF8_C.pdb > 6EF8_trimer.pdb

# 4) Remove, if presnet, a tcl script for editing the concatonated PDB
 rm EditPDB.tcl 2> /dev/null

# 5) Write a TCL script to separate protein and heme portions of the 
#    concatonated PDB into separate PDBs. This script also only takes 
#    enough of chain D to provide the inter-chain ligand to a heme 
#    in chain A.

rm EditPDB.tcl 2> /dev/null
cat >> EditPDB.tcl <<- EOF
mol new 6EF8_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC) or (chain D and resid 1 to 20)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 6) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 7) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 6EF8_trimer_reord.pdb

# 8) Renumber PDB with sequential residue numbering
 pdb_reres 6EF8_trimer_reord.pdb > 6EF8_trimer_reord_reres.pdb

# 9) Change file name
 mv 6EF8_trimer_reord_reres.pdb 6EF8_preped.pdb

#10) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(6EF8_preped.pdb) prep/

#11) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"





