#!/bin/bash

##################################################################################
# This script obtains the relevant PDB, selects what portion of the biological   #
# assembly is desired for modeling, reorders and renumbers the residues so all   #
# the protein chains come before any cofactors with sequential residue numbering #
# and renames for final PDB.                                                     #
#                                                                                #
# The reordering and renumbering of the residues are essential for BioDC to      #
# handle the PDB.                                                                #
##################################################################################

# 1) Fetch biological assembly from PDB
 pdb_fetch -biounit 8D9M > 8D9M.pdb

# 2) Split the PDB by chain
 pdb_splitmodel 8D9M.pdb
 pdb_chain -A 8D9M_1.pdb > 8D9M_A.pdb
 pdb_chain -B 8D9M_2.pdb > 8D9M_B.pdb
 pdb_chain -C 8D9M_3.pdb > 8D9M_C.pdb

# 3) Renumber the residues sequentially
 pdb_reres 8D9M_A.pdb > 8D9M_Areres.pdb
 pdb_reres 8D9M_B.pdb > 8D9M_Breres.pdb
 pdb_reres 8D9M_C.pdb > 8D9M_Creres.pdb

# 4) Concatonate the chains in the linear geometric sequence
 cat 8D9M_Areres.pdb 8D9M_Breres.pdb 8D9M_Creres.pdb > 8D9M_trimer.pdb

# 5) Remove, if presnet, a tcl script for editing the concatonated PDB
 rm EditPDB.tcl 2> /dev/null

# 6) Write a TCL script to separate protein and heme portions of the 
#    concatonated PDB into separate PDBs. This script also only takes 
#    enough of chain D to provide the inter-chain ligand to a heme 
#    in chain A.
cat >> EditPDB.tcl <<- EOF
mol new 8D9M_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 7) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 8) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 8D9M_trimer_reord.pdb

# 9) Renumber PDB with sequential residue numbering
 pdb_reres 8D9M_trimer_reord.pdb > 8D9M_trimer_reord_reres.pdb

#10) Change file name
 mv 8D9M_trimer_reord_reres.pdb 8D9M_preped.pdb

#11) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(8D9M_preped.pdb|Prepare8D9M.sh|prep) prep/

# 5) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"


