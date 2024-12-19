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
 pdb_fetch -biounit 8E5G > 8E5G.pdb

# 2) Split the PDB by chain 
 pdb_splitmodel 8E5G.pdb
 pdb_chain -A 8E5G_1.pdb > 8E5G_A.pdb
 pdb_chain -B 8E5G_2.pdb > 8E5G_B.pdb
 pdb_chain -C 8E5G_3.pdb > 8E5G_C.pdb
 pdb_chain -D 8E5G_4.pdb > 8E5G_D.pdb
 pdb_chain -E 8E5G_5.pdb > 8E5G_E.pdb

# 3) Renumber the residues sequentially
 pdb_reres 8E5G_A.pdb > 8E5G_Areres.pdb
 pdb_reres 8E5G_B.pdb > 8E5G_Breres.pdb
 pdb_reres 8E5G_C.pdb > 8E5G_Creres.pdb
 pdb_reres 8E5G_D.pdb > 8E5G_Dreres.pdb
 pdb_reres 8E5G_E.pdb > 8E5G_Ereres.pdb

# 4) Concatonate the chains in the order of their geometrical linear sequence
 cat 8E5G_Ereres.pdb 8E5G_Dreres.pdb 8E5G_Areres.pdb 8E5G_Creres.pdb 8E5G_Breres.pdb > 8E5G_trimer.pdb

# 5) Remove, if presnet, a tcl script for editing the concatonated PDB
 rm EditPDB.tcl 2> /dev/null

# 6) Write a TCL script to separate protein and heme portions of the 
# concatonated PDB into separate PDBs. 
cat >> EditPDB.tcl <<- EOF
mol new 8E5G_trimer.pdb
set prot [atomselect top "(chain A B C D E and not resname HEC)"]
set HEC  [atomselect top "(chain E and resname HEC and not resid 221) or (chain D and resname HEC) or (chain A and resname HEC) or (chain C and resname HEC) or (chain B and resname HEC and not resid 222 224)"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 7) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 8) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 8E5G_trimer_reord.pdb

# 9) Renumber PDB with sequential residue numbering
 pdb_reres 8E5G_trimer_reord.pdb > 8E5G_trimer_reord_reres.pdb

#10) Add missing OXT atoms and rename file
 pdbfixer 8E5G_trimer_reord_reres.pdb --add-atoms=heavy
 mv output.pdb 8E5G_preped.pdb

#11) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(8E5G_preped.pdb|Prepare8E5G.sh|prep) prep/

#12) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"

