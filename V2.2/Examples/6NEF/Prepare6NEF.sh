#!/bin/bash

##################################################################################
# This script obtains the relevant PDB, selects what portion of the biological   #
# assembly is desired for modeling, reorders and renumbers the residues so all   #
# the protein chains come before any cofactors with sequential residue numbering #
# and renames for final PDB.                                                     #
#                                                                                #
# The reordering and renumbering of the residues are essential for BioDC to      #
# process the PDB.                                                                #
##################################################################################

# 1) Fetch biological assembly from PDB
 pdb_fetch -biounit 6NEF > 6NEF.pdb

# 2) Split the PDB by chain 
 pdb_splitmodel 6NEF.pdb
 pdb_chain -A 6NEF_1.pdb > 6NEF_A.pdb
 pdb_chain -B 6NEF_2.pdb > 6NEF_B.pdb
 pdb_chain -C 6NEF_3.pdb > 6NEF_C.pdb
 pdb_chain -D 6NEF_4.pdb > 6NEF_D.pdb

# 3) Renumber the chains
 pdb_reres 6NEF_A.pdb > 6NEF_Areres.pdb
 pdb_reres 6NEF_B.pdb > 6NEF_Breres.pdb
 pdb_reres 6NEF_C.pdb > 6NEF_Creres.pdb
 pdb_reres 6NEF_D.pdb > 6NEF_Dreres.pdb

# 4) Concatonate the chain in the linear geometric sequence
 cat 6NEF_Dreres.pdb 6NEF_Breres.pdb 6NEF_Areres.pdb 6NEF_Creres.pdb | grep -v "MG    MG"  > 6NEF_trimer.pdb

# 5) Remove, if present, a tcl script for editing the concatonated PDB
 rm EditPDB.tcl 2> /dev/null

# 6) Write a TCL script to separate protein and heme portions of the 
# concatonated PDB into separate PDBs. This script also only takes 
# enough of chain D to provide the inter-chain ligand to a heme 
# in chain A.
cat >> EditPDB.tcl <<- EOF
mol new 6NEF_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC) or (chain D and resid 1 to 20)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 7) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 8) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 6NEF_trimer_reord.pdb

# 9) Renumber PDB with sequential residue numbering
 pdb_reres 6NEF_trimer_reord.pdb > 6NEF_trimer_reord_reres.pdb

#10) Add missing OXT atoms and rename file
 pdbfixer 6NEF_trimer_reord_reres.pdb --add-atoms=heavy
 mv output.pdb 6NEF_preped.pdb

#11) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(6NEF_preped.pdb|Prepare6NEF.sh|prep) prep/

 #12) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"

