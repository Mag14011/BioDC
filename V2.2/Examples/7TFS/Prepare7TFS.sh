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
 pdb_fetch -biounit 7TFS > 7TFS.pdb

# 2) Split the PDB by chain 
 pdb_splitmodel 7TFS.pdb
 pdb_chain -A 7TFS_1.pdb > 7TFS_A.pdb
 pdb_chain -B 7TFS_2.pdb > 7TFS_B.pdb
 pdb_chain -C 7TFS_3.pdb > 7TFS_C.pdb
 pdb_chain -D 7TFS_4.pdb > 7TFS_D.pdb

# 3) Renumber the chains
 pdb_reres 7TFS_A.pdb > 7TFS_Areres.pdb
 pdb_reres 7TFS_B.pdb > 7TFS_Breres.pdb
 pdb_reres 7TFS_C.pdb > 7TFS_Creres.pdb
 pdb_reres 7TFS_D.pdb > 7TFS_Dreres.pdb

# 4) Concatonate the chain in the linear geometric sequence
 cat 7TFS_Areres.pdb 7TFS_Breres.pdb 7TFS_Creres.pdb 7TFS_Dreres.pdb | grep -v "MG    MG"  > 7TFS_trimer.pdb

# 5) Remove, if presnet, a tcl script for editing the concatonated PDB
 rm EditPDB.tcl 2> /dev/null

# 6) Write a TCL script to separate protein and heme portions of the 
# concatonated PDB into separate PDBs. This script also only takes 
# enough of chain D to provide the inter-chain ligand to a heme 
# in chain A.
cat >> EditPDB.tcl <<- EOF
mol new 7TFS_trimer.pdb
set prot [atomselect top "(chain B C D and not resname HEC) or (chain A and resid 1 to 32)"]
set HEC  [atomselect top "chain B C D and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 7) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 8) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 7TFS_trimer_reord.pdb

# 9) Renumber PDB with sequential residue numbering
 pdb_reres 7TFS_trimer_reord.pdb > 7TFS_trimer_reord_reres.pdb

#10) Change file name
 mv 7TFS_trimer_reord_reres.pdb 7TFS_preped.pdb

#11) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(7TFS_preped.pdb|Prepare7TFS.sh|prep) prep/

#12) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"



