#!/bin/bash

##################################################################################
# This script obtains the relevant PDB, selects what portion of the biological   #
# assembly is desired for modeling, reorders and renumbers the residues so all   #
# the protein chains come before any cofactors with sequential residue numbering #
# and renames for final PDB.                                                     #
#                                                                                #
# The reordering and renumbering of the residues are essential for BioDC to      #
# handle the PDB.                                                                #
#                                                                                #
# The disulfide-linked Cys residues do not have to be re-labeled at this stage.  #  
# BioDC will ask if there are disulfide-linked Cys residues and if so, what are  #
# the pairs of residue IDs. However, it is easier to re-label the Cys residues   #
# beforehand, and to note the residue IDs in the prepared PDB submitted to BioDC #
# than to figure out after all the residue re-numbering, what Cys residues the   #
# published IDs now correspond to in the prepared structure for BioDC. The       #
# residue IDs of these Cys residues still need to be given to BioDC so the       # 
# appropriate bond definitions can be specified to TLEaP                         #.
##################################################################################

# 1) Fetch biological assembly from PDB
 pdb_fetch -biounit 8E5F > 8E5F.pdb

# 2) Split the PDB by chain 
 pdb_splitmodel 8E5F.pdb
 pdb_chain -A 8E5F_1.pdb > 8E5F_A.pdb
 pdb_chain -B 8E5F_2.pdb > 8E5F_B.pdb
 pdb_chain -C 8E5F_3.pdb > 8E5F_C.pdb
 pdb_chain -D 8E5F_4.pdb > 8E5F_D.pdb

# 3) Re-name Cys residues involved in 
#    disulfide bonds according to AMBER conventions
rm LabelCYS.tcl 2> /dev/null
cat >> LabelCYS.tcl <<- EOF
mol new 8E5F_A.pdb
set sel [atomselect 0 "resname CYS and resid 62 69 133 163"]
\$sel set resname CYX
set all [atomselect 0 "all"]
\$all writepdb 8E5F_Acyx.pdb

mol new 8E5F_B.pdb
set sel [atomselect 1 "resname CYS and resid 62 69 133 163"]
\$sel set resname CYX
set all [atomselect 1 "all"]
\$all writepdb 8E5F_Bcyx.pdb

mol new 8E5F_C.pdb
set sel [atomselect 2 "resname CYS and resid 62 69 133 163"]
\$sel set resname CYX
set all [atomselect 2 "all"]
\$all writepdb 8E5F_Ccyx.pdb

mol new 8E5F_D.pdb
set sel [atomselect 3 "resname CYS and resid 62 69 133 163"]
\$sel set resname CYX
set all [atomselect 3 "all"]
\$all writepdb 8E5F_Dcyx.pdb

exit
EOF

# 4) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e LabelCYS.tcl > LabelCYS.log

# 5) Renumber residues sequentially
 pdb_reres 8E5F_Acyx.pdb > 8E5F_Areres.pdb
 pdb_reres 8E5F_Bcyx.pdb > 8E5F_Breres.pdb
 pdb_reres 8E5F_Ccyx.pdb > 8E5F_Creres.pdb
 pdb_reres 8E5F_Dcyx.pdb > 8E5F_Dreres.pdb

# 6) Concatonate chains in the linear geometric sequence
 cat 8E5F_Dreres.pdb 8E5F_Areres.pdb 8E5F_Creres.pdb 8E5F_Breres.pdb | egrep -v "CRYST1|END" > 8E5F_trimer.pdb

# 7) Write a TCL script to separate protein and heme portions of the 
#    concatonated PDB into separate PDBs. 
rm EditPDB.tcl 2> /dev/null
cat >> EditPDB.tcl <<- EOF
mol new 8E5F_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC) or (chain D and resid 20 to 60)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

# 8) Run the TCL script with Visual Molecular Dynamics (VMD)
 vmd -e EditPDB.tcl > EditPDB.log

# 9) Concatonate the protein and heme PDBs and remove CRYST1 and END lines
 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 8E5F_trimer_reord.pdb

#10) Renumber PDB with sequential residue numbering
 pdb_reres 8E5F_trimer_reord.pdb > 8E5F_trimer_reord_reres.pdb

#11) Change file name
 mv 8E5F_trimer_reord_reres.pdb 8E5F_preped.pdb

#12) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(8E5F_preped.pdb|Prepare8E5F.sh|prep) prep/

#12) Get residue IDs of the disulfide-linked Cys residues
if [ -f 8E5F_preped.pdb ];then
  count=$(grep "CA  CYX" 8E5F_preped.pdb | wc -l)
  PairNum=$(awk "BEGIN {print $count/2}")

  for i in $(seq 1 1 $count);do
    R[i]=$(grep -m${i} "CA  CYX" 8E5F_preped.pdb | tail -n1 | cut -b 23-26)
  done

  echo " 
  When running the first module of BioDC, please indicate that these are 
  $PairNum disulfide bonds. The residue IDs of the particiapting Cys 
  residues are: "

  echo "  ${R[@]}"

  echo "
  In this particular case, every two residue IDs form a disulfide-linked 
  pair, but this may not be the case for other structures."

fi

#13) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"
