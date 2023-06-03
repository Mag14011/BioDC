#!/bin/bash

 pdb_fetch -biounit 6NEF > 6NEF.pdb
 pdb_splitmodel 6NEF.pdb

 pdb_chain -A 6NEF_1.pdb > 6NEF_A.pdb
 pdb_chain -B 6NEF_2.pdb > 6NEF_B.pdb
 pdb_chain -C 6NEF_3.pdb > 6NEF_C.pdb
 pdb_chain -D 6NEF_4.pdb > 6NEF_D.pdb

 pdb_reres 6NEF_A.pdb > 6NEF_Areres.pdb
 pdb_reres 6NEF_B.pdb > 6NEF_Breres.pdb
 pdb_reres 6NEF_C.pdb > 6NEF_Creres.pdb
 pdb_reres 6NEF_D.pdb > 6NEF_Dreres.pdb

 cat 6NEF_Dreres.pdb 6NEF_Breres.pdb 6NEF_Areres.pdb 6NEF_Creres.pdb | grep -v "MG    MG"  > 6NEF_trimer.pdb

 rm EditPDB.tcl 2> /dev/null

cat >> EditPDB.tcl <<- EOF
mol new 6NEF_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC) or (chain D and resid 1 to 20)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

 vmd -e EditPDB.tcl > EditPDB.log

 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 6NEF_trimer_reord.pdb
 pdb_reres 6NEF_trimer_reord.pdb > 6NEF_trimer_reord_reres.pdb
 mv 6NEF_trimer_reord_reres.pdb 6NEF_preped.pdb

