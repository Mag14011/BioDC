#!/bin/bash

 pdb_fetch -biounit 8E5G > 8E5G.pdb
 pdb_splitmodel 8E5G.pdb

 pdb_chain -A 8E5G_1.pdb > 8E5G_A.pdb
 pdb_chain -B 8E5G_2.pdb > 8E5G_B.pdb
 pdb_chain -C 8E5G_3.pdb > 8E5G_C.pdb
 pdb_chain -D 8E5G_4.pdb > 8E5G_D.pdb
 pdb_chain -E 8E5G_5.pdb > 8E5G_E.pdb

 pdb_reres 8E5G_A.pdb > 8E5G_Areres.pdb
 pdb_reres 8E5G_B.pdb > 8E5G_Breres.pdb
 pdb_reres 8E5G_C.pdb > 8E5G_Creres.pdb
 pdb_reres 8E5G_D.pdb > 8E5G_Dreres.pdb
 pdb_reres 8E5G_E.pdb > 8E5G_Ereres.pdb

 cat 8E5G_Ereres.pdb 8E5G_Dreres.pdb 8E5G_Areres.pdb 8E5G_Creres.pdb 8E5G_Breres.pdb > 8E5G_trimer.pdb

 rm EditPDB.tcl 2> /dev/null

cat >> EditPDB.tcl <<- EOF
mol new 8E5G_trimer.pdb
set prot [atomselect top "(chain A B C D E and not resname HEC)"]
set HEC  [atomselect top "(chain E and resname HEC and not resid 221) or (chain D and resname HEC) or (chain A and resname HEC) or (chain C and resname HEC) or (chain B and resname HEC and not resid 222 224)"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

 vmd -e EditPDB.tcl > EditPDB.log

 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 8E5G_trimer_reord.pdb
 pdb_reres 8E5G_trimer_reord.pdb > 8E5G_trimer_reord_reres.pdb
 mv 8E5G_trimer_reord_reres.pdb 8E5G_preped.pdb

