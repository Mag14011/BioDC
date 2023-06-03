#!/bin/bash

 pdb_fetch -biounit 7TFS > 7TFS.pdb
 pdb_splitmodel 7TFS.pdb

 pdb_chain -A 7TFS_1.pdb > 7TFS_A.pdb
 pdb_chain -B 7TFS_2.pdb > 7TFS_B.pdb
 pdb_chain -C 7TFS_3.pdb > 7TFS_C.pdb
 pdb_chain -D 7TFS_4.pdb > 7TFS_D.pdb

 pdb_reres 7TFS_A.pdb > 7TFS_Areres.pdb
 pdb_reres 7TFS_B.pdb > 7TFS_Breres.pdb
 pdb_reres 7TFS_C.pdb > 7TFS_Creres.pdb
 pdb_reres 7TFS_D.pdb > 7TFS_Dreres.pdb

 cat 7TFS_Areres.pdb 7TFS_Breres.pdb 7TFS_Creres.pdb 7TFS_Dreres.pdb | grep -v "MG    MG"  > 7TFS_trimer.pdb

 rm EditPDB.tcl 2> /dev/null

cat >> EditPDB.tcl <<- EOF
mol new 7TFS_trimer.pdb
set prot [atomselect top "(chain B C D and not resname HEC) or (chain A and resid 1 to 32)"]
set HEC  [atomselect top "chain B C D and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

 vmd -e EditPDB.tcl > EditPDB.log

 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 7TFS_trimer_reord.pdb
 pdb_reres 7TFS_trimer_reord.pdb > 7TFS_trimer_reord_reres.pdb
 mv 7TFS_trimer_reord_reres.pdb 7TFS_preped.pdb


