#!/bin/bash

 pdb_fetch -biounit 8D9M > 8D9M.pdb
 pdb_splitmodel 8D9M.pdb

 pdb_chain -A 8D9M_1.pdb > 8D9M_A.pdb
 pdb_chain -B 8D9M_2.pdb > 8D9M_B.pdb
 pdb_chain -C 8D9M_3.pdb > 8D9M_C.pdb

 pdb_reres 8D9M_A.pdb > 8D9M_Areres.pdb
 pdb_reres 8D9M_B.pdb > 8D9M_Breres.pdb
 pdb_reres 8D9M_C.pdb > 8D9M_Creres.pdb

 cat 8D9M_Areres.pdb 8D9M_Breres.pdb 8D9M_Creres.pdb > 8D9M_trimer.pdb

 rm EditPDB.tcl 2> /dev/null

cat >> EditPDB.tcl <<- EOF
mol new 8D9M_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

 vmd -e EditPDB.tcl > EditPDB.log

 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 8D9M_trimer_reord.pdb
 pdb_reres 8D9M_trimer_reord.pdb > 8D9M_trimer_reord_reres.pdb
 mv 8D9M_trimer_reord_reres.pdb 8D9M_preped.pdb


