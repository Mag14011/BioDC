#!/bin/bash

 pdb_fetch -biounit 6EF8 > 6EF8.pdb

 pdb_splitchain 6EF8.pdb
 cat 6EF8_D.pdb 6EF8_B.pdb 6EF8_A.pdb 6EF8_C.pdb > 6EF8_trimer.pdb

 rm EditPDB.tcl 2> /dev/null

cat >> EditPDB.tcl <<- EOF
mol new 6EF8_trimer.pdb
set prot [atomselect top "(chain A B C and not resname HEC) or (chain D and resid 1 to 20)"]
set HEC  [atomselect top "chain A B C and resname HEC"]
\$prot writepdb prot.pdb
\$HEC writepdb HEC.pdb
exit
EOF

 vmd -e EditPDB.tcl > EditPDB.log

 cat prot.pdb HEC.pdb | egrep -v "CRYST1|END" > 6EF8_trimer_reord.pdb
 pdb_reres 6EF8_trimer_reord.pdb > 6EF8_trimer_reord_reres.pdb
 mv 6EF8_trimer_reord_reres.pdb 6EF8_preped.pdb

