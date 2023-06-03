#!/bin/bash

## Note that VMD, AmberTools, MCCE and the comand-line pdb_tools 
## need to be installed to use this script
## pbs_tools: https://github.com/haddocking/pdb-tools

echo -n " Enter the name of your prmtop file (omit the .prmtop extention)? "
read name

rm DesolvateForMCCE.in 2> /dev/null
cat >> DesolvateForMCCE.in <<-EOF
parm ${name}.prmtop
trajin min.rst7 
strip :WAT,Na+ 
trajout desolv.pdb
EOF

cpptraj -i DesolvateForMCCE.in > cpptraj.log &
wait

rm AMBER2MCCE.tcl 2> /dev/null
cat >> AMBER2MCCE.tcl <<- EOF 
#input
#------------------------------------------------------------------
 set COORDS [lindex \$argv 0]
 set ResInd "ResIndexing.txt"

 mol new \$COORDS
 set INPUT   [open "\$ResInd" r]

 set CONTENT [read -nonewline \$INPUT]
 close \$INPUT

 set RECORDS [split \$CONTENT "\n"]
 foreach REC \$RECORDS {
   set FIELDS [split \$REC " "]
   lassign \$FIELDS proxCys distCys proxHis distHis OrgHemeID NewHemeID 
   puts "\$NewHemeID"

   #Reunite HEM and PRN ---------------------------------------------
   set HEM [atomselect 0 "resname HEH PRN"]
   \$HEM set resname HEM

   #Relabel Propionates ---------------------------------------------
   set CAA   [atomselect top "resid [expr {\$NewHemeID + 1}] and name CA"]
   \$CAA set name CAA
   set CBA   [atomselect top "resid [expr {\$NewHemeID + 1}] and name CB"]
   \$CBA set name CBA
   set CGA   [atomselect top "resid [expr {\$NewHemeID + 1}] and name CG"]
   \$CGA set name CGA
   set O1A   [atomselect top "resid [expr {\$NewHemeID + 1}] and name O1"]
   \$O1A set name O1A
   set O2A   [atomselect top "resid [expr {\$NewHemeID + 1}] and name O2"]
   \$O2A set name O2A

   set CAD   [atomselect top "resid [expr {\$NewHemeID + 2}] and name CA"]
   \$CAD set name CAD
   set CBD   [atomselect top "resid [expr {\$NewHemeID + 2}] and name CB"]
   \$CBD set name CBD
   set CGD   [atomselect top "resid [expr {\$NewHemeID + 2}] and name CG"]
   \$CGD set name CGD
   set O1D   [atomselect top "resid [expr {\$NewHemeID + 2}] and name O1"]
   \$O1D set name O1D
   set O2D   [atomselect top "resid [expr {\$NewHemeID + 2}] and name O2"]
   \$O2D set name O2D
   #-----------------------------------------------------------------

   #Relabel Proximal His---------------------------------------------
   set CB [atomselect top "resname HEM and resid \$NewHemeID and name CB1"]
   \$CB set resid \$proxHis
   \$CB set name CB
   set CG [atomselect top "resname HEM and resid \$NewHemeID and name CG1"]
   \$CG set resid \$proxHis
   \$CG set name CG 
   set ND1 [atomselect top "resname HEM and resid \$NewHemeID and name ND11"]
   \$ND1 set resid \$proxHis
   \$ND1 set name ND1
   set CD2 [atomselect top "resname HEM and resid \$NewHemeID and name CD21"]
   \$CD2 set resid \$proxHis
   \$CD2 set name CD2 
   set CE1 [atomselect top "resname HEM and resid \$NewHemeID and name CE11"]
   \$CE1 set resid \$proxHis
   \$CE1 set name CE1 
   set NE2 [atomselect top "resname HEM and resid \$NewHemeID and name NE21"]
   \$NE2 set resid \$proxHis
   \$NE2 set name NE2
   #-----------------------------------------------------------------

   #Relabel Distal His ----------------------------------------------
   set CB [atomselect top "resname HEM and resid \$NewHemeID and name CB2"]
   \$CB set resid \$distHis 
   \$CB set name CB
   set CG [atomselect top "resname HEM and resid \$NewHemeID and name CG2"]
   \$CG set resid \$distHis
   \$CG set name CG
   set ND1 [atomselect top "resname HEM and resid \$NewHemeID and name ND12"]
   \$ND1 set resid \$distHis
   \$ND1 set name ND1
   set CD2 [atomselect top "resname HEM and resid \$NewHemeID and name CD22"]
   \$CD2 set resid \$distHis
   \$CD2 set name CD2
   set CE1 [atomselect top "resname HEM and resid \$NewHemeID and name CE12"]
   \$CE1 set resid \$distHis
   \$CE1 set name CE1
   set NE2 [atomselect top "resname HEM and resid \$NewHemeID and name NE22"]
   \$NE2 set resid \$distHis
   \$NE2 set name NE2
   #-----------------------------------------------------------------

   #Relabel Cys attached to B ring ----------------------------------         
   set CB [atomselect top "resname HEM and resid \$NewHemeID and name CBB2"]
   \$CB set resid \$proxCys
   \$CB set name CB
   set SG [atomselect top "resname HEM and resid \$NewHemeID and name SGB2"]
   \$SG set resid \$proxCys
   \$SG set name SG
   #-----------------------------------------------------------------         

   #Relabel Distal Cys ----------------------------------------------         
   set CB [atomselect top "resname HEM and resid \$NewHemeID and name CBC1"]
   \$CB set resid \$distCys
   \$CB set name CB
   set SG [atomselect top "resname HEM and resid \$NewHemeID and name SGC1"]
   \$SG set resid \$distCys
   \$SG set name SG
   #-----------------------------------------------------------------         

   #Reunite HEH and PRN ---------------------------------------------         
   set HEM [atomselect top "resid \$NewHemeID [expr {\$NewHemeID + 1}] [expr {\$NewHemeID + 2}]"]
   \$HEM set resid \$NewHemeID
   #-----------------------------------------------------------------         

   #Re-establish HIS and CYS Residues -------------------------------         
   set HIS [atomselect 0 "resname HIO or resid \$proxHis \$distHis"]
   \$HIS set resname HIL
   set CYS [atomselect 0 "resname CYO or resid \$proxCys \$distCys"]
   \$CYS set resname CYL
   #-----------------------------------------------------------------         

   #Rename solvent residues -----------------------------------------         
   #set SOL [atomselect 0 "not protein and not resname HEM"]
   #\$SOL set resname SOL 
   #------------------------------------------------------------------
}

set sel [atomselect 0 "all and noh"]
\$sel writepdb PrepdforMCCE.pdb

exit
EOF

vmd -dispdev text -e AMBER2MCCE.tcl -args desolv.pdb > vmd.log 2> /dev/null &
wait
pdb_sort -R PrepdforMCCE.pdb > PrepdforMCCEWithSortedRes.pdb

rm SubmitMCCE.sb 2> /dev/null
cat >> SubmitMCCE.sb <<-EOF
#!/bin/bash

step1.py PrepdforMCCEWithSortedRes.pdb
wait

step2.py -l 2 -d 8
wait

step3.py -p 6 -d 8
wait

step4.py -i -300 -t eh -d 60
wait

#Analyze output ---------------------------------------------------------------------------------------------------

i=0
while IFS= read -r F1 <&2;do 
  i=\$((\$i+1))
  A[\$i]=\$(echo \$F1 | awk '{printf $6}')
done 2< ResIndexing.txt

len=\${#A[@]}

echo "=============================================================="
rm RXN-MCCE.txt 2> /dev/null
for (( S=1; S<=\$len; S++ ));do 
  data=\$(grep "HEM+X0\${A[\$S]}_" pK.out | awk '{print \$2}')
  echo "\${A[\$S]} \$data"
  echo "\${A[\$S]} \$data" >> RXN-MCCE.txt
  done
echo "=============================================================="

rm DG.txt 2> /dev/null
for (( S=1; S<\$len; S++ ));do
  S1=\$(grep \${A[\$S]}        RXN-MCCE.txt  | awk '{print \$2}')
  S2=\$(grep \${A[\$((\$S+1))]} RXN-MCCE.txt  | awk '{print \$2}')
  D=\$(echo "(-1*((-1*\$S1) + (\$S2)))"|bc -l)
  printf "%s -> %s: %10s %10s %10s \n" "\${A[\$S]}" "\${A[\$((\$S+1))]}" "\$S1" "\$S2" "\$D"
  printf "%s -> %s: %.4f %.4f %.4f \n" "\${A[\$S]}" "\${A[\$((\$S+1))]}" "\$S1" "\$S2" "\$D" >> DG.txt
done
echo "=============================================================="

#Clean Directory --------------------------------------------------------------------------------------------------
cd energies
rm *.opp
cd ../
rm -r acc.atm
rm -r acc.res
rm -r energies
#rm -r fort.38
rm -r null
rm -r param
rm -r respair.lst
rm -r run.prm
rm -r run.prm.record
rm -r vdw0.lst
cd ../

EOF

