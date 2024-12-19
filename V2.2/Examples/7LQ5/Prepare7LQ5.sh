#!/bin/bash

##################################################################################
# This script obtains the relevant PDB, renumbers the residues sequentially and  #
# renames the PDB.                                                               #
#                                                                                #
# Note that the deposited PDB already lists all the hemes after all the protein  #
# residues, so no re-ordering is needed before using BioDC.                      #
##################################################################################

# 1) Fetch biological assembly from PDB
 pdb_fetch -biounit 7LQ5 > 7LQ5.pdb
 
# 2) Renumber PDB with sequential residue numbering
 pdb_reres 7LQ5.pdb > 7LQ5_reres.pdb

# 3) Change file name
 mv 7LQ5_reres.pdb 7LQ5_preped.pdb

# 4) Clean directory
 mkdir -p prep 2> /dev/null
 shopt -s extglob
 mv !(7LQ5_preped.pdb|Prepare7LQ5.sh|prep) prep/

# 5) Launch BioDC
python "/workhorse/BioDCProblem/BioDC/V2.1/BioDCv2.py"

