#!/bin/bash

 pdb_fetch -biounit 7LQ5 > 7LQ5.pdb
 pdb_reres 7LQ5.pdb > 7LQ5_reres.pdb
 mv 7LQ5_reres.pdb 7LQ5_preped.pdb

