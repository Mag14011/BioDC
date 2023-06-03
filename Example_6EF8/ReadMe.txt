Prepare6EF8.sh
  Run as: chmod u+x Prepare6EF8; ./Prepare6EF8

  This BASH script uses the Python pdb_tools and VMD
  to fetch and process PDB 6EF8 for use with BioDC

  Output: 6EF8_preped.pdb

ResIndexing.txt
  A critical file that, for this structure, can be 
  generated automatically by specifying a distance
  threshold of 2.5 Å for identifying the His and 
  Cys residues bonded to a given heme group

6EF8.prmtop and min.rst7
  Structures generated from the Structure Relaxation
  and Preparation module of BioDC 

LinearizedHemeSequence.txt
  The linear sequence of hemes specified to study
  redox conduction.

Lambda.txt, DG.txt, Hda.txt, rates.txt, D.txt
  The estimated reorganization energies, 
  reaction free energies, electronic couplings, 
  Marcus-theory rates, and one-dimensional 
  diffusion coefficient computed for charge flow 
  along the heme chain 

SubunitLength.txt
  Measurement of the Fe-toFe distance from the first
  heme in one subunit to the first heme in the 
  subsequent subunit

output.log
  Output from the interactive BioDC session.
