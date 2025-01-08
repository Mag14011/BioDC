# biodc/core/preparation/select_mutate.py

"""
Module for selecting and performing mutations on protein structures using BioPython.
"""

from pathlib import Path
from typing import Dict, List, Tuple, Union
import logging
from Bio import PDB
from biodc.utils.interaction import InteractionManager

logger = logging.getLogger(__name__)

class StructureMutator:
    def __init__(self, pdb_path: str):
        self.pdb_path = pdb_path
        self.parser = PDB.PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', f"{pdb_path}.pdb")

    def validate_mutation(self, orig_res: str, res_id: str) -> bool:
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if (str(residue.get_id()[1]) == res_id and
                        residue.get_resname() == orig_res):
                        return True
        return False

    def apply_mutations(self, mutations: List[Tuple[str, str, str]]):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.get_id()[1]
                    for orig_res, mut_id, new_res in mutations:
                        if str(res_id) == mut_id and residue.get_resname() == orig_res:
                            backbone_atoms = ['N', 'CA', 'C', 'O']
                            atoms_to_remove = [atom for atom in residue
                                             if atom.get_id() not in backbone_atoms]
                            for atom in atoms_to_remove:
                                residue.detach_child(atom.get_id())
                            residue.resname = new_res

    def save_structure(self, output_path: Union[str, Path]):
        io = PDB.PDBIO()
        io.set_structure(self.structure)
        output_path = str(output_path)
        io.save(output_path)

def select_mutate(pdb: str, launch_dir: Path, input_dict: Dict) -> str:
    interaction_manager = InteractionManager(
        launch_dir=launch_dir,
        input_dict=input_dict
    )

    num_mut = interaction_manager.prompt(
        "NumMut",
        "\n How many residues do you want to mutate?",
        input_type=int
    )

    if num_mut != 0:
        mutator = StructureMutator(pdb)
        mutations = []

        print("""
 For each mutation, please enter the residue three-letter code
 before the mutation, the residue ID, and the three-letter code
 after the mutation, each separated by a space.

 The program will keep the backbone atoms (N, CA, C, O) of the
 original residue, change its name, and later TLEaP will build
 the new sidechain based on the appropriate template in the
 selected force field library.""")

        for idx in range(num_mut):
            while True:
                sel_res = interaction_manager.prompt(
                    f"Mutation_{idx+1}",
                    f"  Mutation {idx+1}"
                )

                mut_res_array = sel_res.split()
                if len(mut_res_array) != 3:
                    print(" Invalid input format. Please provide: ORIGINAL_RES ID NEW_RES")
                    continue

                orig_res, res_id, new_res = mut_res_array
                if not mutator.validate_mutation(orig_res, res_id):
                    print(f" Error: No residue found with ID {res_id} and name {orig_res}")
                    continue

                break

            mutations.append(tuple(mut_res_array))

        print("\n Generating PDB for the mutated structure...")
        try:
            mutator.apply_mutations(mutations)
            mutator.save_structure(launch_dir / "SPR" / "mutated.pdb")
            print(" Structure mutation completed successfully.")
            return "mutated"
        except Exception as e:
            logger.error(f" Error during mutation: {str(e)}")
            raise

    return pdb
