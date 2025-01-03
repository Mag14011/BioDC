import numpy as np
from typing import List, Tuple, Dict, Optional, Set
from sklearn.decomposition import PCA
import argparse
import sys
import json
from pathlib import Path
import networkx as nx
from collections import defaultdict

class PDBProcessor:
    def __init__(self):
        self.heme_atoms = [
            'FE', 'NA', 'C1A', 'C2A', 'C3A', 'C4A', 'CHB', 'C1B', 'NB',
            'C2B', 'C3B', 'C4B', 'CHC', 'C1C', 'NC', 'C2C', 'C3C', 'C4C',
            'CHD', 'C1D', 'ND', 'C2D', 'C3D', 'C4D', 'CHA'
        ]
        self.core_atoms = [
            'NA', 'NB', 'NC', 'ND',
            'C1A', 'C4A', 'C1B', 'C4B', 'C1C', 'C4C', 'C1D', 'C4D'
        ]
        self.distance_cutoff = 20.0

    def read_linearized_sequence(self, filename: str) -> List[Tuple[int, int]]:
        """Read the linearized heme sequence file"""
        pairs = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    _, resid = map(int, line.strip().split())
                    pairs.append(resid)
                except ValueError:
                    continue
        return list(zip(pairs[:-1], pairs[1:]))

    def read_pdb_atoms(self, pdb_file: str) -> Dict[int, Dict[str, np.ndarray]]:
        """Read atom coordinates from PDB file"""
        atoms_dict = {}
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM  ') or line.startswith('HETATM'):
                    try:
                        atom_name = line[12:16].strip()
                        resid = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        if atom_name in self.heme_atoms:
                            if resid not in atoms_dict:
                                atoms_dict[resid] = {}
                            atoms_dict[resid][atom_name] = np.array([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        return atoms_dict

    def calculate_fe_distances(self, atoms_dict: Dict[int, Dict[str, np.ndarray]]) -> Dict[Tuple[int, int], float]:
        """Calculate distances between all pairs of heme FE atoms."""
        distances = {}
        residues = list(atoms_dict.keys())
        
        for i, res1 in enumerate(residues):
            for res2 in residues[i+1:]:
                if 'FE' in atoms_dict[res1] and 'FE' in atoms_dict[res2]:
                    dist = np.linalg.norm(atoms_dict[res1]['FE'] - atoms_dict[res2]['FE'])
                    distances[(res1, res2)] = dist
                    distances[(res2, res1)] = dist
                    
        return distances

    def find_neighbors(self, distances: Dict[Tuple[int, int], float], 
                      cutoff: float) -> Dict[int, Set[int]]:
        """Find neighboring hemes within cutoff distance."""
        neighbors = defaultdict(set)
        for (res1, res2), dist in distances.items():
            if dist <= cutoff:
                neighbors[res1].add(res2)
                neighbors[res2].add(res1)
        return neighbors

    def detect_linear_sequence(self, atoms_dict: Dict[int, Dict[str, np.ndarray]]) -> List[int]:
        """
        Detect linear sequence of hemes by finding termini and building path.
        Assumes linear topology with two termini.
        """
        distances = self.calculate_fe_distances(atoms_dict)
        neighbors = self.find_neighbors(distances, self.distance_cutoff)
        
        termini = [res for res, neighs in neighbors.items() if len(neighs) == 1]
        
        if len(termini) != 2:
            raise ValueError(f"Expected 2 termini in linear sequence, found {len(termini)}. "
                           "Structure may not be linear or distance cutoff may need adjustment.")
        
        sequence = [termini[0]]
        current = termini[0]
        visited = {termini[0]}
        
        while current != termini[1]:
            next_options = neighbors[current] - visited
            if not next_options:
                raise ValueError("Path broken - cannot reach second terminus")
                
            next_heme = min(next_options, 
                          key=lambda x: distances[(current, x)])
            
            sequence.append(next_heme)
            visited.add(next_heme)
            current = next_heme
            
        return sequence

    def detect_branched_sequence(self, atoms_dict: Dict[int, Dict[str, np.ndarray]]) -> List[List[int]]:
        """
        Detect branched sequence of hemes using graph analysis.
        Returns list of branch sequences.
        """
        distances = self.calculate_fe_distances(atoms_dict)
        
        G = nx.Graph()
        for (res1, res2), dist in distances.items():
            if dist <= self.distance_cutoff:
                G.add_edge(res1, res2, weight=dist)
                
        branch_points = [node for node, degree in G.degree() if degree > 2]
        
        if not branch_points:
            raise ValueError("No branch points found. Structure may be linear.")
            
        terminal_nodes = [node for node, degree in G.degree() if degree == 1]
        paths = []
        
        for start in terminal_nodes:
            for end in terminal_nodes:
                if start < end:
                    path = nx.shortest_path(G, start, end, weight='weight')
                    paths.append(path)
                    
        paths.sort(key=len, reverse=True)
        final_paths = []
        covered_edges = set()
        
        for path in paths:
            path_edges = set(zip(path[:-1], path[1:]))
            if not path_edges.issubset(covered_edges):
                final_paths.append(path)
                covered_edges.update(path_edges)
                
        return final_paths

    def detect_sequence(self, atoms_dict: Dict[int, Dict[str, np.ndarray]], 
                       topology: str = 'linear') -> List[int]:
        """
        Detect sequence of hemes using specified topology.
        """
        if topology == 'linear':
            return self.detect_linear_sequence(atoms_dict)
        elif topology == 'branched':
            return self.detect_branched_sequence(atoms_dict)
        else:
            raise ValueError(f"Unknown topology: {topology}")

    def fit_plane_to_heme(self, heme_coords: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit a plane to heme group atoms using PCA.
        """
        coords_list = []
        for atom in self.core_atoms:
            if atom in heme_coords:
                coords_list.append(heme_coords[atom])
        
        if len(coords_list) < 3:
            raise ValueError("Not enough atoms to fit plane")
            
        coords_array = np.array(coords_list)
        center = np.mean(coords_array, axis=0)
        centered_coords = coords_array - center
        
        pca = PCA(n_components=3)
        pca.fit(centered_coords)
        normal_vector = pca.components_[2]
        
        if 'FE' in heme_coords:
            fe_vector = heme_coords['FE'] - center
            if np.dot(normal_vector, fe_vector) < 0:
                normal_vector = -normal_vector
                
        return normal_vector, center

    def calculate_plane_angle(self, coords1: Dict[str, np.ndarray],
                            coords2: Dict[str, np.ndarray]) -> float:
        """Calculate the angle between two heme group planes."""
        normal1, _ = self.fit_plane_to_heme(coords1)
        normal2, _ = self.fit_plane_to_heme(coords2)
        
        cos_angle = np.dot(normal1, normal2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        
        angle_deg = np.degrees(angle_rad)
        if angle_deg > 90:
            angle_deg = 180 - angle_deg
            
        return angle_deg

    def calculate_min_distance(self, coords1: Dict[str, np.ndarray],
                             coords2: Dict[str, np.ndarray]) -> float:
        """Calculate minimum distance between two sets of coordinates"""
        min_dist = float('inf')
        
        for atom1 in self.heme_atoms:
            if atom1 not in coords1:
                continue
            for atom2 in self.heme_atoms:
                if atom2 not in coords2:
                    continue
                    
                dist = np.linalg.norm(coords1[atom1] - coords2[atom2])
                min_dist = min(min_dist, dist)
        
        return min_dist

    def measure_heme_angles(self, pdb_file: str, 
                          sequence_file: Optional[str] = None,
                          topology: Optional[str] = None) -> List[Dict]:
        """Calculate angles between consecutive heme pairs."""
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        if sequence_file:
            residue_pairs = self.read_linearized_sequence(sequence_file)
        else:
            if not topology:
                raise ValueError("Either sequence_file or topology must be specified")
            sequence = self.detect_sequence(atoms_dict, topology)
            if topology == 'linear':
                residue_pairs = list(zip(sequence[:-1], sequence[1:]))
            else:
                residue_pairs = []
                for branch in sequence:
                    residue_pairs.extend(list(zip(branch[:-1], branch[1:])))
        
        angles = []
        for res1, res2 in residue_pairs:
            if res1 in atoms_dict and res2 in atoms_dict:
                try:
                    angle = self.calculate_plane_angle(atoms_dict[res1], atoms_dict[res2])
                    angles.append({'pair': [res1, res2], 'angle': angle})
                except ValueError as e:
                    print(f"Warning: Could not calculate angle for residues {res1}-{res2}: {str(e)}",
                          file=sys.stderr)
        
        return angles

    def measure_avg_dist(self, pdb_file: str,
                        sequence_file: Optional[str] = None,
                        topology: Optional[str] = None) -> List[Dict]:
        """Calculate distances between consecutive heme pairs"""
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        if sequence_file:
            residue_pairs = self.read_linearized_sequence(sequence_file)
        else:
            if not topology:
                raise ValueError("Either sequence_file or topology must be specified")
            sequence = self.detect_sequence(atoms_dict, topology)
            if topology == 'linear':
                residue_pairs = list(zip(sequence[:-1], sequence[1:]))
            else:
                residue_pairs = []
                for branch in sequence:
                    residue_pairs.extend(list(zip(branch[:-1], branch[1:])))
        
        distances = []
        for res1, res2 in residue_pairs:
            if res1 in atoms_dict and res2 in atoms_dict:
                dist = self.calculate_min_distance(atoms_dict[res1], atoms_dict[res2])
                distances.append({'pair': [res1, res2], 'distance': dist})
            
        return distances

    def measure_subunit_length(self, pdb_file: str,
                             sequence_file: Optional[str] = None,
                             topology: Optional[str] = None) -> float:
        """Calculate distance between terminal FE atoms"""
        atoms_dict = self.read_pdb_atoms(pdb_file)
        
        if sequence_file:
            residue_pairs = self.read_linearized_sequence(sequence_file)
            first_res = residue_pairs[0][0]
            second_to_last_res = residue_pairs[-1][0]
        else:
            if not topology or topology != 'linear':
                raise ValueError("Subunit length measurement requires linear topology")
            sequence = self.detect_linear_sequence(atoms_dict)
            first_res = sequence[0]
            second_to_last_res = sequence[-1]
        
        if 'FE' in atoms_dict.get(first_res, {}) and 'FE' in atoms_dict.get(second_to_last_res, {}):
            return np.linalg.norm(
                atoms_dict[first_res]['FE'] - atoms_dict[second_to_last_res]['FE']
            )
        else:
            raise ValueError(f"Missing FE atoms in residues {first_res} or {second_to_last_res}")

    def classify_stacking(self, angle: float) -> str:
        """Classify heme stacking based on angle."""
        if angle < 45.0 or angle > 135.0:
            return "slip-stacked"
        else:
            return "T-stacked"

    def analyze_stacking_statistics(self, angles: List[Dict], distances: List[Dict]) -> Dict:
        """Analyze statistics for different stacking arrangements."""
        # Combine angle and distance data
        stacking_data = {
            'slip-stacked': {'angles': [], 'distances': []},
            'T-stacked': {'angles': [], 'distances': []}
        }
        
        # Match angles with distances and classify
        for angle_data, dist_data in zip(angles, distances):
            if angle_data['pair'] != dist_data['pair']:
                continue  # Skip if pairs don't match
                
            stacking_type = self.classify_stacking(angle_data['angle'])
            stacking_data[stacking_type]['angles'].append(angle_data['angle'])
            stacking_data[stacking_type]['distances'].append(dist_data['distance'])
        
        # Calculate statistics
        stats = {}
        for stype in ['slip-stacked', 'T-stacked']:
            if stacking_data[stype]['angles']:  # Only if we have data
                stats[stype] = {
                    'count': len(stacking_data[stype]['angles']),
                    'mean_angle': np.mean(stacking_data[stype]['angles']),
                    'std_angle': np.std(stacking_data[stype]['angles']),
                    'mean_distance': np.mean(stacking_data[stype]['distances']),
                    'std_distance': np.std(stacking_data[stype]['distances']),
                    'min_distance': np.min(stacking_data[stype]['distances']),
                    'max_distance': np.max(stacking_data[stype]['distances'])
                }
        
        return stats

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process PDB files to analyze heme groups')
    parser.add_argument('pdb_file', type=str, help='Input PDB file')
    parser.add_argument('--sequence-file', '-s', type=str, 
                       help='Input sequence file (optional)')
    parser.add_argument('--topology', '-t', type=str,
                       choices=['linear', 'branched'],
                       help='Topology type for automatic sequence detection (required if no sequence file)')
    parser.add_argument('--output', '-o', type=str, help='Output JSON file (optional)')
    parser.add_argument('--measurements', '-m', type=str, nargs='+',
                       choices=['angles', 'distances', 'length'],
                       default=['angles'],
                       help='Measurements to compute (default: angles)')
    parser.add_argument('--distance-cutoff', '-d', type=float, default=20.0,
                       help='Distance cutoff for neighboring hemes in Angstroms (default: 20.0)')
    return parser.parse_args()

def format_measurement_output(results: Dict) -> str:
    """Format measurement results as human-readable text."""
    output = []

    # Format detected sequence if present
    if 'detected_sequence' in results:
        output.append("Detected Linear Sequence:")
        output.append("-" * 22)
        output.append(" → ".join(str(x) for x in results['detected_sequence']))
        output.append("")

    if 'detected_branches' in results:
        output.append("Detected Branched Structure:")
        output.append("-" * 24)
        for i, branch in enumerate(results['detected_branches'], 1):
            output.append(f"Branch {i}: " + " → ".join(str(x) for x in branch))
        output.append("")

    # Format combined angle and distance data with stacking classification
    if 'angles' in results and 'distances' in results:
        output.append("Heme Pair Analysis:")
        output.append("-" * 17)
        for angle_data, dist_data in zip(results['angles'], results['distances']):
            res1, res2 = angle_data['pair']
            angle = angle_data['angle']
            dist = dist_data['distance']
            stacking = "Slip-stacked" if angle < 45.0 or angle > 135.0 else "T-stacked"
            output.append(f"Hemes {res1:4d}-{res2:<4d}: {angle:6.1f}° | {dist:6.1f} Å | {stacking}")
        output.append("")

    # Format stacking statistics if present
    if 'stacking_stats' in results:
        stats = results['stacking_stats']
        output.append("Stacking Statistics:")
        output.append("-" * 19)
        for stype in ['slip-stacked', 'T-stacked']:
            if stype in stats:
                s = stats[stype]
                output.append(f"{stype.title()} Arrangements:")
                output.append(f"  Count: {s['count']}")
                output.append(f"  Angles: {s['mean_angle']:.1f}° ± {s['std_angle']:.1f}°")
                output.append(f"  Distances: {s['mean_distance']:.1f} ± {s['std_distance']:.1f} Å")
                output.append(f"  Distance range: {s['min_distance']:.1f} - {s['max_distance']:.1f} Å")
                output.append("")

    # Format length if present
    if 'length' in results:
        output.append("Subunit Length:")
        output.append("-" * 14)
        output.append(f"Terminal Fe-Fe distance: {results['length']:.1f} Å")
        output.append("")

    return "\n".join(output)

def main():
    """Main function for command-line execution."""
    args = parse_args()

    if not Path(args.pdb_file).exists():
        print(f"Error: PDB file '{args.pdb_file}' not found", file=sys.stderr)
        sys.exit(1)

    if args.sequence_file and not Path(args.sequence_file).exists():
        print(f"Error: Sequence file '{args.sequence_file}' not found", file=sys.stderr)
        sys.exit(1)

    if not args.sequence_file and not args.topology:
        print("Error: Either --sequence-file or --topology must be specified", file=sys.stderr)
        sys.exit(1)

    processor = PDBProcessor()
    processor.distance_cutoff = args.distance_cutoff
    results = {}

    try:
        # If using automatic detection, include the detected sequence
        if not args.sequence_file:
            atoms_dict = processor.read_pdb_atoms(args.pdb_file)
            detected_sequence = processor.detect_sequence(atoms_dict, args.topology)
            if args.topology == 'linear':
                results['detected_sequence'] = detected_sequence
            else:
                results['detected_branches'] = detected_sequence

        # Perform requested measurements
        if 'angles' in args.measurements:
            angles = processor.measure_heme_angles(
                args.pdb_file,
                sequence_file=args.sequence_file,
                topology=args.topology
            )
            results['angles'] = angles

        if 'distances' in args.measurements:
            distances = processor.measure_avg_dist(
                args.pdb_file,
                sequence_file=args.sequence_file,
                topology=args.topology
            )
            results['distances'] = distances

            # Add stacking analysis if we have both angles and distances
            if 'angles' in results and 'distances' in results:
                results['stacking_stats'] = processor.analyze_stacking_statistics(
                    results['angles'], results['distances']
                )

        if 'length' in args.measurements:
            if args.topology != 'branched':
                length = processor.measure_subunit_length(
                    args.pdb_file,
                    sequence_file=args.sequence_file,
                    topology=args.topology
                )
                results['length'] = length
            else:
                print("Warning: Subunit length measurement not applicable for branched topology",
                      file=sys.stderr)

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

    # Format and output results
    formatted_output = format_measurement_output(results)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(formatted_output)
    else:
        print(formatted_output)

if __name__ == '__main__':
    main()
