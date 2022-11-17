from os.path import join
from pathlib import Path

import mdtraj as mdj
import numpy as np
from glob import glob
from joblib import Parallel, delayed

def extract_coords_and_seq(pdb_path, output_dir):
    traj = mdj.load_pdb(pdb_path)
    top = traj.topology

    # Load carbon alpha atoms of each chain
    c_alpha_coords = {}
    residues = {}
    for chain in top.chains:
        # 'CA' indicates the carbon alpha atoms
        indices = top.select(f'chainid {chain.index} and name CA')

        # Ignore chains that have none
        if len(indices):
            # MDTraj considers a PDB file to be a trajectory with one frame
            coords = traj.xyz[0][indices]
            c_alpha_coords[str(chain.index)] = coords

        # Keep only the amino acid name (don't need the index)
        chain_residues = [r.name[:3] for r in chain.residues]
        residues[str(chain.index)] = chain_residues

    pdb_id = Path(pdb_path).stem
    np.savez(join(output_dir, 'coordinates', f'{pdb_id}_c_alphas.npz'), **c_alpha_coords)
    np.savez(join(output_dir, 'primary_seqs', f'{pdb_id}_residues.npz'), **residues)

    print(f'Saved {pdb_id}\'s c_alpha coordinates and primary sequence to file')

pdb_dir = join('Project', 'PDBs')
pdb_paths = glob(join(pdb_dir, '*.pdb'))
output_dir = join('Project', 'outputs')

n_proc = min(len(pdb_paths), 20)
Parallel(n_jobs=n_proc)(delayed(extract_coords_and_seq)(pdb_path, output_dir) for pdb_path in pdb_paths)