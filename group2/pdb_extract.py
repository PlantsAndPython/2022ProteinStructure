from os.path import join
from pathlib import Path

import numpy as np
from glob import glob
from joblib import Parallel, delayed
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

def extract_coords_and_seq(pdb_path, output_dir):
    pdb_id = Path(pdb_path).stem

    seq_io = SeqIO.parse(pdb_path, 'pdb-seqres')
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_path)

    all_coords = {}
    seqs = []
    for model in structure.get_list():
        for chain, record in zip(model.get_list(), seq_io):
            chain_coords = []

            for residue in chain.get_list():
                try:
                    chain_coords.append(residue['CA'].get_coord())
                except KeyError:
                    # No carbon alpha's were present in the amino acid
                    pass

            if chain_coords:
                all_coords[str(model.get_id())] = {str(chain.get_id()):chain_coords}
                seqs.append(record)

    np.savez(join(output_dir, 'coordinates', f'{pdb_id}_c_alphas.npz'), **all_coords, )
    SeqIO.write(seqs, join(output_dir, 'primary_seqs', f'{pdb_id}.fasta'), 'fasta')

    print(f'Saved {pdb_id}\'s c_alpha coordinates and primary sequence to file')

pdb_dir = join('Project', 'PDBs')
pdb_paths = glob(join(pdb_dir, '*.pdb'))
output_dir = join('Project', 'outputs')

TESTING = True

if TESTING:
    for i in range(20):
        for pdb_path in pdb_paths:
            extract_coords_and_seq(pdb_path, output_dir)
else:
    n_proc = min(len(pdb_paths), 20)
    Parallel(n_jobs=n_proc)(delayed(extract_coords_and_seq)(pdb_path, output_dir) for pdb_path in pdb_paths)