from os.path import join
from pathlib import Path
from glob import glob

import numpy as np
from joblib import Parallel, delayed
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

def extract_coords_and_seq(pdb_path, output_dir):
    """Extract the alpha carbon atomic coordinates and primary sequences in a PDB file and save to
    disk.

    Note: This function expects 'coordinates' and 'primary_seqs' subdirectories to exist within the
    output directory.

    Args:
        pdb_path (str or path): Path to the PDB file
        output_dir (str or path): Path to the output directory
    """
    # Parse the PDB ID from the PDB file name
    pdb_id = Path(pdb_path).stem

    # Load the PDB file
    # 'pdb-seqres' extracts the SEQRES values in the PDB file
    # Source: https://biopython.org/docs/1.75/api/Bio.SeqIO.PdbIO.html
    seq_io = SeqIO.parse(pdb_path, 'pdb-seqres')

    # How to read PDB files with Biopython:
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec182
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_path)

    all_coords = {}
    seqs = []
    # Note: Biopython considers a PDB file to contain the following as nested objects:
    # structure -> models -> chains -> residues -> atoms
    # Source: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec193
    #
    # Additional examples for navigating through Biopython's Structure object:
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec210
    for model in structure.get_list():
        for chain, record in zip(model.get_list(), seq_io):
            chain_coords = []

            for residue in chain.get_list():
                try:
                    # Attempt to load coordinates of the alpha carbon in each residue
                    chain_coords.append(residue['CA'].get_coord())
                except KeyError:
                    # No carbon alpha's were present in the amino acid
                    pass

            # Only store atomic coordinates and primary sequence if there are alpha carbons in the
            # chain
            if chain_coords:
                all_coords[str(model.get_id())] = {str(chain.get_id()):chain_coords}
                seqs.append(record)

                # The 'dbxrefs' variable appears to contain cross-reference values that were
                # extracted from the PDB file. UNP stands for Uniprot.
                print('Database cross-reference values:', record.dbxrefs)

    # Atomic coords. are stored as:
    # {model ID: {chain ID:np.array([atom, (x,y,z)])}}
    np.savez(join(output_dir, 'coordinates', f'{pdb_id}_c_alphas.npz'), **all_coords, )

    # Primary sequences are stored in a single FASTA file
    SeqIO.write(seqs, join(output_dir, 'primary_seqs', f'{pdb_id}.fasta'), 'fasta')

    print(f'Saved {pdb_id}\'s c_alpha coordinates and primary sequence to file')

# This is the path to the unprocessed PDB files
PDB_DIR = join('Project', 'PDBs')

# Get the paths to PDB files in the directory
pdb_paths = glob(join(PDB_DIR, '*.pdb'))

output_dir = join('Project', 'outputs')

# Use this flag to process PDB files:
#   True: in serial -> easier debugging
#   False: in parallel -> faster runtime if multiple cores are available
TESTING = False

if TESTING:
    for pdb_path in pdb_paths:
        extract_coords_and_seq(pdb_path, output_dir)
else:
    # Set the number of processes to use
    n_proc = min(len(pdb_paths), 20)

    # Background and examples of embarrassingly parallel (AKA pleasantly parallel) for loops with
    # Joblib: https://joblib.readthedocs.io/en/latest/parallel.html#common-usage
    Parallel(n_jobs=n_proc)(delayed(extract_coords_and_seq)(pdb_path, output_dir)\
        for pdb_path in pdb_paths)