import glob
import os
import pandas as pd
import pathlib
import logging
import redis
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from tqdm import tqdm
import time

logging.basicConfig(level=logging.INFO)

###############################
# Platform specific variables #
#                             #
# Change to fit to job        #
###############################


KEY_PREFIX = 'perturb_' # Prefix of every job, as appears in the Redis database
CLUSTER = 'CS' # or 'DCC'

if CLUSTER == 'CS':
    cwd = pathlib.Path(__file__).parent.resolve()
    NUM_JOBS_TO_SUBMIT = 10
    PYTHON_EXECUTABLE = '/usr/project/dlab/Users/jaden/mambaforge/envs/tnet2017/bin/python'
    ROOT_DIR = '/usr/project/dlab/Users/jaden/algebraic-topology-biomolecules/perturbations'
    os.system(f'mkdir -p {ROOT_DIR}/slurm-outs')
    SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm-outs/%x-%j-slurm.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci,grisman

source ~/.zshrc
date
hostname
conda activate tnet2017
cd {ROOT_DIR}


    """

    DB = redis.Redis(host='cybermen', port=6379, decode_responses=True, password="topology")
    FOLDERS = glob.glob('/usr/project/dlab/Users/jaden/pdbbind/refined-set/*')
    PDB_NAMES = [f.split('/')[-1] for f in FOLDERS]

    ORIGINAL_PROTEIN_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_protein.pdb'
    LIGAND_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_ligand.mol2'
    PERTURBATION_SAVE_FOLDER = '/hpc/group/donald/yl708/perturbations-data/'

elif CLUSTER == 'DCC':
    raise NotImplementedError
    NUM_JOBS_TO_SUBMIT = 81000
    PYTHON_EXECUTABLE = '/hpc/group/donald/yl708/mambaforge/envs/tnet2017/bin/python'
    ROOT_DIR = '/hpc/group/donald/yl708/TopologyNet-2017/perturbations'
    SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --partition=common-old,scavenger
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm_logs/%x-%j-slurm.out
#SBATCH --mem=2500M

source ~/.bashrc
source ~/.bash_profile
date
hostname
conda activate tnet2017
cd {ROOT_DIR}


    """

    DB = redis.Redis(host='dcc-login-03', port=6379, decode_responses=True, password="topology")
    ORIGINAL_PROTEIN_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_protein.pdb'
    LIGAND_FILE = '/hpc/group/donald/yl708/perturbations-data/1a4k_ligand.mol2'
    PERTURBATION_SAVE_FOLDER = '/hpc/group/donald/yl708/perturbations-data/'

else:
    raise Exception     # Incorrect specification of cluster variable


# Perturbations
def simple_sampler(original_position, epsilon, idx):
    # Translates the original position in all six directions
    if idx >= 6: raise Exception
    new_position = original_position
    new_position[int(idx / 2)] += (idx%2) * epsilon
    return new_position


# PERTURBATION_SAMPLER returns a position based on original_position, idx
PERTURBATION_SAMPLER = lambda original_position, idx: simple_sampler(original_position, 0.1, idx)
PERTURBATIONS_PER_ATOM = 6


#############################
# Pre-execution Tests       #
#############################

# Database connection
DB.set('connection-test', '123')
if DB.get('connection-test') == '123':
    DB.delete('abc')
    logging.info('Database connection successful')
else:
    raise Exception     # Database connection failed


pdbio = PDBIO()
pdbparser = PDBParser(PERMISSIVE=1, QUIET=True)


get_num_atoms = lambda structure: len([a for a in structure.get_atoms()])
get_pdb_name = lambda file: file.split('/')[-1].split('_')[0]

protein_structure = pdbparser.get_structure(ORIGINAL_PROTEIN_FILE, ORIGINAL_PROTEIN_FILE)
n_atoms = get_num_atoms(protein_structure)

print(f'Number of atoms in protein {get_pdb_name(ORIGINAL_PROTEIN_FILE)}: {n_atoms}')

KEY_PREFIX += get_pdb_name(ORIGINAL_PROTEIN_FILE) + '_'

#############################
# Actual logic              #
#############################

# Each entry in the redis database should be a dictionary in the following form

# job_index (key_prefix + job_index = key in database),
# {
#     protein_file: name of pdb file,
#     ligand_file: name of ligand mol2 file,
#     save_folder: where to save the output
#     atom_index: index of atom that is perturbed,
#     perturbation_index: perturbation index (passed in to the perturbation function)
#     attempted: true/false
#     error: true/false
#     finished: true/false
# }

def main(dry_run=False):
    # Initialize database on first run
    if dry_run:
        populate_db()

    # Then submit jobs until either running out of entries or running out of number of jobs to submit
    i = 0

    database_keys = DB.keys(KEY_PREFIX + '*')
    for key in database_keys:
        if i == NUM_JOBS_TO_SUBMIT:
            break
        info = DB.hgetall(key)

        if info['finished'] == 'True' and info['error'] == 'False':
        # if info['attempted'] == 'True':
            continue
        else:
            i += 1
            # submit job for it
            if not dry_run:
                info['attempted'] = 'True'
                DB.hset(key, mapping=info)

                # sbatch run job wrapper
                sbatch_cmd = SBATCH_TEMPLATE + f'\n{PYTHON_EXECUTABLE} {str(pathlib.Path(__file__).parent) + "/job_wrapper.py"} --key {key}'

                # print(sbatch_cmd)
                with open('run.sh', 'w') as f:
                    f.write(sbatch_cmd)

                os.system(f'sbatch --job-name={key} run.sh')

    if dry_run:
        print(f'Number of jobs that would be submitted: {i}')
        time.sleep(5)
    else:
        print(f'Number of jobs submitted: {i}')


def populate_db():
    logging.info('Populating database')
    n_jobs = n_atoms * PERTURBATIONS_PER_ATOM
    keys = [KEY_PREFIX + str(i) for i in range(n_jobs)]

    database_keys = DB.keys()
    for k in tqdm(keys):
        if k in database_keys:
            logging.debug(f"Key {k} already exists in database")
            continue

        # First perturb the original pdb file and save it to a folder with corresponding index
        # Then add the entry to the database
        idx = int(k.split(KEY_PREFIX)[1])
        atom_idx = int(idx / PERTURBATIONS_PER_ATOM)
        perturbation_idx = idx % PERTURBATIONS_PER_ATOM

        # Perturb the original pdb file
        protein_structure = pdbparser.get_structure(ORIGINAL_PROTEIN_FILE, ORIGINAL_PROTEIN_FILE)
        atom = [a for a in protein_structure.get_atoms()][atom_idx]
        atom.set_coord(PERTURBATION_SAMPLER(atom.get_coord(), perturbation_idx))

        # Save the perturbed pdb file
        save_folder = pathlib.Path(f'{PERTURBATION_SAVE_FOLDER}/{KEY_PREFIX}{idx}').mkdir(parents=True, exist_ok=True)
        pdbio.set_structure(protein_structure)
        pdbio.save(f'{PERTURBATION_SAVE_FOLDER}/{KEY_PREFIX}{idx}/{KEY_PREFIX}{idx}.pdb')

        # Add the entry to the database
        DB.hset(k, mapping={
            'protein_file': f'{PERTURBATION_SAVE_FOLDER}/{KEY_PREFIX}{idx}/{KEY_PREFIX}{idx}.pdb',
            'ligand_file': LIGAND_FILE,
            'save_folder': f'{PERTURBATION_SAVE_FOLDER}/{KEY_PREFIX}{idx}',
            'atom_index': atom_idx,
            'perturbation_index': perturbation_idx,
            'attempted': 'False',
            'error': 'False',
            'finished': 'False'
        })



def rebuild_db():
    raise NotImplementedError

def get_db():
    # Pinnacle of OOP
    return DB

if __name__ == '__main__':
    # rebuild_db()
    main(dry_run=True)
    main(dry_run=False)
