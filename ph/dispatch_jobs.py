# Making the csv file
# If running on python2, also need to `pip install future-fstrings pathlib`

import glob
import os
import pandas as pd
import pathlib
import socket
import redis
import numpy as np
import random

# Change these as necessary

KEY_PREFIX = 'ph_' # Prefix of every job, as appears in the Redis database
CLUSTER = 'CS' # or 'DCC'

if CLUSTER == 'CS':
    cwd = pathlib.Path(__file__).parent.resolve()
    CSV_FILE = f'{cwd}/jobs.csv'
    NUM_JOBS_TO_SUBMIT = 5000
    PYTHON_EXECUTABLE = '/usr/project/dlab/Users/jaden/mambaforge/envs/algtop-ph/bin/python'
    ROOT_DIR = '/usr/project/dlab/Users/jaden/algebraic-topology-biomolecules/ph'
    os.system(f'mkdir -p {ROOT_DIR}/slurm-outs')
    SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --time='4:00:00'
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm-outs/%x-%j-slurm.out
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compsci

source ~/.zshrc
date
hostname
conda activate algtop-ph
cd {ROOT_DIR}


    """

    DB = redis.Redis(host='login-01', port=6379, decode_responses=True, password="topology")
    FOLDERS = glob.glob('/usr/project/dlab/Users/jaden/pdbbind/refined-set/*')
    PDB_NAMES = [f.split('/')[-1] for f in FOLDERS]

elif CLUSTER == 'DCC':
    cwd = pathlib.Path(__file__).parent.resolve()
    CSV_FILE = f'{cwd}/jobs.csv'
    NUM_JOBS_TO_SUBMIT = 5
    PYTHON_EXECUTABLE = '/usr/project/dlab/Users/jaden/mambaforge/envs/algtop-ph/bin/python'
    ROOT_DIR = '/work/yl708/algebraic-topology-biomolecules/ph'
    SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --partition=scavenger
#SBATCH --time='1:00:00'
#SBATCH --requeue
#SBATCH --chdir={ROOT_DIR}
#SBATCH --output={ROOT_DIR}/slurm-outs/%x-%j-slurm.out
#SBATCH --mem=4G

source ~/.bashrc
source ~/.bash_profile
date
hostname
conda activate /work/yl708/algebraic-topology-biomolecules/algtop-environment
cd {ROOT_DIR}

module load R/4.1.1-rhel8
module load Matlab

    """

    FOLDERS = glob.glob('/work/yl708/pdbbind/refined-set/*')

else:
    raise Exception # incorrect specification of cluster

get_name = lambda folder: folder.split('/')[-1]
get_db_key = lambda folder: KEY_PREFIX + get_name(folder)

def main():
    # Initialize database
    populate_db()
    folder_it = iter(FOLDERS)

    # Then submit jobs until either running out of entries or running out of number of jobs to submit
    i = 0
    for folder in folder_it:
        if i == NUM_JOBS_TO_SUBMIT:
            break
        info = DB.hgetall(get_db_key(folder))

        if info['finished'] == 'True' and info['error'] == 'False':
            continue
        else:
            # submit job for it
            info['attempted'] = 'True'
            DB.hset(get_db_key(folder), mapping=info)

            # sbatch run job wrapper
            sbatch_cmd = SBATCH_TEMPLATE + f'\n{PYTHON_EXECUTABLE} {str(pathlib.Path(__file__).parent) + "/job_wrapper.py"} --key {get_db_key(folder)}'

            # print(sbatch_cmd)
            with open('run.sh', 'w') as f:
                f.write(sbatch_cmd)

            os.system(f'sbatch --job-name={get_name(folder)} run.sh')

def populate_db():
    # For each key, add its default entry to the database if it does not exist
    list(map(lambda folder: DB.hset(get_db_key(folder), mapping={
        'name': get_name(folder),
        'folder': folder,
        'attempted': 'False',
        'finished': 'False',
        'error': 'False'
    }), filter(lambda folder: not DB.exists(get_db_key(folder)), FOLDERS)))

# def remake_df(csv_file):
#     # Remakes df based on the folder
#     print('Remaking csv file')
#     os.system(f'rm {csv_file}')
#     df = get_df(csv_file)
#     n_fin = 0
#     n_err = 0
#     for idx in range(len(df)):
#         folder_name = df.at[idx, 'folder']
#         if len(glob.glob(folder_name + '/*')) == 6:
#             df.at[idx, 'attempted'] = False
#             df.at[idx, 'finished'] = False
#             df.at[idx, 'error'] = False
#         elif len(glob.glob(folder_name + '/*')) == 20:
#             df.at[idx, 'attempted'] = True
#             df.at[idx, 'finished'] = True
#             df.at[idx, 'error'] = False
#             n_fin += 1
#         else:
#             df.at[idx, 'attempted'] = True
#             df.at[idx, 'finished'] = False
#             df.at[idx, 'error'] = True
#             n_err += 1

#     df.to_csv(csv_file, index=False)
#     print(f'Number finished: {n_fin}; number erred: {n_err}; number not started: {len(df) - n_fin - n_err}')


def get_db():
    # Pinnacle of OOP
    return DB

if __name__ == '__main__':
    main()
