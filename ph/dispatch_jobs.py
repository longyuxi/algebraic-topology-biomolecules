# -*- coding: future_fstrings -*-
# Making the csv file
# If running on python2, also need to `pip install future-fstrings pathlib`

import glob
import os
import pandas as pd
import pathlib
import socket

# Change these as necessary
cwd = pathlib.Path(__file__).parent.resolve()
CSV_FILE = f'{cwd}/jobs.csv'
NUM_JOBS_TO_SUBMIT = 5000
PYTHON_EXECUTABLE = '/work/yl708/algebraic-topology-biomolecules/algtop-environment/bin/python'
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

def main():
    # This function shouldn't need to be changed

    df = get_df(CSV_FILE)
    # Dispatcher
    import numpy as np
    import random

    # randomly sample an entry that has not been attempted until all has been
    # attempted or reached NUM_JOBS_TO_SUBMIT
    i = 0
    while True:
        if i == NUM_JOBS_TO_SUBMIT:
            break
        i += 1

        valid_entries = df.loc[((df['attempted'] == True) | df['attempted'].isnull()) & df['error'] == False]

        if len(valid_entries) == 0:
            break

        idx = np.random.choice(valid_entries.index)
        df.loc[idx, 'attempted'] = True
        df.to_csv(CSV_FILE, index=False)

        # sbatch run job wrapper
        sbatch_cmd = SBATCH_TEMPLATE + f'\n{PYTHON_EXECUTABLE} {str(pathlib.Path(__file__).parent) + "/job_wrapper.py"} --csv {CSV_FILE} --idx {idx}'

        # print(sbatch_cmd)
        with open('run.sh', 'w') as f:
            f.write(sbatch_cmd)

        os.system(f'sbatch --job-name={idx} run.sh')
        # os.system('rm run.sh')


def get_df(csv_file):
    if not os.path.exists(csv_file):
        df = pd.DataFrame()

        # Keep these three columns - they are needed for the job script
        df['attempted'] = [False] * len(df)
        df['finished'] = [False] * len(df)
        df['error'] = [False] * len(df)

        # Customize all the other columns as fit for the job
        if socket.gethostname() == '1080-ubuntu':
            folders = glob.glob('/home/longyuxi/Documents/mount/scPDB/*')
        else:
            folders = glob.glob('/work/yl708/pdbbind/refined-set/*')
        names = [f.split('/')[-1] for f in folders]
        df['name'] = names
        df['folder'] = folders

        df.to_csv(csv_file, index=False)
    else:
        df = pd.read_csv(csv_file)

    return df

if __name__ == '__main__':
    main()
