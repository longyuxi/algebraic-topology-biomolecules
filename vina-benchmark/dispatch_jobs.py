# Making the csv file

import glob
import os
import pandas as pd
import pathlib
import socket

cwd = pathlib.Path(__file__).parent.resolve()
CSV_FILE = f'{cwd}/jobs.csv'
NUM_FOLDERS_TO_RUN = 3
SBATCH_TEMPLATE = f"""#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yl708@duke.edu     # Where to send mail
#SBATCH --partition=scavenger
#SBATCH --exclusive
#SBATCH --time='7-0'
#SBATCH --chdir='/work/yl708/algebraic-topology-biomolecules/vina-benchmark'
#SBATCH --mem=4g

source ~/.bashrc
source ~/.bash_profile

cd /work/yl708/algebraic-topology-biomolecules/vina-benchmark
conda activate vina

date
hostname
"""

def main():
    df = get_df(CSV_FILE)
    # Dispatcher
    import numpy as np
    import random

    # randomly sample an entry that has not been attempted until all has been
    # attempted or reached NUM_FOLDERS_TO_RUN
    i = 0
    while True:
        if i == NUM_FOLDERS_TO_RUN:
            break
        i += 1

        valid_entries = df.loc[df['attempted'] == False]
        if len(valid_entries) == 0:
            break

        idx = random.choice(valid_entries.index)
        df.at[idx, 'attempted'] = True
        df.to_csv(CSV_FILE)

        # sbatch run job wrapper
        sbatch_cmd = SBATCH_TEMPLATE + f'\npython {pathlib.Path(__file__).parent + "dispatch_jobs.py"}'
        print(sbatch_cmd)

def get_df(csv_file) -> pd.DataFrame:
    if not os.path.exists(csv_file):
        df = pd.DataFrame(columns=['name', 'folder', 'attempted', 'finished'])

        if socket.gethostname() == '1080-ubuntu':
            folders = glob.glob('/home/longyuxi/Documents/mount/scPDB/*')
        else:
            folders = glob.glob('/work/yl708/scPDB/*')
        names = [f.split('/')[-1] for f in folders]

        df['name'] = names
        df['folder'] = folders
        df['attempted'] = False
        df['finished'] = False
        df['error'] = False

        df.to_csv(csv_file)
    else:
        df = pd.read_csv(csv_file)

if __name__ == '__main__':
    main()