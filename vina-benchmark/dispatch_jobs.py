# Making the csv file

import glob
import os
import pandas as pd

CSV_FILE = '/home/longyuxi/Documents/algebraic-topology-biomolecules/vina-benchmark/jobs.csv'

if not os.path.exists(CSV_FILE):
    df = pd.DataFrame(columns=['name', 'folder', 'attempted', 'finished'])

    folders = glob.glob('/home/longyuxi/Documents/mount/scPDB/*')
    names = [f.split('/')[-1] for f in folders]

    df['name'] = names
    df['folder'] = folders
    df['attempted'] = False
    df['finished'] = False
    df['error'] = False

    df.to_csv(CSV_FILE)
else:
    df = pd.read_csv(CSV_FILE)


# Dispatcher
from vinascpdb import run_on_folder
import numpy as np
import random

NUM_FOLDERS_TO_RUN = 100

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

    folder = df.at[idx, 'folder']

    df.at[idx, 'attempted'] = True
    df.to_csv(CSV_FILE)

    try:
        run_on_folder(folder)
    except:
        df.at[idx, 'error'] = True
        df.to_csv(CSV_FILE)

    df.at[idx, 'finished'] = True
    df.to_csv(CSV_FILE)
