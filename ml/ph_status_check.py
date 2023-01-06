"""Standalone script that generates a ph_status.csv, which contains the status of persistent homology run on each folder
"""

# %%
import glob
import os
import pandas as pd
import pathlib
import socket

def load_pdbbind_data_index(index_filename: str) -> pd.DataFrame:
    index = pd.read_csv(index_filename, delim_whitespace=True, skiprows=6, names=['PDB code', "resolution", "release year", "-logKd/Ki", "Kd/Ki", "slashes", "reference", "ligand name"])

    index.drop(columns='slashes', inplace=True)
    index['ligand name'] = index.apply(lambda row:  row['ligand name'][1:][:-1], axis=1)

    return index

def get_df(csv_file) -> pd.DataFrame:
    """Creates a dataframe for storing persistent homology status
    """
    if not os.path.exists(csv_file):
        df = pd.DataFrame()

        # Customize all the other columns as fit for the job
        if socket.gethostname() == '1080-ubuntu':
            folders = glob.glob('/home/longyuxi/Documents/mount/pdbbind-dataset/refined-set/*')
        else:
            raise NotImplementedError
        names = [f.split('/')[-1] for f in folders]
        df['name'] = names
        df['folder'] = folders

        # Keep these three columns - they are needed for the job script
        df['attempted'] = False
        df['finished'] = False
        df['error'] = False

        df.to_csv(csv_file, index=False)
    else:
        df = pd.read_csv(csv_file)

    return df

# %%
# Iterate through all the folders to update the pd dataframe

def main():
    if socket.gethostname() == '1080-ubuntu':
        index_location = '/home/longyuxi/Documents/mount/pdbbind-dataset/index/INDEX_refined_data.2020'
    else:
        raise NotImplementedError
    pdbbind_df = load_pdbbind_data_index(index_location)

    CSV_FILE = 'ph_status.csv'
    df = get_df(CSV_FILE)

    for i in range(len(df)):
        binding_affinity = pdbbind_df.loc[pdbbind_df['PDB code'] == df.at[i, 'name']].iloc[0]['-logKd/Ki']
        df.at[i, '-logKd/Ki'] = float(binding_affinity)
        numfiles = len(glob.glob(df.at[i, 'folder'] + '/*'))
        if numfiles == 20:
            df.at[i, 'attempted'] = True
            df.at[i, 'finished'] = True
            df.at[i, 'error'] = False
        elif numfiles == 4:
            df.at[i, 'attempted'] = False
            df.at[i, 'finished'] = False
            df.at[i, 'error'] = False
        else:
            df.at[i, 'attempted'] = True
            df.at[i, 'finished'] = True
            df.at[i, 'error'] = True


    df.to_csv(CSV_FILE, index=False)

    finished_entries = df.loc[(df['attempted'] == True) & (df['finished'] == True) & (df['error'] == False)]

    print('num folders finished:', len(finished_entries))

if __name__ == '__main__':
    main()