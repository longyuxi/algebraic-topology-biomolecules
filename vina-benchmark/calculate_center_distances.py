# %%
import pandas as pd
from dispatch_jobs import CSV_FILE
import pybel
from utils import calculate_center_distance

df = pd.read_csv(CSV_FILE)

finished_entries = df.loc[(df['finished'] == True) & (df['error'] == False)]

mol_dists = {}
# Center of binding site rmsd
for i in finished_entries.index:
    entry_name = finished_entries.at[i, 'name']
    entry_folder = finished_entries.at[i, 'folder']

    # Get ground truth molecule
    gt = list(pybel.readfile('pdbqt', entry_folder + '/ligand.pdbqt'))
    assert len(gt) == 1
    gt = gt[0]

    curr_dists = []
    for prediction in pybel.readfile('pdbqt', entry_folder + '/ligand_vina_out.pdbqt'):
        curr_dists.append(calculate_center_distance(gt, prediction))

    mol_dists[entry_name] = curr_dists

# %%
min_mol_dists = {}
for key in mol_dists:
    min_mol_dists[key] = min(mol_dists[key])

# %%
import matplotlib.pyplot as plt
import numpy as np

plt.hist(min_mol_dists.values(), bins=np.arange(0, 20, 0.5))
plt.title('COM distances')
plt.savefig('com_dist.jpg')


# %%



