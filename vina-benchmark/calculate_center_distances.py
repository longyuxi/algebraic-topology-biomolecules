# %%
import pandas as pd
from dispatch_jobs import CSV_FILE
import pybel
from utils import calculate_center_distance, calculate_bbox

df = pd.read_csv(CSV_FILE)

finished_entries = df.loc[(df['attempted'] == True) & (df['finished'] == True) & (df['error'] == False)]

search_bboxes = {}
ground_truth_bboxes = {}
for i in finished_entries.index:
    entry_name = finished_entries.at[i, 'name']
    entry_folder = finished_entries.at[i, 'folder']

    search_bbox_file = entry_folder + '/vina_box.txt'

    with open(search_bbox_file) as f:
        r = f.read().splitlines()

    r = [float(e.split('=')[-1]) for e in r]
    search_bboxes[entry_name] = r

    mols = list(pybel.readfile('pdbqt', entry_folder + '/ligand.pdbqt'))
    assert len(mols) == 1
    mol = mols[0]
    gt_bbox_center, gt_bbox_size =  calculate_bbox(mol)
    ground_truth_bboxes[entry_name] = [*gt_bbox_center, *gt_bbox_size]

# %%
def bbox_contains(ligand_bbox, search_space_bbox):
    out = True
    for dim in range(0, 3):
        if ligand_bbox[dim] + ligand_bbox[dim + 3] / 2 > search_space_bbox[dim] + search_space_bbox[dim + 3] / 2:
            out = False
        if ligand_bbox[dim] - ligand_bbox[dim + 3] / 2 < search_space_bbox[dim] - search_space_bbox[dim + 3] / 2:
            out = False
    return out

mol_dists = {}
# Center of binding site rmsd
for i in finished_entries.index:
    try:
        entry_name = finished_entries.at[i, 'name']
        entry_folder = finished_entries.at[i, 'folder']
        l = entry_name
        if not bbox_contains(ground_truth_bboxes[l], search_bboxes[l]):
            continue

        # Get ground truth molecule
        gt = list(pybel.readfile('pdbqt', entry_folder + '/ligand.pdbqt'))
        assert len(gt) == 1
        gt = gt[0]

        curr_dists = []
        for prediction in pybel.readfile('pdbqt', entry_folder + '/ligand_vina_out.pdbqt'):
            curr_dists.append(calculate_center_distance(gt, prediction))

        mol_dists[entry_name] = curr_dists
    except:
        pass

# %%
min_mol_dists = {}
for key in mol_dists:
    min_mol_dists[key] = min(mol_dists[key])

first_mol_dists = {}
for key in mol_dists:
    first_mol_dists[key] = mol_dists[key][0]

print(len(min_mol_dists))
# %%
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = (10, 20)
fig, axs = plt.subplots(2, 1)
axs = iter(axs.flatten())

# Distribution of RMSD first try
ax = next(axs)
ax.hist(first_mol_dists.values(), bins=np.arange(0, 20, 0.5))
ax.set_title(f'COM dist of first prediction. n={len(min_mol_dists)}')

ax = next(axs)
ax.hist(min_mol_dists.values(), bins=np.arange(0, 20, 0.5))
ax.set_title(f'Min COM distances among all.')
plt.savefig('com_dist.jpg')


# %%



