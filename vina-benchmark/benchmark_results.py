# %%
import pandas as pd
from dispatch_jobs import CSV_FILE

df = pd.read_csv(CSV_FILE)

finished_entries = df.loc[(df['attempted'] == True) & (df['finished'] == True) & (df['error'] == False)]

# %%
import matplotlib.pyplot as plt

# %%
from utils import calculate_bbox
import pybel

# Out of these far mismatches, how many are because the search space doesn't contain it?
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

# %%
# Collecting RMSD from each successful dock

rmsds = {}
for i in finished_entries.index:
    entry_name = finished_entries.at[i, 'name']
    entry_folder = finished_entries.at[i, 'folder']
    l = entry_name
    if not bbox_contains(ground_truth_bboxes[l], search_bboxes[l]):
        continue

    rmsd_file = entry_folder + '/vina_prediction_rmsds.txt'
    with open(rmsd_file) as f:
        r = f.read().splitlines()

    r = [float(e) for e in r]
    rmsds[entry_name] = r


plt.rcParams["figure.figsize"] = (10, 20)
fig, axs = plt.subplots(3, 1)
axs = iter(axs.flatten())

# Distribution of RMSD first try
ax = next(axs)
print(ax)
first_try_rmsds = []
for key in rmsds:
    first_try_rmsds.append(rmsds[key][0])

first_try_distrib, _, _ = ax.hist(first_try_rmsds, bins=range(30), edgecolor="black")
ax.set_title(f'Distribution of RMSD first try, n={len(rmsds)}')

# Distribution of min RMSD among first three tries
ax = next(axs)
first_three_tries_rmsds = []
for key in rmsds:
    first_three_tries_rmsds.append(min(rmsds[key][:3]))

ax.set_title('Distribution of min RMSD in three tries')
first_three_tries_distrib, _, _ = ax.hist(first_three_tries_rmsds, bins=range(30), edgecolor="black")

# Distribution of min RMSD among all tries
ax = next(axs)
all_tries_rmsds = []
for key in rmsds:
    all_tries_rmsds.append(min(rmsds[key]))

ax.set_title('Distribution of min RMSD in all tries')
all_tries_distrib, _, _ = ax.hist(all_tries_rmsds, bins=range(30), edgecolor="black")

plt.savefig('results.jpg')

