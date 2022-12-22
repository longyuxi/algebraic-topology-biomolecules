import os
from utils import calculate_bbox, calculate_rmsd
import pybel
import argparse

def run_on_folder(folder):
    os.chdir(folder)
    adfr_root = '/home/longyuxi/ADFRsuite-1.0/bin'

    # Prepare the receptor
    assert os.system(f'obabel protein.mol2 -O protein.pdb') == 0
    assert os.system(f'{adfr_root}/prepare_receptor -r protein.pdb -o protein.pdbqt') == 0

    # Prepare the ligand
    assert os.system('mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt') == 0

    # Calculate bounding box

    mols = list(pybel.readfile('mol2', 'protein.mol2'))
    assert len(mols) == 1
    mol = mols[0]
    bbox_center, bbox_lengths = calculate_bbox(mol)

    with open('vina_box.txt', 'w') as txt:
        txt.write(f'center_x = {bbox_center[0]}\ncenter_y = {bbox_center[1]}\ncenter_z = {bbox_center[2]}\nsize_x = 20.0\nsize_y = 20.0\nsize_z = 20.0')

    # Run Vina
    cmd = f'vina --receptor protein.pdbqt --ligand ligand.pdbqt \
        --config vina_box.txt \
        --exhaustiveness=32 --out ligand_vina_out.pdbqt'

    # cmd = f'vina --receptor protein.pdbqt --ligand ligand.pdbqt \
    #        --autobox \
    #        --exhaustiveness=32 --out ligand_vina_out.pdbqt'

    assert os.system(cmd) == 0

    # Compare against grouth truth
    gt = list(pybel.readfile('pdbqt', 'ligand.pdbqt'))
    assert len(gt) == 1
    gt = gt[0]

    # print(len(list(pybel.readfile('pdbqt', 'ligand_vina_out.pdbqt'))))
    with open('vina_prediction_rmsds.txt', 'w') as f:
        for prediction in pybel.readfile('pdbqt', 'ligand_vina_out.pdbqt'):
            f.write(f'{calculate_rmsd(gt, prediction)}\n')

if __name__ == '__main__':
    # When called from command line, take folder
    parser = argparse.ArgumentParser()
    parser.add_argument('folder') # e.g. /home/longyuxi/Documents/algebraic-topology-biomolecules/vina-benchmark/scpdb-example/1a2b_1
    input_args = parser.parse_args()
    run_on_folder(input_args.folder)