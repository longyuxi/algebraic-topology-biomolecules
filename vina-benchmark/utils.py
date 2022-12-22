import pybel
import numpy as np

def calculate_bbox(mol: pybel.Molecule):
    """Generates bounding box for a Molecule

    Parameters
    ----------
    mol : pybel.Molecule
        The input molecule.

    Returns
    -------
    tuple of ndarray
        bbox_center, bbox_size
    """
    coords = []

    for atom in mol.atoms:
        coords.append(np.array(atom.coords))

    coords = np.array(coords)
    maxes = np.max(coords, axis=0)
    mins = np.min(coords, axis=0)

    bbox_center = (maxes + mins) / 2
    bbox_size = (maxes - mins)

    return bbox_center, bbox_size

def calculate_rmsd(mol1: pybel.Molecule, mol2: pybel.Molecule):
    """Calculates RMSD between two molecules without alignment

    Parameters
    ----------
    mol1 : pybel.Molecule
        molecule 1
    mol2 : pybel.Molecule
        molecule2

    Returns
    -------
    float
        RMSD
    """
    assert len(mol1.atoms) == len(mol2.atoms)

    atomwise_distances = []

    for idx in range(len(mol1.atoms)):
        atom1 = mol1.atoms[idx]
        atom2 = mol2.atoms[idx]
        assert atom1.type == atom2.type

        atom1_coords = np.array(atom1.coords)
        atom2_coords = np.array(atom2.coords)

        distance_squared = np.sum((atom1_coords - atom2_coords) ** 2)

        atomwise_distances.append(distance_squared)

    rmsd = np.sqrt(np.mean(np.array(atomwise_distances)))

    return rmsd
