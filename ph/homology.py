# Dependency---------------------------------------------------------------------------------------
# Javaplex                                                                                        |
# download matlab-examples from https://github.com/appliedtopology/javaplex/releases/tag/4.3.1    |
# unzip, rename the folder as javeplex and put it in the same folder where this script is         |
#-------------------------------------------------------------------------------------------------|
# R-TDA                                                                                           |
# install TDA package in R and properly set the R library path                                    |
#-------------------------------------------------------------------------------------------------|
# Python packages                                                                                 |
# Numpy, Pickle                                                                                   |
#--------------------------------------------------------------------------------------------------

import os
import argparse
import TopBio.ReadFile.ReadMOL2 as ReadMOL2
import TopBio.ReadFile.ReadPDB as ReadPDB
import TopBio.ReadFile.ReadPQR as ReadPQR
import TopBio.PersistentHomology.PHSmallMolecule as PHSmallMolecule
import TopBio.PersistentHomology.PHComplex as PHComplex
import TopBio.Feature.LigandFeature as LigandFeature
import TopBio.Feature.ComplexFeature as ComplexFeature

def calculate(working_dir, ligand_name, protein_name, Clean = True):
    print('parameters:', working_dir, ligand_name, protein_name)
    # raise Exception

    a = ReadMOL2.SmallMolecule(ligand_name,working_dir)
    PHSmallMolecule.Level1_Rips(a, ligand_name, working_dir)

    PHSmallMolecule.Alpha(a, ligand_name, working_dir)
    ReadPDB.get_pdb_structure(a, 50.0, protein_name, working_dir)
    ReadPQR.get_pqr_structure(ligand_name, protein_name,working_dir)
    PHComplex.Interaction_Rips(50.0, protein_name, working_dir)
    PHComplex.Electrostatics_Rips(16.0, protein_name, working_dir)
    # PHComplex.Alpha(12.0, protein_name, working_dir)
    PHComplex.Alpha(50.0, protein_name, working_dir)

    # TopBP-ML
    LigandFeature.GenerateFeature_alpha(ligand_name, working_dir)
    LigandFeature.GenerateFeature_level1(ligand_name, working_dir)
    ComplexFeature.GenerateFeature_interaction_ML(protein_name, working_dir, 'all')
    ComplexFeature.GenerateFeature_electrostatics_ML(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'carbon')
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'heavy')

    # TopBP-DL
    ComplexFeature.GenerateFeature_interaction_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_electrostatics_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_2DCNN(protein_name, working_dir)

    # TopVS-ML
    LigandFeature.GenerateFeature_alpha(ligand_name, working_dir)
    LigandFeature.GenerateFeature_level1(ligand_name, working_dir)
    ComplexFeature.GenerateFeature_interaction_ML(protein_name, working_dir, 'pair')
    ComplexFeature.GenerateFeature_electrostatics_ML(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'carbon')
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'heavy')

    # TopVS-DL
    LigandFeature.GenerateFeature_alpha(ligand_name, working_dir)
    LigandFeature.GenerateFeature_level1(ligand_name, working_dir)
    ComplexFeature.GenerateFeature_interaction_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_electrostatics_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_1DCNN(protein_name, working_dir)

    # All features
    LigandFeature.GenerateFeature_alpha(ligand_name, working_dir)
    LigandFeature.GenerateFeature_level1(ligand_name, working_dir)
    ComplexFeature.GenerateFeature_interaction_ML(protein_name, working_dir, typ='all')
    ComplexFeature.GenerateFeature_interaction_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_electrostatics_ML(protein_name, working_dir)
    ComplexFeature.GenerateFeature_electrostatics_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'carbon')
    ComplexFeature.GenerateFeature_alpha_ML(protein_name, working_dir, 'heavy')
    ComplexFeature.GenerateFeature_alpha_1DCNN(protein_name, working_dir)
    ComplexFeature.GenerateFeature_alpha_2DCNN(protein_name, working_dir)

    if Clean:
        os.system('rm '+working_dir+'/*.PH')
        os.system('rm '+working_dir+'/*.csv')
        os.system('rm '+working_dir+'/*.pts')
        os.system('rm '+working_dir+'/*.bds')
        os.system('rm '+working_dir+'/tmp.out')

def test_example():
    working_dir = '/work/yl708/pdbbind/refined-set/4k7i/'
    ligand_name = '4k7i_ligand'
    protein_name = '4k7i_protein'
    calculate(working_dir, ligand_name, protein_name, Clean = True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir', type=str, help='directory containing the pdb file', default='/work/yl708/pdbbind/refined-set/4k7i/')
    parser.add_argument('--ligand_name', type=str, help='ligand file name without the suffix', default='4k7i_ligand')
    parser.add_argument('--protein_name', type=str, help='protein file name without the suffix', default='4k7i_protein')

    input_args = parser.parse_args()

    calculate(input_args.working_dir, input_args.ligand_name, input_args.protein_name)
