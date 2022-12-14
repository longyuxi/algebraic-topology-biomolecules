"""Converts each pdb file to pqr file, as used by Cang et al. 2018
"""
import glob
root = '/work/yl708/pdbbind/refined-set/'

pdb_files = glob.glob(root + '*/*.pdb')

for pdb_file in pdb_files:
    pqr_file = pdb_file[:pdb_file.rfind('.')] + '.pqr'
    # example command:
    # pdb2pqr30 -ff=AMBER /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pqr
    command = 'pdb2pqr30 -ff=AMBER ' + pdb_file + ' ' + pqr_file

    os.system(command)
