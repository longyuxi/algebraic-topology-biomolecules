"""Converts each pdb file to pqr file, as used by Cang et al. 2018
"""
import glob
import os
root = '/work/yl708/pdbbind/refined-set/'
pdb2pqr_executable = '/work/yl708/algebraic-topology-biomolecules/utils/pdb2pqr-linux-bin64-2.1.1/pdb2pqr'

pdb_files = glob.glob(root + '*/*.pdb')

for pdb_file in pdb_files:
    pqr_file = pdb_file[:pdb_file.rfind('.')] + '.pqr'
    # remove the REMARK line from the pdb file as it seems to confuse the utility

    command = 'grep -v "REMARK" ' + pdb_file + ' > tmpfile && rm ' + pdb_file + ' && mv tmpfile ' + pdb_file
    # e.g. grep -v "REMARK" /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb > tmpfile && rm /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb && mv tmpfile /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb
    # print(command)
    os.system(command)

    # example command:
    # pdb2pqr30 -ff=AMBER /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pqr
    command = pdb2pqr_executable + ' --ff=AMBER ' + pdb_file + ' ' + pqr_file

    # print(command)
    os.system(command)
