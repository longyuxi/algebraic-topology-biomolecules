"""Converts each pdb file to pqr file, as used by Cang et al. 2018
"""
import glob
import os
from tqdm import tqdm
root = '/work/yl708/pdbbind/refined-set/'
pdb2pqr_executable = '/work/yl708/algebraic-topology-biomolecules/ph/utils/pdb2pqr-linux-bin64-2.1.0/pdb2pqr'

pdb_files = glob.glob(root + '*/*.pdb')
cmds = ''

for pdb_file in tqdm(pdb_files):
    pqr_file = pdb_file[:pdb_file.rfind('.')] + '.pqr'
    tmpfile = pdb_file[:pdb_file.rfind('.')] + '.tmp'
    # remove the REMARK line from the pdb file as it seems to confuse the utility

    command1 = 'grep -v "REMARK" ' + pdb_file + ' > ' + tmpfile + ' && rm ' + pdb_file + ' && mv ' + tmpfile + ' ' + pdb_file
    # e.g. grep -v "REMARK" /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb > tmpfile && rm /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb && mv tmpfile /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb
    # os.system(command1)

    # example command:
    # pdb2pqr30 -ff=AMBER /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pdb /work/yl708/pdbbind/refined-set/4k7i/4k7i_protein.pqr
    command2 = pdb2pqr_executable + ' --ff=AMBER ' + pdb_file + ' ' + pqr_file

    # os.system(command2)
    cmds += command1 + ' ; ' + command2 + '\n'

with open('cmds.tmp', 'w') as f:
    f.write(cmds)
