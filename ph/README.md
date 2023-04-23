# Persistence Homology for Biomolecules

Currently a fork of code provided in *Representability of algebraic topology for biomolecules in machine learning based scoring and virtual screening* by *Zixuan Cang, Lin Mu, Guo-Wei Wei* ([PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005929)).

## This folder: Calculates the persistent homology of PDB files

Workflow:
1.  Make a folder elsewhere. Download and extract the *refined set* and the *index files* from [PDBBind](http://www.pdbbind.org.cn/download.php). In the folder for *refined set*, delete the non-PDB folders. Find them via `ls -l | grep -E -v '\s[0-9a-z]{4}$'`.
2. Download [pdb2pqr v2.1.0](https://github.com/Electrostatics/pdb2pqr/releases/tag/v2.1.0) and use `utils/make_pqrs.py` to convert things to `.pqr`.
3.



