# Persistence Homology for Biomolecules

Currently a fork of code provided in *Representability of algebraic topology for biomolecules in machine learning based scoring and virtual screening* by *Zixuan Cang, Lin Mu, Guo-Wei Wei* ([PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005929)).

## This folder: Calculates the persistent homology of PDB files

Workflow:
1.  **Downloading and preprocessing the PDBBind dataset**:
    1. Make a folder elsewhere. Download and extract the *refined set* and the *index files* from [PDBBind](http://www.pdbbind.org.cn/download.php).
    2. In the folder for *refined set*, delete the non-PDB folders. Find them via `ls -l | grep -E -v '\s[0-9a-z]{4}$'`.
    3. **Preprocessing the Download [pdb2pqr v2.1.0](https://github.com/Electrostatics/pdb2pqr/releases/tag/v2.1.0) and use `utils/make_pqrs.py` to convert things to `.pqr`.
2. **Starting a Redis database**: The Redis database is used for managing jobs submitted to the SLURM batch cluster.
    1. Make and install Redis via [https://redis.io/docs/getting-started/installation/install-redis-from-source/].
    2. Optionally, add the `src` folder of Redis to path.
    3. Create a `redis.conf` file somewhere and set a default password by putting e.g. `requirepass topology` in that file.
    4. Start the redis server on a host and adjust the `DB` constant in `dispatch_jobs.py` accordingly.
3. **Test run**: Test run with `test_example.py`. If successful, the folder being operated upon should end up with 20 files. No more, no less.
4. **Full run**: `dispatch_jobs.py`



