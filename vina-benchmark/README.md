## Description

This repository benchmarks the performance of AutoDock Vina on the scPDB dataset.

## Software dependencies:

1. [vina](https://vina.scripps.edu/downloads/)
2. meeko (pip)
    - Also need to install scipy numpy pandas
3. Openbabel with Python bindings
    ```
    wget https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-3-1-1.tar.gz
    tar -xzvf openbabel-3-1-1.tar.gz
    cd openbabel-openbabel-3-1-1
    cd ..
    mv openbabel-openbabel-3-1-1 ob
    rm openbabel-3-1-1.tar.gz
    mkdir build
    cd build
    cmake ../ob -DPYTHON_BINDINGS=ON
    make -j8
    sudo make install
    ```
    and then add `export PYTHONPATH=/usr/local/lib:$PYTHONPATH` to `~/.zshrc` or `~/.bashrc`.

## Description of scripts in this folder

1. The `run_on_folder()` function in `vinascpdb.py` is the primary file to run. This function can also be called in bash shell by calling `python vinascpdb.py --folder/to/scpdb/such_as_1a2b_1`. This function generates files into the folder being run on.

2. `dispatch_jobs` accesses a dataframe stored in `jobs.csv` (and creates it if it doesn't exist), which stores the running status of jobs. `dispatch_jobs` then calls `job_wrapper` to run jobs.

3. `benchmark_results` and `calculate_center_distances` plot results.

## To run this repository

1. Download the scPDB dataset.
2. Install all the software dependencies.
3. Change the paths appropriately. I have specified all the paths variables, so they can be easily changed. Run `rg longyuxi` or `rg yl708` to identify these variables.
4. You are ready to run the aforementioned scripts!