# cuomics
CUDA powered GWAS tool on GPU cluster

# Script to build cuomics from source

### Build from Source

To install cuomics from source, ensure the dependencies are met and follow the steps below:

* Clone the repository and submodules

```sh
CUOMICS_HOME=$(pwd)/cuomics
git clone https://github.com/prasunanand/cuomics.git $CUOMICS_HOME
cd CUOMICS_HOME
git submodule update --init --remote --recursive
```

* Create the conda development environment `cuomics_dev`

```sh
# create the conda environment (assuming in base `cuomics` directory)
conda env create --name cuomics_dev --file conda/environments/cuomics_dev.yml
# activate the environment
source activate cuomics_dev
```

* Install PyVCF

```sh
pip install PyVCF
```

# LICENSE

This software is distributed under the [Apache License 2.0](LICENSE).

Copyright © 2019, Modak Analytics Ltd.
