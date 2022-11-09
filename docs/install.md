# Installation Instructions

## Step 1: Create a conda environment, Install dependent third-party software

```bash
wget https://github.com/caijingtao1993/CellCOSMO_test/edit/main/dep_tools.txt
# dep_tools.txt file 
# bioconda::star=2.7.10a            # STAR
# bioconda::picard=2.18.17          # picard
# bioconda::subread=2.0.1           # featureCounts
# bioconda::ucsc-gtftogenepred=377  # gtfToGenePred
# bioconda::samtools=1.12           # samtools

ENV_NAME=CellCosmo
conda create -n $ENV_NAME -y --file dep_tools.txt

# You can also use mamba installation environment
conda install mamba -y
mamba create -n $ENV_NAME -y --file dep_tools.txt
```

## Step 2: Installing CellCosmo software
```bash
# Download the release package
wget {todo after release the source code}
pip install cell_cosmo-1.0.1.tar.gz

# or you can build with source code
git clone https://github.com/caijingtao1993/CellCOSMO_test.git
python setup.py sdist
pip install dist/cell_cosmo-*.tar.gz
```
