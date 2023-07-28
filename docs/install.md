# Installation Instructions

## Step 1: Create a conda environment, Install dependent third-party software

```bash
wget https://github.com/caijingtao1993/CellCOSMO_test/edit/main/dep_tools.txt
# dep_tools.txt file
conda-forge::python=3.10
bioconda::star=2.7.10b              # STAR
bioconda::subread=2.0.1             # featureCounts
bioconda::picard=2.18.17            # picard
bioconda::ucsc-gtftogenepred=447    # gtfToGenePred
bioconda::samtools=1.12             # samtools       

ENV_NAME=CellCosmo
conda create -n $ENV_NAME -y --file dep_tools.txt

# You can also use mamba installation environment
conda install mamba -y
mamba create -n $ENV_NAME -y --file dep_tools.txt
```

## Step 2: Installing CellCosmo software
```bash
# Download the release package
wget https://github.com/10KGenomics/CellCOSMO/releases/download/v1.0.7/cell_cosmo-1.0.7.tar.gz
pip install cell_cosmo-1.0.7.tar.gz

# or you can build with source code
git clone https://github.com/caijingtao1993/CellCOSMO.git
python setup.py sdist
pip install dist/cell_cosmo-*.tar.gz
```
