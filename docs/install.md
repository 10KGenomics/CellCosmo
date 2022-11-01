# 安装说明

## 首先创建环境，并且安装依赖的第三方软件

```bash

wget https://github.com/caijingtao1993/CellCOSMO_test/blob/a9607c295a2c13d77d216d88b3b7c3044047836d/dep_tools.txt
# dep_tools.txt file 
# bioconda::star=2.6.1b				      # STAR
# bioconda::picard=2.18.17          # picard
# bioconda::subread=2.0.1           # featureCounts
# bioconda::ucsc-gtftogenepred=377  # gtfToGenePred

ENV_NAME=CellCosmo
conda create -n $ENV_NAME -y --file dep_tools.txt

# 你也可以使用mamba安装环境
conda install mamba -y
mamba create -n $ENV_NAME -y --file dep_tools.txt
```

## 安装本软件
```bash
# 下载release 包
wget {todo after release the source code}
pip install cell_cosmo-1.0.1.tar.gz

# or you can build with source code
git clone xxx
python setup.py sdist
pip install dist/cell_cosmo-*.tar.gz
```
