# 安装说明

```bash
# 创建环境并安装所需要的第三方软件
ENV_NAME=CellCosmo
conda create -n $ENV_NAME -y python=3.9
conda activate $ENV_NAME
conda install star=2.6.1b -c bioconda -y        # STAR
conda install picard=2.18.17 -c bioconda -y     # picard
conda install subread=2.0.1 -c bioconda -y      # featureCounts
conda install ucsc-gtftogenepred -c bioconda -y # gtfToGenePred
# 获取并安装软件
pip install cell_cosmo-1.0.1.tar.gz

```
