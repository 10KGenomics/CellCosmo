# Quick start
本软件用来分析RNA 单细胞分析

# Usage Example
## 首先构建索引

- 构建人的参考基因组

```bash
mkdir hg38_ensembl_index
cd hg38_ensembl_index
wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz
CellCosmo rna mkref --gtf Homo_sapiens.GRCh38.99.gtf \
--genome-name humo \
--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gene-name-as-name2
```
- 构建鼠的参考基因组

```
mkdir mm10_ensembl_index
cd mm10_ensembl_index
wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz
CellCosmo rna mkref --gtf Mus_musculus.GRCm38.99.gtf \
--genome-name humo \
--fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
--gene-name-as-name2
```

## 运行流程

