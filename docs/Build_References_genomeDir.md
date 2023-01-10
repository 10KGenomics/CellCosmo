# Make a rna genomeDir

- Construction of human reference genome
Remarks：If there is gene_name in the gtf file, you can add the "--gene-name-as-name2" parameter in the "CellCosmo rna mkref" directive. If there is only "gene_id", you will not add the parameter.

```bash
mkdir GRCh38_ensembl_index
cd GRCh38_ensembl_index

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

CellCosmo rna mkref --gtf Homo_sapiens.GRCh38.99.gtf \
--genome-name GRCh38 \
--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gene-name-as-name2
```
- Construction of mouse reference genome

```
mkdir GRCm38_ensembl_index
cd GRCm38_ensembl_index

wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz

CellCosmo rna mkref --gtf Mus_musculus.GRCm38.99.gtf \
--genome-name GRCm38 \
--fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
--gene-name-as-name2
```

