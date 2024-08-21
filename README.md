# Introduction
CellCosmo is an open source and bioinfomatics analysis pipeline to high-throughput analysis Perseus Series single-cell RNA datasets.
* Propose
   * A collection of bioinfomatics analysis pipelines to process single cell sequencing data generated with Perseus Series single-cell RNA products.
* Language
   * Python3(>=3.9.*)
   * R scripts
* Hardware/Software environment
   * x86-64 compatible processors.
   * require at least 30GB of RAM and 4 CPU.
   * centos 7.x 64-bit operating system (Linux kernel 3.10.0, compatible with higher software and hardware configuration).
# Installation
Installation tutorial manual [here](docs/install.md)

# Start Running CellCosmo
## 1. Build References  genomeDir for Homo sapiens or Mus musculus
Build References genomeDir tutorial manual [here](docs/Build_References_genomeDir.md)

## 2.Generate configuration file(rna_pipeline.cfg)
Run the following command to get the cfg file:
```bash
CellCosmo rna pipeline -t starsolo -g 
```
* Detailed configuration documentation[`rna_starsolo_pipeline.cfg`](docs/rna/rna_pipeline.cfg)
* If already established References  genomeDir and Kit version1 default, You can apply the simplified version of configuration documentation[`rna_pipeline_sv.cfg`](docs/rna/rna_pipeline_sv.cfg)
## 3.Run pipeline
```bash
CellCosmo rna pipeline -t starsolo -c rna_pipeline.cfg
```

## 4.Detailed docs can be found in manual
* (1) Sample sequencing data preparation  
  [CellCosmo rna sample](docs/rna/sample.md)  
* (2) Rna alignment
  [CellCosmo rna starsolo](docs/rna/starsolo.md)
* (3) Cell cluster and web report result output  
  [CellCosmo rna analysis](docs/rna/analysis.md)  

# Support
The officially supported release binaries are available at: (http://www.10kgenomics.com/)

