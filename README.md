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
[Build References genomeDir](docs/Build_References_genomeDir.md)

## 2.Generate configuration file(rna_pipeline.cfg)

生成配置文件

`CellCosmo rna pipeline -g`

[配置文件(`rna_pipeline.cfg`)说明](docs/rna/pipeline_config.md)
## 3.Run pipeline
```bash
CellCosmo rna -c rna_pipeline.cfg
```

## 4.Detailed docs can be found in manual
[sample](docs/rna/sample.md)  
[barcode](docs/rna/barcode.md)  
[cutadapt](docs/rna/cutadapt.md)  
[star](docs/rna/star.md)  
[featureCounts](docs/rna/featureCounts.md)  
[count](docs/rna/count.md)  
[analysis](docs/rna/analysis.md)  

# Support
The officially supported release binaries are available at: (http://www.10kgenomics.com/)

