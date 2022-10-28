# 10K scRNA-seq Data Analysis 
SingleCell RNA-seq Data analysis Software
## Introduction
CellCosmo is an open source and bioinfomatics analysis pipeline to high-throughput analysis Perseus Series single-cell RNA datasets.

## Installation
......(包括安装python版本、第三方开源软件、必备其他一些包)

## Start Running CellCosmo
#### 1. Build References  genomeDir for Homo sapiens and Mus musculus

#### 2.Generate configuration file(rna_pipeline.cfg)

生成配置文件

`CellCosmo rna pipeline -g`

[配置文件(`rna_pipeline.cfg`)说明](docs/rna/pipeline_config.md)
#### 3.Run pipeline
```bash
CellCosmo rna -c rna_pipeline.cfg
```

#### 4.Detailed docs can be found in manual
[sample](docs/rna/sample.md)
[barcode](docs/rna/barcode.md)
[cutadapt](docs/rna/cutadapt.md)
[star](docs/rna/star.md)
[featureCounts](docs/rna/featureCounts.md)
[count](docs/rna/count.md)
[analysis](docs/rna/analysis.md)


