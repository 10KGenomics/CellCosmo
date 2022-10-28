# 10K scRNA-seq Data Analysis 
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


刘经理，我粗略整理了一下需要上传github的材料和说明，主要：
1.关于READNME.MD文件的构建，主要包括软件的介绍、软件安装、软件的运行、软件的每个模块说明；
2.关于一些关键性结果文件的说明，如：
Sample_summary.tsv 、Sample_raw_feature_bc_matrix、Sample_filtered_feature_bc_matrix、
Sample_count_detail.tsv、Sample_counts.tsv
3.可以参考10X、华大、新格元的，或者其他，相信这方面您比我有经验的多，后面我们也打算设置一个链接，
跳到我们公司官网，官网上我们后期自己也会做个软件介绍和教程；
4.关于一些非必要中间文件的删除，：
20220722CYM2
20220722CYM2.summary
20220722CYM2_1.fq.gz
20220722CYM2_2.fq.gz
20220722CYM2_Aligned.out.bam
20220722CYM2_Aligned.sortedByCoord.out.bam
20220722CYM2_Aligned.sortedByCoord.out.bam.bai
20220722CYM2_Aligned.sortedByCoord.out.bam.featureCounts.bam
20220722CYM2_clean_2.fq.gz
20220722CYM2_count_detail.raw.tsv
20220722CYM2_downsample.tsv
20220722CYM2_no_polyt1.fq.gz
20220722CYM2_no_polyt2.fq.gz
cutadapt.log
stat.txt