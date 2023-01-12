## Features
rna analysis pipeline

## Output
- {SampleName}_1.fq.gz
- {SampleName}_2.fq.gz
- {SampleName}_clean_2.fq.gz
- cutadapt.log
- {SampleName}_Aligned.out.bam
- {SampleName}_SJ.out.tab
- {SampleName}_Log.progress.out
- {SampleName}_Log.out
- {SampleName}_Log.final.out
- {SampleName}_Aligned.sortedByCoord.out.bam
- {SampleName}_Aligned.sortedByCoord.out.bam.bai
- {SampleName}_region.log
- {SampleName}.summary
- {SampleName}
- {SampleName}_Aligned.sortedByCoord.out.bam.featureCounts.bam
- {SampleName}_name_sorted.bam
- {SampleName}_count_detail.raw.tsv
- {SampleName}_barcode_umi_less_than_5_filter_27269_data.tsv
- {SampleName}_barcode_itd_nbr_pairs.tsv
- {SampleName}_barcode_itd_nbr_pairs.raw.tsv
- {SampleName}_count_detail.tsv
- {SampleName}_counts.tsv
- {SampleName}_raw_feature_bc_matrix/
- {SampleName}_filtered_feature_bc_matrix/
- {SampleName}_downsample.tsv
- {SampleName}_analysis
- {SampleName}_tsne_coord.tsv
- {SampleName}_umap_coord.tsv
- stat.txt
- {SampleName}_summary.tsv

Note: See the description of each module for the above output file notes.

## Arguments
```bash
Usage: CellCosmo rna pipeline [OPTIONS]

Options:
  -g, --gen-config-tpl  Generate config template file.
  -c, --config TEXT     Config file for the pipeline.
  -h, --help            Show this message and exit.
```

## Detail

#### 1. Generate configuration file(rna_pipeline.cfg)

```bash
CellCosmo rna pipeline -g
# this command generate config template in current dir,
# which name is rna_pipeline.cfg
```

- configuration documentation, see [rna_pipeline.cfg](rna_pipeline.cfg)

#### 2. Run pipeline

```bash
CellCosmo rna pipeline -c rna_pipeline.cfg
```
