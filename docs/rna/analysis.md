## Features
- Cell clustering with Scanpy.
- Calculate the marker genes of each cluster.

## Output
- `{sample}_tsne_coord.tsv` t-SNE coordinates and clustering information.
- `{sample}_umap_coord.tsv` UMAP  coordinates and clustering information.
- `{sample}_analysis` The analysis partial results are saved in this folder, include PCA, t-SNE, UMAP, diffGene expression and clustering. 

## Arguments
---
```bash
Usage: CellCosmo rna analysis [OPTIONS]

Options:
  --genomeDir TEXT      Genome directory.  [required]
  --matrix-file TEXT    Matrix_10X directory from step count.  [required]
  -s, --sample TEXT     Sample name.  [required]
  -o, --outdir TEXT     Output directory.  [required]
  -t, --thread INTEGER  Thread to use.  [default: 4]
  -h, --help            Show this message and exit.
```
