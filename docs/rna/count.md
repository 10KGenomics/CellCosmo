## Features
- Cell-calling: Distinguish valid cell barcodes from background barcodes. 
- Generate expression matrix.

## Output
- `{sample}_raw_feature_bc_matrix` The expression matrix of all detected barcodes in [Matrix Market Exchange Formats](https://math.nist.gov/MatrixMarket/formats.html).
- `{sample}_filtered_feature_bc_matrix` The expression matrix of valid cell barcodes in Matrix Market Exchange Formats. 
- `{sample}_count_detail.tsv` 4 columns: 
    - barcode  
    - gene ID
    - UMI count  
    - read_count  
- `{sample}_counts.tsv` 6 columns:
    - Barcode: barcode sequence
    - readcount: read count of each barcode
    - UMI2: read count with reads per UMI >= 2 for each barcode
    - UMI: UMI count for each barcode
    - geneID: gene count for each barcode
    - mark: valid cell barcode or backgound barcode.
        `CB` valid cell barcode
        `UB` background barcode  
- `{sample}_downsample.tsv` Subset a fraction of reads and calculate median gene number and sequencing saturation.

## Arguments
---
```bash
Usage: CellCosmo rna count [OPTIONS]

Options:
  --genomeDir TEXT               Genome directory.  [required]
  -b, --bam TEXT                 BAM file from feature_counts.  [required]
  --expected-cell-num INTEGER    Expected cell number.  [default: 3000]
  --cell-calling-method TEXT     Choose from [`auto`, `EmptyDrops_CR`]
                                 [default: EmptyDrops_CR]
  --n-umi-filter INTEGER         when correct barcode and umi, the umis count
                                 less than this value will be discard.set this
                                 params will accelerated running speed but
                                 some data will be discarded.  [default: 0]
  --barcode-correct-limit FLOAT  when barcode
                                 correct,low_umis_count/high_umis_count need
                                 less than this value.if set 1,merge low to
                                 high for all match case.  [default: 0.01]
  --umi-correct-limit FLOAT      when umi correct,low count/high count need
                                 less than this value.if set 1,merge low to
                                 high for all match case.  [default: 0.1]
  --force-cell-num TEXT          Force the cell number to be this number.
  -s, --sample TEXT              Sample name.  [required]
  -o, --outdir TEXT              Output directory.  [required]
  -t, --thread INTEGER           Thread to use.  [default: 4]
  -h, --help                     Show this message and exit.
```
