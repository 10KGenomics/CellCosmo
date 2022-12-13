
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
