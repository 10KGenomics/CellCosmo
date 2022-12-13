```bash
Usage: CellCosmo rna featureCounts [OPTIONS]

Options:
  -s, --sample TEXT            Sample name.  [required]
  --input TEXT                 BAM file path.  [required]
  -o, --outdir TEXT            Output directory.  [required]
  --gtf-type TEXT              Specify feature type (exon or gene) in GTF annotation
                               [default: exon]
  --genomeDir TEXT             Genome directory.  [required]
  -t, --thread INTEGER         Thread to use.  [default: 4]
  -d, --debug                  If this argument is used, may output additional
                               file for debugging.
  --feature-counts-param TEXT  Additional parameters for the  called software.
                               Need to be enclosed in quotation marks. For
                               example:  `--{software}_param "--param1 value1
                               --param2 value2"`.
  -h, --help                   Show this message and exit.
```
