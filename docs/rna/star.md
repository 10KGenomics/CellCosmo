## Features
---

## Output
---

## Arguments
```bash
Usage: CellCosmo rna star [OPTIONS]

Options:
  --fq TEXT                       R2 fastq file.  [required]
  --consensus-fq                  A indicator that the input fastq has been
                                  consensused.
  --genomeDir TEXT                Genome directory.  [required]
  --out-unmapped                  Output unmapped reads.
  --out-filter-match-n-min INTEGER
                                  Alignment will be output only if the number
                                  of matched basesis higher than or equal to
                                  this value.  [default: 0]
  --out-filter-multimap-n-max INTEGER
                                  How many places are allowed to match a read
                                  at most.  [default: 1]
  --star-mem INTEGER              Maximum memory that STAR can use.  [default:
                                  30]
  --star-param TEXT               Additional parameters for the called
                                  software. Need to be enclosed in quotation
                                  marks. For example: `--{software}_param "--
                                  param1 value1 --param2 value2"`.
  -o, --outdir TEXT               Output directory.  [required]
  -s, --sample TEXT               Sample name.  [required]
  -t, --thread INTEGER            Thread to use.  [default: 4]
  -h, --help                      Show this message and exit.
```
