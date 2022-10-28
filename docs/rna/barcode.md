## Features
---

## Output
---

## Arguments
---
```bash
Usage: CellCosmo rna barcode [OPTIONS]

  split barcode

Options:
  --fq1 TEXT                      Fastq file containing library information,
                                  the file suffix must be .fastq.gz or .fq.gz
                                  [required]
  --fq2 TEXT                      Another fastq file that paired with fq1
                                  [required]
  -s, --sample TEXT               Sample name.  [required]
  -n, --chemistry-name TEXT       Chemistry name,If this parameter is
                                  specified, the program will load the
                                  chemistry information corresponding to the
                                  name from the preset database
  -c, --chemistry-config TEXT     Chemistry configuration file path, the same
                                  level directory contains chemistry
                                  information such as barcode library and link
                                  string
  -p, --pattern TEXT              The pattern of sequences,eg:`C5U3L15U3C6U3L6C5T30`,
                                  The number after the letter represents the number of bases
                                  - `C`: barcode
                                  - `L`: linker
                                  - `U`: UMI
                                  - `T`: poly T
  --use-link-valid-reads          Validate reads using link sequences
  --use-barcode-valid-reads       Validate reads using barcode sequences
  --use-polyt-valid-reads         Validate reads using polyt sequences
  --allow-link-diff-num INTEGER   mismatch number with link sequence,this
                                  parameter only takes effect when `--use-
                                  link-valid-reads` is specified  [default: 2]
  --allow-barcode-diff-num INTEGER
                                  mismatch number with barcode library,this
                                  parameter is invalid when `--use-barcode-
                                  valid-reads` is specified  [default: 1]
  --low-qual INTEGER              The barcode and UMI whose phred value are
                                  lower than --low-qual will be regarded as
                                  low-quality bases.  [default: 0]
  --low-num INTEGER               The maximum number of low qual bases allowed
                                  in barcode and UMI.  [default: 2]
  --polyt-rate FLOAT              The proportion of T bases in polyT that need
                                  to be satisfied,range is [0,1].When
                                  specified as 0, the limit of polyT is
                                  ignored;When specified as 1, it means all T
                                  bases in polyT;When the specified pattern
                                  contains T configuration, the default value
                                  of this parameter is 2/3, otherwise the
                                  default value is 0  [default: 0.7]
  --gzip                          Output fastq file in compressed format
  --output-r1                     Output the R1 sequence corresponding to the
                                  valid sequence
  -o, --outdir TEXT               Output directory.  [required]
  -t, --thread INTEGER            Thread to use.  [default: 4]
  -h, --help                      Show this message and exit.
```