## Features
- Alignment the sequencing R2 RNA reads to the reference genome with STAR software.

## Arguments
```bash
Usage: CellCosmo rna starsolo [OPTIONS]

Options:
  --fq1 TEXT                       R1 fastq file. Multiple files are separated by comma.  [required]
  --fq2 TEXT                       R2 fastq file. Multiple files are separated by comma.  [required]
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
  --consensus-fq                  A indicator that the input fastq has been
                                  consensused.
  --genomeDir TEXT                Genome directory.  [required]
  --star-mem INTEGER              Maximum memory that STAR can use.  [default:
                                  30]
  --star-param TEXT               Additional parameters for the called
                                  software. Need to be enclosed in quotation
                                  marks. For example: `--star_param "--param1
                                  value1 --param2 value2"`.
  --out-filter-match-n-min INTEGER
                                  Alignment will be output only if the number
                                  of matched basesis higher than or equal to
                                  this value.  [default: 50]
  --sam_attributes TEXT           Additional attributes(other than NH HI nM AS
                                  CR UR CB UB GX GN ) to be added to SAM file
                                  [default: NH HI nM AS CR UR CB UB GX GN ]
  --solo-features TEXT            The same as the soloFeatures argument in
                                  STARsolo  [default: Gene GeneFull_Ex50pAS]
  --out_sam_type TEXT             type of SAM/BAM output  [default: BAM
                                  SortedByCoordinate]
  --solo_cell_filter_method TEXT  Cellcalling Method  [default: EmptyDrops_CR]
  --solo_cell_filter_n_expect INTEGER
                                  Expect_num  [default: 3000]
  --solo_cell_filter_args TEXT    solo cell filter args  [default: 0.99 10
                                  45000 90000 500 0.01 20000 0.001 10000]
  --adapter_3p TEXT               Adapter sequence to clip from 3 prime.
                                  Multiple sequences are seperated by space
                                  [default: AAAAAAAAAAAA]
  --out-unmapped                  Output unmapped reads.
  -o, --outdir TEXT               Output directory.  [required]
  -s, --sample TEXT               Sample name.  [required]
  -t, --thread INTEGER            Thread to use.  [default: 4]
  -h, --help                      Show this message and exit.
```
