## Features
- Alignment the sequencing R2 RNA reads to the reference genome with STAR software.

## Output
- `{sample}_Aligned.sortedByCoord.out.bam` The uniquely aligned reads, sorted by coordinate, in BAM format.
- `{sample}_Log.final.out` A summary of mapping statistics for the sample.
- `{sample}_Log.out` A running log from STAR, withinformation about the run.
- `{sample}_Log.progress.out` Job progress with the number of processed reads, % of mapped reads etc., updated every ~1 minute.
- `{sample}_SJ.out.tab` High confidence collapsed splice junctions in tab-delimited format. Only junctions supported by uniquely mapping reads are reported.

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
