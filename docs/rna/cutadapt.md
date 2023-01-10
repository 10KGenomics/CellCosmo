## Features
- Trim adapters in R2 reads with cutadapt.

## Output
- cutadapt.log Cutadapt output log file.
- `{sample}_clean_2.fq.gz` R2 reads file after trimming adapters.

## Arguments
```bash
Usage: CellCosmo rna cutadapt [OPTIONS]

Options:
  --fq TEXT                 Required. R2 reads from step Barcode.  [required]
  --gzip                    Output gzipped fastq files.
  --adapter-fasta TEXT      Additional adapter fasta file.
  --minimum-length INTEGER  Discard processed reads that are shorter than
                            LENGTH.  [default: 20]
  --nextseq-trim INTEGER    Quality trimming of reads using two-color
                            chemistry (NextSeq). Some Illumina instruments use
                            a two-color chemistry to encode the four bases.
                            This includes the NextSeq and the NovaSeq. In
                            those instruments, a `dark cycle` (with no
                            detected color) encodes a G. However, dark cycles
                            also occur when sequencing `falls off` the end of
                            the fragment. The read then contains a run of
                            high-quality, but incorrect `G` calls at its 3'
                            end.  [default: 20]
  --overlap INTEGER         Since Cutadapt allows partial matches between the
                            read and the adapter sequence, short matches can
                            occur by chance, leading to erroneously trimmed
                            bases. For example, roughly 0.25 of all reads end
                            with a base that is identical to the first base of
                            the adapter. To reduce the number of falsely
                            trimmed bases, the alignment algorithm requires
                            that at least {overlap} bases match between
                            adapter and read.  [default: 10]
  --insert INTEGER          Read2 insert length.  [default: 150]
  --cutadapt-param TEXT     Other cutadapt parameters. For example,
                            --cutadapt_param '-g AAA'
  -o, --outdir TEXT         Output directory.  [required]
  -s, --sample TEXT         Sample name.  [required]
  -t, --thread INTEGER      Thread to use.  [default: 4]
  -h, --help                Show this message and exit.
```
