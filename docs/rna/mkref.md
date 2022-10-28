## Features
---
通过参考基因组文件和gtf文件构建索引
## Output
---
`genome.config` config file
`{sample}.refFlat` gtfToGenePred output
`chrLength.txt` STAR output
`chrNameLength.txt` STAR output
`chrName.txt` STAR output
`chrStart.txt` STAR output
`exonGeTrInfo.tab` STAR output
`exonInfo.tab` STAR output
`geneInfo.tab` STAR output
`Genome` STAR output
`genomeParameters.txt` STAR output
`Log.out` STAR output
`SA` STAR output
`SAindex` STAR output
`sjdbInfo.txt` STAR output
`sjdbList.fromGTF.out.tab` STAR output
`sjdbList.out.tab` STAR output
`transcriptInfo.tab` STAR output


## Arguments
---
```bash
Usage: CellCosmo rna mkref [OPTIONS]

Options:
  --genome-name TEXT             genome name.  [required]
  --fasta TEXT                   Genome fasta file. Use absolute path or
                                 relative path to `genomeDir`.  [required]
  --gtf TEXT                     Genome gtf file. Use absolute path or
                                 relative path to `genomeDir`.  [required]
  --mt-gene-list TEXT            Mitochondria gene list file.  Use absolute
                                 path or relative path to `genomeDir`. It is a
                                 plain text file with one gene per line. If
                                 not provided, will use `MT-` and `mt-` to
                                 determine mitochondria genes.
  --genomeSAindexNbases INTEGER  STAR param: genomeSAindexNbases.  [default:
                                 14]
  --gene-name-as-name2           control the param of `-geneNameAsName2` in
                                 gtfToGenePred
  -t, --thread INTEGER           Threads to use.  [default: 6]
  -h, --help                     Show this message and exit.
```
