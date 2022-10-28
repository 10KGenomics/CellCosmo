## Features
---
rna analysis pipeline
## Output
---

## Arguments

```bash
Usage: CellCosmo rna pipeline [OPTIONS]

Options:
  -g, --gen-config-tpl  Generate config template file.
  -c, --config TEXT     Config file for the pipeline.
  -h, --help            Show this message and exit.
```

## Detail

1. 生成配置命令

```bash
CellCosmo rna pipeline -g
# this command generate config template in current dir,
# which name is rna_pipeline.cfg
```

1. 关于配置文件说明，详见[]

1. 按制定配置信息运行全部命令

```bash
CellCosmo rna pipeline -c rna_pipeline.cfg
```