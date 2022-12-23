#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : GenomeUtil.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""

GENOME_CONFIG = 'genome.config'  # 当前命令运行路径


def get_file_path(raw_file_path, genomeDir):
    if raw_file_path is None or str(raw_file_path) == 'None':
        return None

    file_path = raw_file_path
    if not file_path.startswith('/'):
        file_path = f'{genomeDir}/{file_path}'

    return file_path


def parse_dir(genome_dir, files=()):
    import configparser
    files = ('fasta',) + files

    config_file = f'{genome_dir}/{GENOME_CONFIG}'
    config = configparser.ConfigParser()
    with open(config_file, encoding='utf-8') as f:
        config.read_file(f)
        genome = dict(config['genome'])

        for entry in files:
            if entry not in genome:
                raise ValueError(f'{entry} not in {config_file}')
            genome[entry] = get_file_path(genome[entry], genome_dir)

    return genome


def parse_rna_dir(genome_dir):
    return parse_dir(genome_dir, files=('gtf', 'refflat', 'mt_gene_list'))
