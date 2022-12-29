#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : MKRef4STAR.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import logging
import subprocess
from cell_cosmo.output_runner.base_mkref_runner import MKRefBase
from cell_cosmo.util import runtime

logger = logging.getLogger(__name__)


class MKRef4RNA(MKRefBase):
    __genome_type__ = "rna"
    __meta_files__ = ['gtf', 'mt_gene_list', 'refflat']
    __meta_non_files__ = ['genomeSAindexNbases']

    def __init__(self, **kwargs):
        kwargs.update({'refflat': f'{kwargs.get("genome_name")}.refFlat'})
        self.use_gene_name_as_name2 = kwargs.get("gene_name_as_name2", False)
        super(MKRef4RNA, self).__init__(**kwargs)

    @runtime(__name__)
    def run(self):
        super().run()
        self.build_star_index()
        self.build_refflat()

    @runtime(f"{__name__}.build_star_index")
    def build_star_index(self):
        cmd = (
            f'STAR \\\n'
            f'--runMode genomeGenerate \\\n'
            f'--runThreadN {self.thread} \\\n'
            f'--genomeDir ./ \\\n'
            f'--genomeFastaFiles {self.fasta} \\\n'
            f'--sjdbGTFfile {self.get("gtf")} \\\n'
            f'--sjdbOverhang 100 \\\n'
        )
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @runtime(__name__)
    def build_refflat(self):
        gene_name_alias = "-geneNameAsName2" if self.use_gene_name_as_name2 else ""
        cmd = (
            f'gtfToGenePred -genePredExt {gene_name_alias} \\\n'
            f'{self.get("gtf")} /dev/stdout | \\\n'
            'awk \'{print $12"\\t"$1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10}\' \\\n'
            f'> {self.get("refflat")} \\\n'
        )
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)