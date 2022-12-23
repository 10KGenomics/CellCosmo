#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : GTFDictUtil.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import os
import re
import logging
from collections import Counter
from cell_cosmo.util import runtime, reader

logger = logging.getLogger(__name__)


class GTFDict(dict):
    """
    key: gene_id
    value: gene_name
    If the key does not exist, return key. This is to avoid the error:
        The gtf file contains one exon lines with a gene_id, but do not contain a gene line with the same gene_id. FeatureCounts
        work correctly under this condition, but the gene_id will not appear in the Gtf_dict.
    """

    def __init__(self, gtf_file):
        super().__init__()
        self.gtf_file = gtf_file
        self.load_gtf()

    @runtime(f"{__name__}.load_gtf")
    def load_gtf(self):
        """
        get gene_id:gene_name from gtf file
            - one gene_name with multiple gene_id: "_{count}" will be added to gene_name.
            - one gene_id with multiple gene_name: error.
            - duplicated (gene_name, gene_id): ignore duplicated records and print a warning.
            - no gene_name: gene_id will be used as gene_name.

        Returns:
            {gene_id: gene_name} dict
        """

        gene_id_pattern = re.compile(r'gene_id "(\S+)";')
        gene_name_pattern = re.compile(r'gene_name "(\S+)"')
        id_name = {}
        c = Counter()
        for line, *_ in reader(self.gtf_file, ignore_test_env=True):
            if not line.strip():
                continue
            if line.startswith('#'):
                continue
            tabs = line.split('\t')
            gtf_type, attributes = tabs[2], tabs[-1]
            if gtf_type == 'gene':
                gene_id = gene_id_pattern.findall(attributes)[-1]
                gene_names = gene_name_pattern.findall(attributes)
                if not gene_names:
                    gene_name = gene_id
                else:
                    gene_name = gene_names[-1]
                c[gene_name] += 1
                if c[gene_name] > 1:
                    if gene_id in id_name:
                        assert id_name[gene_id] == gene_name, (
                            'one gene_id with multiple gene_name '
                            f'gene_id: {gene_id}, '
                            f'gene_name this line: {gene_name}'
                            f'gene_name previous line: {id_name[gene_id]}'
                        )
                        logger.warning(
                            'duplicated (gene_id, gene_name)'
                            f'gene_id: {gene_id}, '
                            f'gene_name {gene_name}'
                        )
                        c[gene_name] -= 1
                    else:
                        gene_name = f'{gene_name}_{c[gene_name]}'
                id_name[gene_id] = gene_name
        self.update(id_name)

    def __getitem__(self, key):
        """if key not exist, return key"""
        return dict.get(self, key, key)
