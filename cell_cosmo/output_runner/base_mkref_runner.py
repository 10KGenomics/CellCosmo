#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_mkref_runner.py
@Time       : 2022/06/08
@Version    : 1.0
@Desc       : None
"""
import configparser
import sys
from cell_cosmo.output_runner.base_runner import BaseRunner
from cell_cosmo.util import GenomeUtil


class MKRefBase(BaseRunner):
    __genome_type__ = None
    __meta_files__ = None
    __meta_non_files__ = None

    def __init__(self, genome_name, fasta, thread=4, dry_run=False, **kwargs):
        super(MKRefBase, self).__init__()
        # genome_name, fasta,
        self._check_init_set()
        self.fasta = fasta  # required
        self.genome_name = genome_name  # required
        self.meta_files = ['fasta'] + list(self.__meta_files__)
        self.meta_non_files = ['genome_name'] + list(self.__meta_non_files__)
        self.dry_run = dry_run
        self.thread = thread
        self.genome_type = self.__genome_type__
        self.config_dict = {}
        self.config_file = GenomeUtil.GENOME_CONFIG  # 会输出到当前运行目录
        self._init(**kwargs)

    def _check_init_set(self):
        for meta in ('genome_type', 'meta_files', 'meta_non_files'):
            meta_name = f"__{meta}__"
            if self.get(meta_name) is None:
                print(f"must set {meta_name} in {self.__class__.__name__}")
                raise Exception(f"must set {meta_name} in {self.__class__.__name__}")
            # TODO meta 属性需要为可迭代对象,验证？

    def _init(self, **kwargs):
        self.config_dict['genome_type'] = self.genome_type
        for entry in self.meta_files + self.meta_non_files:
            if entry in {"fasta", "genome_name"}:
                value = self.get(entry)
                print(f'{entry} : {value}')
                self.config_dict[entry] = value
            elif entry not in kwargs:
                # 正常情况不会发生
                print(f"{entry} not set!!!{self.__class__.__name__}")
                raise Exception(f"{entry} not set!!!{self.__class__.__name__}")
            else:
                value = str(kwargs.pop(entry))
                print(f'{entry} : {value}')
                setattr(self, entry, value)
                self.config_dict[entry] = value

    def run(self):
        if self.dry_run:
            sys.exit(0)

    def _write_config(self):
        config = configparser.ConfigParser()
        config['genome'] = self.config_dict

        with open(self.config_file, 'w') as config_handle:
            config.write(config_handle)

    def __exit__(self, *args, **kwargs):
        self._write_config()
