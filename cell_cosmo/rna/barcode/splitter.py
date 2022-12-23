#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : splitter.py
@Time       : 2022/05/27
@Version    : 1.0
@Desc       : None
"""
import sys

from cell_cosmo.rna.barcode.valid import BaseOutput, ValidCfg, ValidPolyt, ValidLink, ValidBarcode, ValidQual
from cell_cosmo.rna.barcode import parse_pattern as pp
from collections import Counter
import logging

logger = logging.getLogger(__name__)


class BarcodeSplitter(BaseOutput):
    __filename__ = ""

    def __init__(self, chemistry_name, chemistry_config, pattern, **kwargs):
        """

        :param chemistry_name: Library name
        :param chemistry_config: Library config file path
        :param pattern: Library location mode string,eg:C5U3L15U3C6U3L6C5T30
        """
        super(BarcodeSplitter, self).__init__(need_valid=True, need_output=True, **kwargs)
        # 优先从文库名中获取相关信息
        if chemistry_name is not None:
            if chemistry_config is not None:
                logger.warning("`--chemistry-config` param invalid when specify `--chemistry-name`")
            # TODO 待完成文库标准化再开发使用 library_name 参数
            self.library_name = chemistry_name
            self.config_file = ""  # 从数据库中获取该信息,
            logger.info("chemistry database is developing,please wait!!")
            sys.exit(1)
        else:
            self.library_name = pattern
            self.config_file = chemistry_config
        self.pattern_str = pattern

        # init
        self.qual_counter_barcode = Counter()
        self.qual_counter_read = Counter()
        self.qual_counter_umi = Counter()
        self.num_total = 0

        self.pattern = pp.PatternParser(self.pattern_str)
        self.c_len = self.pattern.get_pattern_len(pp.C)
        self.checker4polyt = ValidPolyt(pattern=self.pattern, **kwargs)
        self.checker4link = ValidLink(pattern=self.pattern, **kwargs)
        self.checker4barcode = ValidBarcode(pattern=self.pattern, **kwargs)
        self.checker4qual = ValidQual(pattern=self.pattern, **kwargs)
        self.valid_chain = self._init_valid_chain()

    def _init_valid_chain(self):
        vs = []
        if self.checker4polyt.need_valid:
            vs.append(self.checker4polyt)
            logger.info(f"user polyT valid reads,"
                        f"T base >= {self.checker4polyt.n_allow}")
        if self.checker4barcode.need_valid or self.checker4link.need_valid:
            if self.config_file is None:
                raise Exception(
                    "use barcode or link valid reads,"
                    "`--chemistry-name` or `--chemistry-config` must be specified")
            self.v_cfg = ValidCfg(config_file=self.config_file)
            # 初始化barcode和valid的验证集
            if self.checker4barcode.need_valid:
                vs.append(self.checker4barcode)
                logger.info(f"user barcode valid reads,allow mismatch "
                            f"base num is {self.checker4barcode.n_allow}")
                self.checker4barcode.files = self.v_cfg.barcode_files
            if self.checker4link.need_valid:
                vs.append(self.checker4link)
                logger.info(f"user link valid reads,allow mismatch base "
                            f"num is {self.checker4link.n_allow}")
                self.checker4link.files = self.v_cfg.link_files
        if self.checker4qual.need_valid:
            vs.append(self.checker4qual)
            logger.info(f"user qual valid reads,"
                        f"low_num is {self.checker4qual.low_num},"
                        f"low_qual is {self.checker4qual.low_qual}")
        return vs

    def get_new_header(self, r1):
        # out_fq2.write(f'@{cb}_{umi}_{self.total_num}\n{seq2}\n+\n{qual2}\n')
        b = self.pattern.get_seq_str(r1[1], pp.C)
        umi = self.pattern.get_seq_str(r1[1], pp.U)
        qual_b = self.pattern.get_seq_str(r1[3], pp.C)
        qual_umi = self.pattern.get_seq_str(r1[3], pp.U)
        self.qual_counter_barcode.update(qual_b)
        self.qual_counter_umi.update(qual_umi)
        if self.checker4barcode.need_valid and \
                b not in self.checker4barcode.barcode_set:
            # 对barcode进行校正
            b = self.checker4barcode.barcode_mismatch_dit[b]
        return f"@{b}_{umi}_{self.num_total}"

    def set_read_qual(self, r2):
        self.qual_counter_read.update(r2[3])

    def process_reads_pair(self, r1, r2):
        self.num_total += 1
        if all([not v.need_valid for v in self.valid_chain]):
            # 不验证任何条件,则所有序列都是有效的
            new_header = self.get_new_header(r1)
            self.set_read_qual(r2)
            new_r1 = (new_header, r1[1], r1[2], r1[3])
            new_r2 = (new_header, r2[1], r2[2], r2[3])
            self._write(new_r1, new_r2, True)
            return True
        for v in self.valid_chain:
            if not v.valid(r1, r2):
                return False
        new_header = self.get_new_header(r1)
        self.set_read_qual(r2)
        new_r1 = (new_header, r1[1], r1[2], r1[3])
        new_r2 = (new_header, r2[1], r2[2], r2[3])
        self._write(new_r1, new_r2, True)
        return True

    def close_all(self):
        for v in self.valid_chain:
            v._close()
        self._close()
