#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Args.py
@Time       : 2024/07/06
@Version    : 1.0
@Desc       : None
"""
from pathlib import Path
from dataclasses import dataclass, field, InitVar
from . import KeyName


@dataclass
class ValidParam:
    keyname: KeyName
    is_valid: bool
    is_output: bool
    n_allow: int
    prefix: InitVar[str]
    suffix: InitVar[str]
    # 包含%d的format，用于输出文件名称
    f1: str = field(init=False)
    f2: str = field(init=False)
    # temp_dir: str = field(init=False)
    # temp_format: str = field(init=False)
    id: int = field(init=False)
    name: str = field(init=False)

    def __post_init__(self, prefix, suffix):
        assert suffix == ".gz", "当前版本要求必须以gz格式输出"
        f = Path(prefix).absolute()
        p, s, = f.parent, f.name
        sniff = "" if self.keyname == KeyName.clean_num else self.keyname.name

        self.f1 = f"{p}/{s}_{sniff}1.fq{suffix}"
        self.f2 = f"{p}/{s}_{sniff}2.fq{suffix}"

        # 对于temp目录，需要提前创建
        # temp_dir = p / self.keyname.name
        # if not temp_dir.exists():
        #     temp_dir.mkdir()
        # self.temp_dir = str(temp_dir)
        # # 前补0，不应该写太多文件，需控制文件chunk的规模不能太小
        # self.temp_format = f"{temp_dir}/part-%06d_{s}_{sniff}%d.fq{suffix}"

        self.name = self.keyname.name
        self.id = self.keyname.value


class ValidPolytParam(ValidParam):
    def __init__(self, prefix, suffix, polyt_length, **kwargs):
        self.polyt_rate = kwargs.pop("polyt_rate", 0.7)
        self.polyt_length = polyt_length
        self.n_allow_for_polyt = polyt_length * self.polyt_rate
        self.use_polyt_valid_reads = kwargs.pop("use_polyt_valid_reads", False)
        self.output_no_polyt_reads = kwargs.pop("output_no_polyt_reads", True)
        super(ValidPolytParam, self).__init__(
            keyname=KeyName.no_polyt,
            is_valid=self.use_polyt_valid_reads,
            is_output=self.output_no_polyt_reads,
            n_allow=self.n_allow_for_polyt,
            prefix=prefix,
            suffix=suffix
        )


class ValidBarcodeParam(ValidParam):
    def __init__(self, prefix, suffix, **kwargs):
        self.n_allow_for_barcode = kwargs.pop("n_allow_for_barcode", 1)
        self.use_barcode_valid_reads = kwargs.pop("use_barcode_valid_reads", False)
        self.output_no_barcode_reads = kwargs.pop("output_no_barcode_reads", True)
        super(ValidBarcodeParam, self).__init__(
            keyname=KeyName.no_barcode,
            is_valid=self.use_barcode_valid_reads,
            is_output=self.output_no_barcode_reads,
            n_allow=self.n_allow_for_barcode,
            prefix=prefix,
            suffix=suffix
        )


class ValidLinkParam(ValidParam):
    def __init__(self, prefix, suffix, **kwargs):
        self.n_allow_for_link = kwargs.pop("n_allow_for_link", 2)
        self.use_link_valid_reads = kwargs.pop("use_link_valid_reads", False)
        self.output_no_linker_reads = kwargs.pop("output_no_linker_reads", True)
        super(ValidLinkParam, self).__init__(
            keyname=KeyName.no_linker,
            is_valid=self.use_link_valid_reads,
            is_output=self.output_no_linker_reads,
            n_allow=self.n_allow_for_link,
            prefix=prefix,
            suffix=suffix
        )


class ValidQualParam(ValidParam):
    def __init__(self, prefix, suffix, **kwargs):
        self.low_qual = kwargs.pop("low_qual", 0)
        self.low_num = kwargs.pop("low_num", 0)
        self.qual_offset = kwargs.pop("qual_offset", 33)
        self.n_allow_for_qual = self.low_num
        self.use_qual_valid_reads = self.low_qual != 0
        self.output_low_qual_reads = kwargs.pop("output_low_qual_reads", True)
        super(ValidQualParam, self).__init__(
            keyname=KeyName.low_qual,
            is_valid=self.use_qual_valid_reads,
            is_output=self.output_low_qual_reads,
            n_allow=self.n_allow_for_qual,
            prefix=prefix,
            suffix=suffix
        )


class CleanedReadsParam(ValidParam):
    def __init__(self, prefix, suffix):
        super(CleanedReadsParam, self).__init__(
            keyname=KeyName.clean_num,
            is_valid=True,
            is_output=True,
            n_allow=-1,
            prefix=prefix,
            suffix=suffix
        )
