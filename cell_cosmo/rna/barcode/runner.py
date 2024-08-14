#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : barcode_ck.py
@Time       : 2022/05/26
@Version    : 1.0
@Desc       : None
"""
import os
import logging
from multiprocessing import Pool, Manager
from cell_cosmo.output_runner import BaseReportRunner
from cell_cosmo.util import runtime, fmt_number, get_threads
from cell_cosmo.tools.chemistry import LibraryInfo
from .stat_info import StatInfo
from .validators import Validators
from .reads_processor import reads_processor
from .yield_batch_data import yield_batch_data
from functools import partial

logger = logging.getLogger(__name__)


class BarcodeSplitterRunner(BaseReportRunner):
    _STEP_NAME = "barcode"
    _DISPLAY_TITLE = "Demultiplexing"

    def __init__(self, fq1, fq2, pattern=None, chemistry_config=None, chemistry_name=None, **kwargs):
        super(BarcodeSplitterRunner, self).__init__(**kwargs)

        self.fq1 = fq1
        self.fq2 = fq2
        self.batch_size = kwargs.pop("batch_size", 1000 * 1000 * 100)
        # self.cache_limits = kwargs.pop("cache_limits", 500)
        # 这里要求程序运行线程必须大于5(有5个写进程),这里限制为8
        self.thread = get_threads(kwargs.pop("thread", 8), min_limit=8)
        self.library_info = LibraryInfo(
            chemistry_name, chemistry_config, pattern,
            barcode_n_mismatch=kwargs.get("n_allow_for_barcode", 1),
            link_n_mismatch=kwargs.get("n_allow_for_link", 2)
        )

        self.validators = Validators(
            n_polyt=self.library_info.n_T,
            batch_size=self.batch_size, **kwargs
        )
        # 用于统计的属性
        self.stat_info = StatInfo()

    def update_stat_info(self, stat_info: StatInfo):
        self.stat_info.update(stat_info)

    @runtime(__name__)
    def run(self):
        with Pool(self.thread) as pool:
            res = pool.imap(
                partial(reads_processor, self.validators, self.library_info),
                yield_batch_data(self.fq1, self.fq2, batch_size=self.batch_size))

            for stat_info in res:
                self.update_stat_info(stat_info)

        # 运行完成，清除状态文件
        if os.path.exists(self.validators.temp_state_file):
            os.remove(self.validators.temp_state_file)

    @runtime(f"{__name__}.collect_matrix")
    def collect_matrix(self):
        num_total = self.stat_info.total_num
        clean_num = self.stat_info.clean_num
        qual_counter_umi = self.stat_info.qual_counter_umi
        qual_counter_read = self.stat_info.qual_counter_read
        qual_counter_barcode = self.stat_info.qual_counter_barcode

        # print("clean_num", clean_num)
        # print("num_total", num_total)
        # print("checker4polyt", self.stat_info.num_for_no_polyt)
        # print("checker4barcode", self.stat_info.num_for_no_barcode)
        # print("checker4link", self.stat_info.num_for_no_link)
        # print("checker4low_qual", self.stat_info.num_for_low_qual)

        params_polyt = self.validators.params_polyt
        params_link = self.validators.params_link
        params_barcode = self.validators.params_barcode
        params_qual = self.validators.params_qual
        if params_polyt.use_polyt_valid_reads:
            logger.info(f"no polyT reads number : {self.stat_info.num_for_no_polyt}")
        else:
            logger.info(f"not valid polyT。")
        if params_link.use_link_valid_reads:
            logger.info(f"no_linker : {self.stat_info.num_for_no_link}")
        else:
            logger.info(f"not valid linker。")
        if params_barcode.use_barcode_valid_reads:
            logger.info(f"no_barcode : {self.stat_info.num_for_no_barcode}")
        else:
            logger.info(f"not valid barcode。")
        if params_qual.use_qual_valid_reads:
            logger.info(f"low qual reads number: {self.stat_info.num_for_low_qual}")
        else:
            logger.info(f"not valid qual。")
        if clean_num == 0:
            # TODO ensure param name `chemistry`?
            raise Exception(
                'no valid reads found! please check the --chemistry parameter.')

        BarcodesQ30 = sum([qual_counter_barcode[k] for k in qual_counter_barcode
                           if k >= chr(30 + 33)]) / float(sum(qual_counter_barcode.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        UMIsQ30 = sum([qual_counter_umi[k] for k in qual_counter_umi
                       if k >= chr(30 + 33)]) / float(sum(qual_counter_umi.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'

        ReadQ30 = sum(
            [qual_counter_read[k] for k in qual_counter_read
             if k >= chr(30 + 33)]
        ) / float(sum(qual_counter_read.values())) * 100
        ReadQ30 = round(ReadQ30, 2)
        ReadQ30_display = f'{ReadQ30}%'

        logger.info(f"Raw Reads {num_total}")
        logger.info(f"Valid Reads {round(clean_num / num_total * 100, 2)}%")
        logger.info(f"Q30 of Barcodes {BarcodesQ30}, display {BarcodesQ30_display}")
        logger.info(f"Q30 of Read {ReadQ30}, display {ReadQ30_display}")
        logger.info(f"Q30 of UMIs {UMIsQ30}, display {UMIsQ30_display}")
        self.add_metric(
            name='Number of Raw Reads',
            value=num_total,
            help_info='Total number of read pairs from Fastq file'
        )
        self.add_metric(
            name='Valid Library Reads',
            value=clean_num,
            total=num_total,
            help_info='Fraction of reads that conform to the expected library structure'
        )
        self.add_metric(
            name='Q30 of Barcode',
            value=BarcodesQ30,
            display=BarcodesQ30_display,
            help_info='Fraction of barcode bases with quality scores over Q30',
        )
        self.add_metric(
            name='Q30 of UMI',
            value=UMIsQ30,
            display=UMIsQ30_display,
            help_info='Fraction of UMI bases with quality scores over Q30',
        )
        self.add_metric(
            name='Q30 of RNA Read',
            value=ReadQ30,
            display=ReadQ30_display,
            help_info='Fraction of RNA read bases with quality scores over Q30',
        )
