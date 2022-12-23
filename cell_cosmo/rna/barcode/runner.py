#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : barcode_ck.py
@Time       : 2022/05/26
@Version    : 1.0
@Desc       : None
"""
from cell_cosmo.output_runner import BaseReportRunner
from cell_cosmo.rna.barcode.splitter import BarcodeSplitter
from cell_cosmo.util import runtime, fmt_number, reader
import logging

logger = logging.getLogger(__name__)


class BarcodeSplitterRunner(BaseReportRunner):
    _STEP_NAME = "barcode"
    _DISPLAY_TITLE = "Demultiplexing"

    def __init__(self, fq1, fq2, pattern=None, chemistry_config=None, chemistry_name=None, **kwargs):
        super(BarcodeSplitterRunner, self).__init__(**kwargs)
        if all(x is None for x in [pattern, chemistry_config, chemistry_name]):
            raise Exception("Please specify parameter chemistry_name, "
                            "or specify parameters pattern and chemistry_config path")

        self.barcode_splitter = BarcodeSplitter(
            chemistry_name=chemistry_name,
            chemistry_config=chemistry_config,
            pattern=pattern, **kwargs)

        self.fq1_list = fq1.split(",")
        self.fq2_list = fq2.split(",")
        self.fq_num = len(self.fq1_list)
        if self.fq_num != len(self.fq2_list):
            raise Exception('fq1 and fq2 must be the same file number!')

    def _run(self, fq1, fq2):
        for e1, e2 in zip(reader(fq1), reader(fq2)):
            # header1, seq1, _, qual1 = e1
            # header2, seq2, _, qual2 = e2
            self.barcode_splitter.process_reads_pair(e1, e2)

    @runtime(__name__)
    def run(self):
        for i in range(self.fq_num):
            self._run(self.fq1_list[i], self.fq2_list[i])
            logger.info(f"{self.fq1_list[i]} finished.")
        self.barcode_splitter.close_all()

    @runtime(f"{__name__}.collect_matrix")
    def collect_matrix(self):
        bs = self.barcode_splitter
        num_total = bs.num_total
        clean_num = bs.num
        logger.info(f"processed reads: {fmt_number(num_total)}. "
                    f"valid reads: {fmt_number(bs.num)}.")
        logger.info(f"no polyT reads number : {bs.checker4polyt.num}")
        logger.info(f"no_linker : {bs.checker4link.num}")
        logger.info(f"no_barcode : {bs.checker4barcode.num}")
        logger.info(f"low qual reads number: {bs.checker4qual.num}")
        if clean_num == 0:
            # TODO ensure param name `chemistry`?
            raise Exception(
                'no valid reads found! please check the --chemistry parameter.')

        BarcodesQ30 = sum([bs.qual_counter_barcode[k] for k in bs.qual_counter_barcode
                           if k >= bs.checker4qual.ord2chr(
                30)]) / float(sum(bs.qual_counter_barcode.values())) * 100
        BarcodesQ30 = round(BarcodesQ30, 2)
        BarcodesQ30_display = f'{BarcodesQ30}%'
        UMIsQ30 = sum([bs.qual_counter_umi[k] for k in bs.qual_counter_umi
                       if k >= bs.checker4qual.ord2chr(
                30)]) / float(sum(bs.qual_counter_umi.values())) * 100
        UMIsQ30 = round(UMIsQ30, 2)
        UMIsQ30_display = f'{UMIsQ30}%'

        ReadQ30 = sum(
            [bs.qual_counter_read[k] for k in bs.qual_counter_read
             if k >= bs.checker4qual.ord2chr(30)]
        ) / float(sum(bs.qual_counter_read.values())) * 100
        ReadQ30 = round(ReadQ30, 2)
        ReadQ30_display = f'{ReadQ30}%'

        logger.info(f"Raw Reads {bs.num_total}")
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

