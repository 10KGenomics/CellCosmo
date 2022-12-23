#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : cutadapt_test.py
@Time       : 2022/06/06
@Version    : 1.0
@Desc       : None
"""
import re
import logging
import subprocess
from cell_cosmo.util import reader, runtime
from cell_cosmo.output_runner import BaseReportRunner

ADAPTER = ['polyT=A{18}', 'p5=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC']
logger = logging.getLogger(__name__)

LOG_METRICS_TITLE = (
    'Total reads processed',
    'Reads with adapters',
    'Reads that were too short',
    'Reads written (passing filters)',
    'Total basepairs processed',
    'Quality-trimmed',
    'Total written (filtered)',
)


def read_cutadapt_log(cutadapt_log):
    """
    Returns:
        dict: {
            'Total reads processed': int, 'Reads with adapters': int, 'Reads that were too short': int,
            'Reads written (passing filters)': int, 'Total base pairs processed': int,
            'Quality-trimmed': int, 'Total written (filtered)': int
        }
    """
    metrics_dict = {}
    remove_strs = [r',', r' bp', r'\(.*\)']

    for _line_index, line in enumerate(cutadapt_log.split('\n')):
        line = line.strip()
        if not line:
            continue
        attr = line.split(":")
        if attr[0] in LOG_METRICS_TITLE:
            title, number = attr[0], attr[1]
            number = attr[1].strip()
            for remove_str in remove_strs:
                number = re.sub(pattern=remove_str, repl='', string=number)
            metrics_dict[title] = int(number)

    return metrics_dict


class Cutadapt(BaseReportRunner):
    _STEP_NAME = "cutadapt"
    _DISPLAY_TITLE = "Trimming"

    def __init__(self, minimum_length, nextseq_trim, overlap, insert, cutadapt_param, fq, **kwargs):
        super(Cutadapt, self).__init__(**kwargs)
        # TODO 需要由报告对象管理的参数
        # self.outdir = kwargs.get("outdir")
        # self.sample = kwargs.get("sample")
        # self.assay = kwargs.get("assay")  # ??
        # self.thread = kwargs.get("thread")
        # self.debug = kwargs.get("debug")
        # self.out_prefix = f'{self.outdir}/{self.sample}'
        # self.display_title = kwargs.get("display_title")
        # -------------------------------
        # TODO 待验证的入参
        self.minimum_length = minimum_length
        self.nextseq_trim = nextseq_trim
        self.overlap = overlap
        self.insert = insert
        self.cutadapt_param = cutadapt_param
        self.fq = fq
        # --------------------------
        self.adapter_args = self.read_adapter_fasta(kwargs.get("adapter_fasta"))
        self.adapter_args += ADAPTER
        suffix = ".gz" if kwargs.get("gzip", False) else ""
        self.out_fq2 = f'{self.outdir}/{self.sample}_clean_2.fq{suffix}'
        self.cutadapt_log = None
        self.cutadapt_log_file = f'{self.outdir}/cutadapt.log'

    def read_adapter_fasta(self, adapter_fasta):
        """

        :param adapter_fasta:
        :return: ['adapter1=AAA','adapter2=BBB']
        """
        adapter_args = []
        if adapter_fasta and adapter_fasta != 'None':
            for name, seq in reader(adapter_fasta):
                adapter_args.append(f'{name}={seq}')
        return adapter_args

    @runtime(__name__)
    def run(self):
        adapter_args_str = " ".join(['-a ' + adapter for adapter in self.adapter_args])
        cmd = (
            'cutadapt '
            f'{adapter_args_str} '
            f'-n {len(self.adapter_args)} '
            f'-j {self.thread} '
            f'-m {self.minimum_length} '
            f'--nextseq-trim={self.nextseq_trim} '
            f'--overlap {self.overlap} '
            f'-l {self.insert} '
            f'-o {self.out_fq2} '
            f'{self.cutadapt_param} '
            f'{self.fq} '
        )
        logger.info(cmd)
        # need encoding argument to return str
        results = subprocess.run(
            cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE,
            encoding='utf-8', check=True, shell=True
        )
        self.cutadapt_log = results.stdout
        with open(self.cutadapt_log_file, 'w') as f:
            f.write(self.cutadapt_log)

    def collect_matrix(self):
        # Total reads processed:...Total written (filtered):
        metrics_dict = read_cutadapt_log(self.cutadapt_log)

        total_reads = metrics_dict['Total reads processed']
        reads_with_adapters = metrics_dict['Reads with adapters']
        reads_too_short = metrics_dict['Reads that were too short']
        reads_written = metrics_dict['Reads written (passing filters)']
        total_base_pairs = metrics_dict['Total basepairs processed']
        quality_trimmed = metrics_dict['Quality-trimmed']
        base_pairs_written = metrics_dict['Total written (filtered)']

        self.add_metric(
            name='RNA Reads with Adapters',
            value=reads_with_adapters,
            total=total_reads,
            help_info='reads number with sequencing adapters or with poly A'
        )
        self.add_metric(
            name='Too Short Reads',
            value=reads_too_short,
            total=total_reads,
            help_info='reads number with length less than 20bp after trimming'
        )
        self.add_metric(
            name='Reads Written(passing filters)',
            value=reads_written,
            total=total_reads,
            help_info='reads number pass filtering'
        )
        self.add_metric(
            name='Total Base pairs Processed',
            value=total_base_pairs,
            help_info='total processed base pairs'
        )
        self.add_metric(
            name='Base Pairs Quality-Trimmed',
            value=quality_trimmed,
            total=total_base_pairs,
            help_info='Base pairs removed from the end of the read whose '
                      'quality is smaller than the given threshold'
        )
        self.add_metric(
            name='Base Pairs Written(filtered)',
            value=base_pairs_written,
            total=total_base_pairs,
            help_info='base pairs pass filtering'
        )
