#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Pipeline.py
@Time       : 2022/06/15
@Version    : 1.0
@Desc       : None
"""
import os
import time
import logging
import subprocess
from cell_cosmo.util import GenomeUtil, PathUtil
from cell_cosmo.rna.PipelineMeta import (
    MetaBarcode as S2,
    MetaCount as S6,
    MetaCutAdapt as S3,
    MetaFeatureCounts as S5,
    MetaMkRef as S0,
    MetaSample as S1,
    MetaStar as S4,
    MetaAnalysis as S7,
)
from cell_cosmo.rna.PipelineMeta.CMDMeta import CMDMeta, CMDBase, PipeLineConfigParser

logger = logging.getLogger(__name__)


class Pipeline:
    def __init__(self, config_path):
        self._cmds = []
        self.cfg = PipeLineConfigParser(config_path)
        # 将输出放到下一层,使真正的输出目录(上一层)中可以输出html
        self.outdir = os.path.join(self.cfg.outdir, "result")
        self.sample = self.cfg.sample
        self.genomeDir = self.cfg.genomeDir
        self._init_cmdline()

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    def _add(self, cb: CMDBase):
        is_test_env = os.getenv('CELL_COSMO_TEST_ENV_CMD_TEST_MODE', False)
        cmd = cb.getcmd(self.cfg)
        if is_test_env and str(is_test_env).upper() == 'TEST':
            self._cmds.append('python -m ' + cmd)
        else:
            self._cmds.append(cmd)

    def _init_cmdline(self):
        cfg = self.cfg
        # 确定是否需要构建索引
        if not os.path.exists(os.path.join(cfg.genomeDir, GenomeUtil.GENOME_CONFIG)):
            logger.info(f"the pipeline will run star index first,and genomeDir is: {cfg.genomeDir}")
            self._cmds.append(S0().getcmd(cfg))

        common = dict(outdir=self.outdir, sample=self.sample)
        stp1 = S1(fq1=cfg.fq1, **common)
        self._add(stp1)

        stp2 = S2(fq1=cfg.fq1, fq2=cfg.fq2, **common)
        self._add(stp2)

        stp3 = S3(fq=self._get_file_path("_2.fq", stp2.gzip), **common)
        self._add(stp3)

        common.update(genomeDir=self.genomeDir)
        stp4 = S4(fq=self._get_file_path("_clean_2.fq", stp3.gzip), **common)
        self._add(stp4)

        stp5 = S5(inp=self._get_file_path("_Aligned.sortedByCoord.out.bam", None), **common)
        self._add(stp5)
        stp6 = S6(bam=self._get_file_path("_name_sorted.bam", None), **common)
        self._add(stp6)
        stp7 = S7(matrix_file=self._get_file_path("_filtered_feature_bc_matrix", None), **common)
        self._add(stp7)

    def run(self):
        for cmd in self._cmds:

            logger.info(f"Run: {cmd}")
            if "CellCosmo rna mkref" in str(cmd):
                # 此时需要切换目录:
                with PathUtil.chdir(self.genomeDir):
                    subprocess.check_call(cmd, shell=True)
            else:
                subprocess.check_call(cmd, shell=True)

    def _get_file_path(self, sub_name, gzip):
        suffix = ".gz" if gzip else ""
        return os.path.join(self.outdir, f"{self.sample}{sub_name}{suffix}")
