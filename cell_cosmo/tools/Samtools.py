#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : Samtools.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import logging
import subprocess
import pysam
from cell_cosmo.util import runtime
from cell_cosmo.util.GTFDictUtil import GTFDict

logger = logging.getLogger(__name__)


class Samtools():
    def __init__(self, in_bam, out_bam, threads=1, debug=False):
        self.in_bam = in_bam
        self.out_bam = out_bam
        self.threads = threads
        self.temp_sam_file = f"{self.out_bam}_sam.temp"
        self.debug = debug

    @runtime(f"{__name__}.sort")
    def samtools_sort(self, in_file, out_file, by='coord'):
        cmd = f"samtools sort {in_file} -o {out_file} --threads {self.threads}"
        if by == "name":
            cmd += " -n"
        logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    @runtime(f"{__name__}.index")
    def samtools_index(self, in_file):
        cmd = f"samtools index {in_file}"
        logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    def sort_bam(self, by='coord'):
        """sort in_bam"""
        self.samtools_sort(self.in_bam, self.out_bam, by=by)

    def index_bam(self):
        """index out_bam"""
        self.samtools_index(self.out_bam)

    @runtime(f"{__name__}.add_tag")
    def add_tag(self, gtf_file):
        """
        - CB cell barcode
        - UB UMI
        - GN gene name
        - GX gene id
        """
        gtf_dict = GTFDict(gtf_file)
        save = pysam.set_verbosity(0)
        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            pysam.set_verbosity(save)
            header = original_bam.header
            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    attr = read.query_name.split('_')
                    barcode = attr[0]
                    umi = attr[1]
                    read.set_tag(tag='CB', value=barcode, value_type='Z')
                    read.set_tag(tag='UB', value=umi, value_type='Z')
                    # assign to some gene
                    if read.has_tag('XT'):
                        gene_id = read.get_tag('XT')
                        # if multi-mapping reads are included in original bam,
                        # there are multiple gene_ids
                        if ',' in gene_id:
                            gene_name = [gtf_dict[i] for i in gene_id.split(',')]
                            gene_name = ','.join(gene_name)
                        else:
                            gene_name = gtf_dict[gene_id]
                        read.set_tag(tag='GN', value=gene_name, value_type='Z')
                        read.set_tag(tag='GX', value=gene_id, value_type='Z')
                    temp_sam.write(read)

    @runtime(f"{__name__}.add_rg")
    def add_RG(self, barcodes):
        """
        barcodes list
        """
        with pysam.AlignmentFile(self.in_bam, "rb") as original_bam:
            header = original_bam.header.to_dict()
            header['RG'] = []
            for index, barcode in enumerate(barcodes):
                header['RG'].append({
                    'ID': barcode,
                    'SM': index + 1,
                })

            with pysam.AlignmentFile(self.temp_sam_file, "w", header=header) as temp_sam:
                for read in original_bam:
                    read.set_tag(tag='RG', value=read.get_tag('CB'), value_type='Z')
                    temp_sam.write(read)

    def temp_sam2bam(self, by=None):
        self.samtools_sort(self.temp_sam_file, self.out_bam, by=by)
        self.rm_temp_sam()

    def rm_temp_sam(self):
        cmd = f"rm {self.temp_sam_file}"
        subprocess.check_call(cmd, shell=True)
