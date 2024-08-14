# #!/usr/bin/env python
# # -*- coding:utf-8 -*-
# """
# @Author     : ice-melt@outlook.com
# @File       : base_starsolo_runner.py
# @Time       : 2024/07/09
# @Version    : 1.0
# @Desc       : None
# """
# import re
# import logging
# import subprocess
# from cell_cosmo.util import runtime
# from cell_cosmo.tools import utils
# from cell_cosmo.util import GenomeUtil
# from cell_cosmo.output_runner import BaseReportRunner
# from cell_cosmo.tools.chemistry import LibraryInfo
#
# logger = logging.getLogger(__name__)
#
#
# def get_solo_pattern(library_info: LibraryInfo) -> (str, str, str):
#     """
#     Returns:
#         solo_type
#         cb_str
#         umi_str
#     """
#     if len(library_info.pattern_U) != 1:
#         raise Exception(f"Error: Wrong pattern:{library_info.pattern_str}. \n "
#                         f"Solution: fix pattern so that UMI only have 1 position.\n")
#
#     ul, ur = library_info.pattern_U[0]
#     umi_len = ur - ul
#
#     if len(library_info.pattern_C) == 1:
#         solo_type = 'CB_UMI_Simple'
#         l, r = library_info.pattern_C[0]
#         cb_start = l + 1
#         cb_len = r - l
#         umi_start = ul + 1
#         cb_str = f'--soloCBstart {cb_start} --soloCBlen {cb_len} '
#         umi_str = f'--soloUMIstart {umi_start} --soloUMIlen {umi_len} '
#     else:
#         solo_type = 'CB_UMI_Complex'
#         cb_pos = ' '.join([f'0_{s}_0_{e - 1}' for s, e in library_info.pattern_C])
#         umi_pos = f'0_{ul}_0_{ur - 1}'
#         cb_str = f'--soloCBposition {cb_pos} '
#         umi_str = f'--soloUMIposition {umi_pos} --soloUMIlen {umi_len} '
#     return solo_type, cb_str, umi_str
#
#
# class BaseSTARRunner(BaseReportRunner):
#     """
#     base class for STAR
#     """
#     _STEP_NAME = None
#     _DISPLAY_TITLE = None
#
#     def __init__(self, genomeDir, fq1, fq2, out_unmapped,
#                  outFilterMatchNmin, consensus_fq,
#                  STAR_param, samtools_index_param=None, **kwargs):
#         # 删除了 outFilterMultimapNmax 参数
#         super(BaseSTARRunner, self).__init__(**kwargs)
#
#         self.fq1 = fq1
#         self.fq2 = fq2
#         self.chemistry_name = kwargs.get("chemistry_name", "default")
#         self.library_info = LibraryInfo(
#             library_name=self.chemistry_name,
#             auto_init_db=False
#         )
#         self.solo_type, self.cb_str, self.umi_str = get_solo_pattern(self.library_info)
#         self.whitelist_str = " ".join(self.library_info.barcode_files)
#
#         self.genomeDir = genomeDir
#         self.out_unmapped = out_unmapped
#
#         self.soloFeatures = kwargs.get("soloFeatures")
#         self.SAM_attributes = kwargs.get("SAM_attributes")
#         self.soloCellFilter = kwargs.get("soloCellFilter")
#
#         self.outFilterMatchNmin = int(outFilterMatchNmin)
#         # self.multi_max = int(outFilterMultimapNmax)
#         self.STAR_param = STAR_param
#         self.samtools_index_param = samtools_index_param
#         self.consensus_fq = consensus_fq
#
#         # parse
#         self.genome = GenomeUtil.parse_dir(self.genomeDir)
#         self.stat_prefix = 'Reads'
#         if self.consensus_fq:
#             self.stat_prefix = 'UMIs'
#
#         # out
#         self.outPrefix = f'{self.outdir}/{self.sample}_'
#         self.STAR_map_log = f'{self.outPrefix}Log.final.out'
#         self.unsort_STAR_bam = f'{self.outPrefix}Aligned.out.bam'
#         self.STAR_bam = f'{self.outPrefix}Aligned.sortedByCoord.out.bam'
#
#         # output files
#         self.solo_out_dir = f'{self.outdir}/{self.sample}_Solo.out/'
#         solo_dir = f'{self.outdir}/{self.sample}_Solo.out/GeneFull_Ex50pAS'
#         self.raw_matrix = f'{solo_dir}/raw'
#         self.filtered_matrix = f'{solo_dir}/filtered'
#         self.summary_file = f'{solo_dir}/Summary.csv'
#         self.outs = [self.raw_matrix, self.filtered_matrix, self.STAR_bam]
#
#     @runtime(f'{__name__}.STAR')
#     def STAR(self):
#
#         cmd = [
#             'STAR',
#             '--genomeDir', self.genomeDir,
#             '--readFilesIn', self.fq2 + ' ' + self.fq1,
#             f'--soloCBwhitelist {self.whitelist_str}',
#             f'--soloCellFilter {self.soloCellFilter}',
#             '--outFileNamePrefix', self.outPrefix,
#             '--runThreadN', str(self.thread),
#             '--outFilterMatchNmin', str(self.outFilterMatchNmin),
#             f'--soloFeatures {self.soloFeatures}',
#             f'--outSAMattributes {self.SAM_attributes}',
#             f'--soloType {self.solo_type}',
#             self.cb_str,  # https://github.com/alexdobin/STAR/issues/1607
#             self.umi_str,
#             # length of the barcode read 1 equal to sum of soloCBlen+soloUMIlen 0 not defined, do not check
#             # only one match in whitelist with 1 mismatched base allowed
#             '--soloCBmatchWLtype 1MM',
#             # controls sort by Coordinate or not,(Unsorted，SortedByCoordinate)
#             '--outSAMtype BAM SortedByCoordinate',
#             # cell reads stats https://github.com/alexdobin/STAR/issues/1501
#             '--soloCellReadStats Standard',
#             '--soloBarcodeReadLength 0',
#             # '--outFilterMultimapNmax', str(self.multi_max),
#             # '--soloBarcodeMate', '0',
#             # '--clip5pNbases','74 0',
#
#         ]
#         if self.fq1[-3:] == ".gz":
#             # 苹果电脑 zcat 异常
#             cmd += ['--readFilesCommand', 'zcat']
#         # if self.out_unmapped:
#         #     cmd += ['--outReadsUnmapped', 'Fastx']
#         print(cmd)
#         cmd = ' '.join(cmd)
#         if self.STAR_param:
#             # STAR_param 异常添加双引号或单引号的处理
#             self.STAR_param = str(self.STAR_param).lstrip(
#                 "'").lstrip('"').rstrip("'").rstrip('"')
#             cmd += (" " + self.STAR_param)
#         logger.info(cmd)
#         subprocess.check_call(cmd, shell=True)
#
#     def run(self):
#         self.STAR()
#         # TODO 数据收集
#         # self.get_star_metrics()
#         # self.sort_bam()
#         self.index_bam()
#
#     @runtime(f'{__name__}.sort_bam')
#     def sort_bam(self):
#         utils.sort_bam(
#             self.unsort_STAR_bam,
#             self.STAR_bam,
#             threads=self.thread,
#         )
#
#     @runtime(f'{__name__}.index_bam')
#     def index_bam(self):
#         utils.index_bam(self.STAR_bam, self.samtools_index_param)
#
#     def collect_matrix(self):
#         """collect matrix"""
#
#         with open(self.STAR_map_log, 'r') as map_log:
#             # number amd percent
#             unique_reads_list = []
#             multi_reads_list = []
#             total_reads = 0
#             for line in map_log:
#                 if line.strip() == '':
#                     continue
#                 if re.search(r'Uniquely mapped reads', line):
#                     unique_reads_list.append(line.strip().split()[-1])
#                 if re.search(r'of reads mapped to too many loci', line):
#                     multi_reads_list.append(line.strip().split()[-1])
#                 if re.search(r'Number of input reads', line):
#                     total_reads = int(line.strip().split()[-1])
#
#         unique_reads = int(unique_reads_list[0])
#         multi_reads = int(multi_reads_list[0])
#
#         self.add_metric(
#             name='Genome',
#             value=self.genome['genome_name'],
#         )
#         self.add_metric(
#             name=f'Uniquely Mapped {self.stat_prefix}',
#             value=unique_reads,
#             total=total_reads,
#             # help_info='reads that mapped uniquely to the genome'
#             help_info='Fraction of reads that mapped uniquely to the genome'
#         )
#         self.add_metric(
#             name=f'Multi-Mapped {self.stat_prefix}',
#             value=multi_reads,
#             total=total_reads,
#             # help_info='reads that mapped to multiple locations in the genome'
#             help_info='Fraction of Reads that mapped to multiple locations in the genome'
#         )
