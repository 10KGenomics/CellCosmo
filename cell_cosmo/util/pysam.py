#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : pysam.py
@Time       : 2022/06/10
@Version    : 1.0
@Desc       :
这个文件并非为了替代pysam而开发的高效代码,
本代码仅仅为了使本项目在windows下运行时可以正常调试运行
"""


class xAMRowObj:
    """sam or bam one line obj"""

    def __init__(self, line: str):
        if type(line) == bytes:
            line = bytes.decode(line)
        (
            self._query_name,  # 比对的序列名称
            self._flag,  # Bwise FLAG（表明比对类型：paring，strand，mate strand等） 例如：99
            self.rename,  # 比对上的参考序列名 例如：NC_000075.6
            self.pos,  # 比对上的参考序列名 例如：NC_000075.6
            self.MAPQ,  # MAPQ 比对质量 例如：60
            self.CIGAR,  # Extended CIGAR string（操作符：MIDNSHP）比对结果信息；匹配碱基数，可变剪接等 例如：87M
            self.MRNM,  # MRNM 相匹配的另外一条序列，比对上的参考序列名 例如：=
            self.MPOS,  # MPOS 1-Based leftmost Mate Position （相比于MRNM列来讲意思和POS差不多） 例如：124057667
            self.ISIZE,  # ISIZE 插入片段长度 例如：200
            self.SEQ,  # SEQ 和参考序列在同一个链上比对的序列（若比对结果在负义链上，则序列是其反向重复序列，反向互补序列） 例如：ATTACTTGGCTGCT
            self.QUAL,  # QUAL 比对序列的质量（ASCII-33=Phred base quality）reads碱基质量值 例如：-8CCCGFCCCF7@E-
            *self.tags
        ) = line.split("\t")
        self.name = self._query_name
        self._tag_dict = [t.split(":", maxsplit=1) for t in self.tags]

    @property
    def query_name(self) -> str:
        return self._query_name

    def set_tag(self, tag, value, value_type):
        pass

    def has_tag(self, tag):
        return tag in self._tag_dict

    def get_tag(self, tag):
        assert self.has_tag(tag)
        return self._tag_dict[tag]


class AlignmentFile:
    def __init__(self, bam, mode="r", header=None):
        # 需要打开一个文件流
        self.bam = bam
        self.header = header
        self._fh = open(bam, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __iter__(self):
        # 需要从文件流中1行1行读取
        yield xAMRowObj(self._fh.readline())

    def write(self, line: xAMRowObj):
        pass


def main():
    from itertools import groupby
    sam_file = AlignmentFile(r"E:\result\20220303HXY6_name_sorted.bam", "rb")
    for g, h in groupby(sam_file):
        print(g, h)
    with AlignmentFile("in_bam", "rb") as original_bam:
        header = original_bam.header
        with AlignmentFile("temp_sam_file", "w", header=header) as temp_sam:
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
