#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : correct_barcode.py
@Time       : 2023/07/19
@Version    : 1.0
@Desc       : None
"""
import logging
import pandas as pd
from cell_cosmo.util.BarcodeCorrectUtil import ConstNS
from cell_cosmo.util.BarcodeCorrectUtil.BaseOut import BaseOut
from cell_cosmo.util.BarcodeCorrectUtil.correct_umis import CorrectUMIs
from cell_cosmo.util.BarcodeCorrectUtil.FindAllItdNbrPairs import FindAllItdNbrPairs
from cell_cosmo.util.BarcodeCorrectUtil.get_correct_dict import get_correct_dict

logger = logging.getLogger(__name__)


class CorrectBarcodeUMI:
    miss_pos = "miss%d"

    def __init__(self, out_util: BaseOut, df: pd.DataFrame, thread=1):
        self.full_df = df
        self.out_util = out_util
        self.barcode_len, self.umi_len = self._get_len()

        # TODO 差异的碱基数,当前算法中由于需要根据单个异变的位置,寻找itd-nbr对
        # todo 由于目前算法不支持2bp差异的itd-nbr对的搜索,所以暂时屏蔽相关参数的用户控制
        self.diff_num = 1
        self.thread = thread

    # 获取长度，并进行必要的验证
    def _get_len(self):
        # 对于full_data数据的入库对列名进行检查，以便后续操作都在预期的schema下进行
        columns = set(ConstNS.get_1_cols())
        barcode, umi = ConstNS.barcode, ConstNS.umi
        for col in self.full_df.columns.tolist():
            assert col in columns, f"the column {col} not in file header!"
        # 对 barcode 和 umi length 进行检查
        res = self.full_df[barcode].map(lambda x: len(x)).unique()
        assert res.size == 1, f"the barcode has inconsistent length,please check it!!"
        barcode_len = int(res[0])  # todo 5<=len<=30,是否需要验证一下？
        res = self.full_df[umi].map(lambda x: len(x)).unique()
        assert res.size == 1, f"the umi has inconsistent length,please check it!!"
        umi_len = int(res[0])
        return barcode_len, umi_len

    def get_barcode_matrix(self, n_umi_filter=20):
        barcode, count = ConstNS.barcode, ConstNS.count  # 将列名提取出来，避免后续修改出错
        barcode_len = self.barcode_len
        stat_df = self.full_df.groupby(barcode, as_index=False).agg({count: "sum"})
        # 按需求逻辑， barcode umi count <= 20,这些数据将被去除
        le_df = stat_df[stat_df[count] <= n_umi_filter]
        df = stat_df[stat_df[count] > n_umi_filter].copy()
        self.out_util.to_csv(f"barcode_umi_less_than_{n_umi_filter}_filter_{le_df.shape[0]}_data.tsv", le_df)

        logger.info(f"there is {le_df.shape[0]} barcode umi less than {n_umi_filter}")

        # 构建缺失列，以便后续分析
        remove_cols = [0, barcode_len + 1]
        # 过滤0
        ndf = pd.concat([df, df[barcode].str.split('', expand=True)], axis=1)
        ndf.drop(remove_cols, axis=1, inplace=True)
        # 将缺失序列提前构建
        for i in range(1, barcode_len + 1):
            cols = [c for c in range(1, barcode_len + 1) if c != i]
            ndf[self.miss_pos % i] = ndf[list(cols)].apply("".join, axis=1)
        return ndf


def _correct_barcode(full_df: pd.DataFrame, nbr2itd: dict) -> pd.DataFrame:
    def _map_correct(x):
        return nbr2itd.get(x, x)

    barcode = ConstNS.barcode
    gene_id = ConstNS.gene_id
    umi = ConstNS.umi
    count = ConstNS.count

    logger.info("start to correct barcode ... ")
    # 该日志必须在full_df 没有修改时打印
    logger.info(f"before barcode correct, n_total/n_barcode is "
                f"{full_df.shape[0]}/{full_df[barcode].unique().shape[0]}")
    full_df[barcode] = full_df[barcode].map(_map_correct)
    group_by_cols = [barcode, gene_id, umi]
    full_data2 = full_df.groupby(group_by_cols)[count].sum().reset_index()
    logger.info(f"after barcode correct, n_total/n_barcode is "
                f"{full_data2.shape[0]}/{full_data2[barcode].unique().shape[0]}")

    logger.info("correct barcode complete! ")
    return full_data2


def correct_barcode_and_umi(
        df: pd.DataFrame, outdir: str, sample: str,
        n_umi_filter=20, filter_limit=0.01, percent=0.1,
        thread=4
) -> pd.DataFrame:
    # 使用文件数据库来处理，避免多进程间的数据传递，同时避免输出更多临时数据到文件中
    # todo 输入参数已经修改了，注意修改， 另外添加一些参数的用户控制
    out_prefix = f"{outdir}/{sample}"
    out_util = BaseOut(outdir, sample)
    cbu = CorrectBarcodeUMI(out_util, df)
    barcode_df = cbu.get_barcode_matrix(n_umi_filter=n_umi_filter)
    # todo 如果 barcode_df 没有数据怎么办法？
    barcode = ConstNS.barcode
    full_df = df.loc[df[barcode].isin(barcode_df[barcode]), :].copy()

    # 构造barcode的缺失序列
    finder = FindAllItdNbrPairs(outdir, sample, barcode_df, barcode_len=cbu.barcode_len, thread=thread)
    itd_nbr_df = finder.search_to_db()
    # # 数据输出到文件
    # itd_nbr_df.to_csv(f"{outdir}/{sample}_barcode_itd_nbr.tsv", sep="\t", index=False)

    # 过滤无效的itd-nbr对，得到校正字典
    nbr2itd = get_correct_dict(itd_nbr_df, out_util, filter_limit=filter_limit)
    barcode_correct_df = _correct_barcode(full_df, nbr2itd)

    all_correct_df = CorrectUMIs(barcode_correct_df,
                                 out_util=out_util,

                                 percent=percent,
                                 thread=thread).correct()
    return all_correct_df


if __name__ == '__main__':
    _out = "/mnt/data/cellcosmo/xxx"
    _df = pd.read_csv(f"{_out}/test.csv")
    _rdf = correct_barcode_and_umi(_df, _out, "s1", n_umi_filter=5)
    _rdf.to_csv(f"{_out}/abc.csv", sep="\t")
