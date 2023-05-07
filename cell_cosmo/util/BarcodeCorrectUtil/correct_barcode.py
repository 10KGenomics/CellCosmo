#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : correct_barcode.py
@Time       : 2022/07/13
@Version    : 1.0
@Desc       : None
"""
import pandas as pd
from cell_cosmo.util.BarcodeCorrectUtil.db_util import DBUtil
from cell_cosmo.util.BarcodeCorrectUtil.find_all_itd_nbr_pairs import FindAllItdNbrPairs
from cell_cosmo.util.BarcodeCorrectUtil.correct_umis import CorrectUMIs
from cell_cosmo.util.BarcodeCorrectUtil.get_correct_dict import get_correct_dict
import logging

logger = logging.getLogger(__name__)


def _correct_barcode(db_util: DBUtil, full_df: pd.DataFrame, nbr2itd: dict) -> pd.DataFrame:
    def _map_correct(x):
        return nbr2itd.get(x, x)

    logger.info("start to correct barcode ... ")
    # 该日志必须在full_df 没有修改时打印
    logger.info(f"before barcode correct, n_total/n_barcode is "
                f"{full_df.shape[0]}/{full_df[DBUtil.barcode].unique().shape[0]}")
    full_df[DBUtil.barcode] = full_df[DBUtil.barcode].map(_map_correct)
    group_by_cols = [DBUtil.barcode, DBUtil.gene_id, DBUtil.umi]
    full_data2 = full_df.groupby(group_by_cols)[DBUtil.count].sum().reset_index()
    # 将barcode校正后的detail_count数据入库保存
    DBUtil.CountDetailAll2.df2db(db_util, full_data2)

    logger.info(f"after barcode correct, n_total/n_barcode is "
                f"{full_data2.shape[0]}/{full_data2[DBUtil.barcode].unique().shape[0]}")

    logger.info("correct barcode complete! ")
    return full_data2


def correct_barcode_and_umi(
        df: pd.DataFrame, outdir: str, sample: str,
        n_umi_filter=20, filter_limit=0.01, percent=0.1,
        thread=4
) -> pd.DataFrame:
    # 使用文件数据库来处理，避免多进程间的数据传递，同时避免输出更多临时数据到文件中
    out_prefix = f"{outdir}/{sample}"  # todo 输入参数已经修改了，注意修改， 另外添加一些参数的用户控制
    db_util = DBUtil(outdir, sample)
    # 将输入数据入库
    DBUtil.CountDetailAll.df2db(db_util, df)
    # 检查
    barcode_len, umi_len = DBUtil.CountDetailAll.check_df(db_util, df)

    barcode_df = DBUtil.CountDetailAll.get_barcode_matrix(
        db_util,
        barcode_len,
        n_umi_filter=n_umi_filter
    )
    full_df = df.loc[df[DBUtil.barcode].isin(barcode_df[DBUtil.barcode]), :].copy()
    # 构造barcode的缺失序列
    finder = FindAllItdNbrPairs(db_util, barcode_df, barcode_len=barcode_len, thread=thread)
    finder.search_to_db()
    # 过滤无效的itd-nbr对，得到校正字典
    nbr2itd = get_correct_dict(db_util, filter_limit=filter_limit)

    barcode_correct_df = _correct_barcode(db_util, full_df, nbr2itd)
    all_correct_df = CorrectUMIs(db_util, barcode_correct_df,
                                 output=f"{out_prefix}_count_detail.tsv",
                                 percent=percent,
                                 thread=thread).correct()
    return all_correct_df
