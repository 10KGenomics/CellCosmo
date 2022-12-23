#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : find_all_itd_nbr_pairs.py
@Time       : 2022/07/07
@Version    : 1.0
@Desc       : None
"""
import itertools
import pandas as pd
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
# import Levenshtein
from cell_cosmo.util.BarcodeCorrectUtil.db_util import DBUtil
import logging

logger = logging.getLogger(__name__)


def multiprocess_callback(pbar: tqdm, res):
    pbar.update()


class FindAllItdNbrPairs:
    # modify: del params: diff_num=1, barcode_diff_limit=0.1,
    def __init__(self, db_util: DBUtil, df: pd.DataFrame, barcode_len: int, thread=4):
        self.db_util = db_util
        self.barcode_len = barcode_len  # TODO 超参数,1:barcode的长度都需要保持一致;2: 5<=len<=30
        self.df = df

        # TODO 差异的碱基数,当前算法中由于需要根据单个异变的位置,寻找itd-nbr对
        # todo 由于目前算法不支持2bp差异的itd-nbr对的搜索,所以暂时屏蔽相关参数的用户控制
        self.diff_num = 1
        self.thread = thread

    def write_error(self, e):
        DBUtil.ErrorLog.write_log(self.db_util, msg="[barcode correct error]:" + str(e))

    def _find_synthesis_error_pair(self, pos):
        """寻找合成错误"""
        # 将列名提取出来，避免后续修改出错

        miss_fmt = DBUtil.miss_pos
        miss_col = miss_fmt % pos
        miss_col_last = miss_fmt % self.barcode_len

        df1 = self.df[[DBUtil.barcode, DBUtil.count, miss_col, pos]].rename({
            DBUtil.barcode: DBUtil.intended_barcode,
            DBUtil.count: DBUtil.intended_size,
            pos: DBUtil.intended_base
        }, axis=1).set_index(miss_col)
        df2 = self.df[[DBUtil.barcode, DBUtil.count, miss_col_last]].rename({
            DBUtil.barcode: DBUtil.neighbor_barcode,
            DBUtil.count: DBUtil.neighbor_size,
        }, axis=1).set_index(miss_col_last)

        df3 = df1.join(df2, how='inner')
        df4 = df3[df3[DBUtil.intended_barcode] != df3[DBUtil.neighbor_barcode]].copy()
        df4[DBUtil.position] = pos
        df4[DBUtil.neighbor_base] = "-"
        df4 = df4[[
            DBUtil.intended_barcode,
            DBUtil.neighbor_barcode,
            DBUtil.intended_size,
            DBUtil.neighbor_size,
            DBUtil.position,
            DBUtil.intended_base,
            DBUtil.neighbor_base
        ]]
        DBUtil.ItdNbr.df2db(self.db_util, df4)

    def _find_substitution_error_pair(self, pos):
        # 将列名提取出来，避免后续修改出错
        miss_fmt = DBUtil.miss_pos
        miss_col = miss_fmt % pos
        miss_col_last = miss_fmt % self.barcode_len
        base_idx = int(pos) - 1
        result = []
        for key, ddf in self.df.groupby(miss_col):
            if ddf.shape[0] <= 1:
                # 没有与其差异 diff_num bp 的barcode
                continue
            # orient in {'split', 'records', 'index', 'columns', 'values', 'table'}
            group_data = ddf[[DBUtil.barcode, DBUtil.count]].to_dict(orient='split')["data"]
            group_data = sorted(group_data, key=lambda x: int(x[1]), reverse=True)
            # 将umi count最大的暂定为intended seq
            max_count = max([int(d[1]) for d in group_data])
            itd_list = [d for d in group_data if int(d[1]) == max_count]
            nbr_list = [d for d in group_data if int(d[1]) != max_count]
            # if len(nbr_list)==0, the list of itd_list can not build itd-nbr pair
            for (itd, n_itd), (nbr, n_nbr) in itertools.product(itd_list, nbr_list):
                result.append([
                    itd, nbr, n_itd,
                    n_nbr, pos,
                    itd[base_idx], nbr[base_idx]
                ])
        df = pd.DataFrame(data=result, columns=[
            DBUtil.intended_barcode,
            DBUtil.neighbor_barcode,
            DBUtil.intended_size,
            DBUtil.neighbor_size,
            DBUtil.position,
            DBUtil.intended_base,
            DBUtil.neighbor_base
        ])

        DBUtil.ItdNbr.df2db(self.db_util, df)

    def _find_all_pairs(self, pos):
        if pos != self.barcode_len:
            # 找缺失对,缺失位置不是最后一位都与最后一位缺失的列进行分组
            self._find_synthesis_error_pair(pos)
        self._find_substitution_error_pair(pos)
        return

    def search_to_db(self):
        logger.debug(f"barcode has {self.barcode_len}bp")
        barcode_site = set([i for i in range(1, self.barcode_len + 1)])
        pool = Pool(self.thread)
        # todo diff num must be 1
        tasks = list(itertools.combinations(barcode_site, self.barcode_len - self.diff_num))
        pbar = tqdm(total=len(tasks))
        pbar.set_description('step1: find all intended-neighbor pair in barcode')

        for cols in tasks:
            # pos start with 1
            pos = list(barcode_site - set(cols))[0]
            pool.apply_async(
                self._find_all_pairs, (pos,),
                callback=partial(multiprocess_callback, pbar),
                error_callback=self.write_error)
        pool.close()
        pool.join()
        pbar.close()
