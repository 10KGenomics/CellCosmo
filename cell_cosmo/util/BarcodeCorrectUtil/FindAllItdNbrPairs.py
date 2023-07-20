#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : find.py
@Time       : 2023/07/19
@Version    : 1.0
@Desc       : None
"""
import glob
import logging
import itertools
import os

import pandas as pd
from tqdm import tqdm
from functools import partial
from multiprocessing import Pool
from cell_cosmo.util.BarcodeCorrectUtil import ConstNS

# import Levenshtein

logger = logging.getLogger(__name__)

miss_pos = "miss%s"


def multiprocess_callback(pbar: tqdm, res):
    pbar.update()


class FindAllItdNbrPairs:
    # modify: del params: diff_num=1, barcode_diff_limit=0.1,
    def __init__(self, outdir, sample, df: pd.DataFrame, barcode_len: int, thread=4):
        self.outdir = outdir
        self.sample = sample
        self.temp_file_output_fmt = f"{outdir}/.{sample}_%s_%s_cache"  # {[synthesis,substitution]}_{position}

        # self.db_util = db_util
        self.barcode_len = barcode_len  # TODO 超参数,1:barcode的长度都需要保持一致;2: 5<=len<=30
        self.df = df

        # TODO 差异的碱基数,当前算法中由于需要根据单个异变的位置,寻找itd-nbr对
        # todo 由于目前算法不支持2bp差异的itd-nbr对的搜索,所以暂时屏蔽相关参数的用户控制
        self.diff_num = 1
        self.thread = thread

    def write_error(self, e):
        # logger.exception(e)
        logger.error(f"[barcode correct error]:{str(e)}")
        # DBUtil.ErrorLog.write_log(self.db_util, msg="[barcode correct error]:" + str(e))

    def _find_synthesis_error_pair(self, pos):
        """寻找合成错误"""
        # 将列名提取出来，避免后续修改出错

        miss_fmt = miss_pos
        miss_col = miss_fmt % str(pos)
        miss_col_last = miss_fmt % self.barcode_len

        barcode, count = ConstNS.barcode, ConstNS.count
        intended_barcode = ConstNS.intended_barcode
        intended_size = ConstNS.intended_size
        intended_base = ConstNS.intended_base
        neighbor_barcode = ConstNS.neighbor_barcode
        neighbor_size = ConstNS.neighbor_size
        neighbor_base = ConstNS.neighbor_base
        position = ConstNS.position
        df1 = self.df[[barcode, count, miss_col, int(pos)]].rename({
            barcode: intended_barcode,
            count: intended_size,
            int(pos): intended_base
        }, axis=1).set_index(miss_col)
        df2 = self.df[[barcode, count, miss_col_last]].rename({
            barcode: neighbor_barcode,
            count: neighbor_size,
        }, axis=1).set_index(miss_col_last)

        df3 = df1.join(df2, how='inner')
        df4 = df3[df3[intended_barcode] != df3[neighbor_barcode]].copy()
        df4[position] = str(pos)
        df4[neighbor_base] = "-"
        df4 = df4[[
            intended_barcode,
            neighbor_barcode,
            intended_size,
            neighbor_size,
            position,
            intended_base,
            neighbor_base
        ]]
        # write to csv

        filepath = self.temp_file_output_fmt % ("synthesis", str(pos))
        df4.to_csv(filepath, index=False, sep="\t")
        # DBUtil.ItdNbr.df2db(self.db_util, df4)

    def _find_substitution_error_pair(self, pos):
        # 将列名提取出来，避免后续修改出错
        miss_fmt = miss_pos
        miss_col = miss_fmt % str(pos)
        barcode = ConstNS.barcode
        count = ConstNS.count
        intended_barcode = ConstNS.intended_barcode
        neighbor_barcode = ConstNS.neighbor_barcode
        intended_size = ConstNS.intended_size
        neighbor_size = ConstNS.neighbor_size
        intended_base = ConstNS.intended_base
        neighbor_base = ConstNS.neighbor_base
        position = ConstNS.position

        miss_col_last = miss_fmt % self.barcode_len
        base_idx = int(pos) - 1
        result = []
        for key, ddf in self.df.groupby(miss_col):
            if ddf.shape[0] <= 1:
                # 没有与其差异 diff_num bp 的barcode
                continue
            # orient in {'split', 'records', 'index', 'columns', 'values', 'table'}
            group_data = ddf[[barcode, count]].to_dict(orient='split')["data"]
            group_data = sorted(group_data, key=lambda x: int(x[1]), reverse=True)
            # 将umi count最大的暂定为intended seq
            max_count = max([int(d[1]) for d in group_data])
            itd_list = [d for d in group_data if int(d[1]) == max_count]
            nbr_list = [d for d in group_data if int(d[1]) != max_count]
            # if len(nbr_list)==0, the list of itd_list can not build itd-nbr pair
            for (itd, n_itd), (nbr, n_nbr) in itertools.product(itd_list, nbr_list):
                result.append([
                    itd, nbr, n_itd,
                    n_nbr, str(pos),
                    itd[base_idx], nbr[base_idx]
                ])
        df = pd.DataFrame(data=result, columns=[
            intended_barcode,
            neighbor_barcode,
            intended_size,
            neighbor_size,
            position,
            intended_base,
            neighbor_base
        ])
        filepath = self.temp_file_output_fmt % ("substitution", str(pos))
        df.to_csv(filepath, index=False, sep="\t")
        # DBUtil.ItdNbr.df2db(self.db_util, df)

    def _find_all_pairs(self, pos: int):
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
        # 合并数据
        fmt = self.temp_file_output_fmt.replace("%s", "*")

        frames = []
        for f in glob.glob(fmt):
            a_df = pd.read_csv(f, sep="\t")
            frames.append(a_df)
            os.remove(f)  # 删除临时文件
        cct_df = pd.concat(frames)
        pbar.close()
        return cct_df


if __name__ == '__main__':
    _out = "/mnt/data/cellcosmo/xxx"
    _df = pd.read_csv(f"{_out}/s1_extract_sub_table.csv", sep='\t')
    fa_inp = FindAllItdNbrPairs(_out, "s1", _df, 24, thread=4)
    fa_inp.search_to_db()
