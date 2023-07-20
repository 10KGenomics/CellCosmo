#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : get_correct_dict.py
@Time       : 2023/07/19
@Version    : 1.0
@Desc       : None
"""
import abc
import logging
import pandas as pd
from tqdm import tqdm
from cell_cosmo.util.BarcodeCorrectUtil import ConstNS
from cell_cosmo.util.BarcodeCorrectUtil.BaseOut import BaseOut

# from cell_cosmo.util.BarcodeCorrectUtil import filter
# from cell_cosmo.util.BarcodeCorrectUtil.db_util import DBUtil

logger = logging.getLogger(__name__)

intended_base = ConstNS.intended_base
neighbor_base = ConstNS.neighbor_base
intended_size = ConstNS.intended_size
neighbor_size = ConstNS.neighbor_size
intended_barcode = ConstNS.intended_barcode
neighbor_barcode = ConstNS.neighbor_barcode
filter_reason = ConstNS.filter_reason
position = ConstNS.position

percent = ConstNS.percent
neighbors = ConstNS.neighbors
n_neighbors = ConstNS.n_neighbors
neighbor_size_total = ConstNS.neighbor_size_total


class FilterBase(metaclass=abc.ABCMeta):
    rank = "rank"

    def __init__(self, pbar: tqdm):
        self.pbar = pbar
        self.processing = ""
        self.reason_str = ""

    @abc.abstractmethod
    def _run(self, df):
        pass

    def run(self, df):
        self.pbar.set_postfix(processing=self.processing)
        res = self._run(df)
        self.pbar.update(1)
        return res  # None or df


class FilterItdSizeLeNbrSize(FilterBase):
    def __init__(self, pbar: tqdm):
        super(FilterItdSizeLeNbrSize, self).__init__(pbar)
        self.processing = "filter itd count < nbr count"
        self.reason_str = 'filter_by_step1: itd size <= nbr size'

    def _run(self, df):
        data_slice = df[intended_size] <= df[neighbor_size]
        df.loc[data_slice, filter_reason] = self.reason_str


class FilterDuplicateSameBaseInItd(FilterBase):
    def __init__(self, pbar: tqdm):
        super(FilterDuplicateSameBaseInItd, self).__init__(pbar)
        self.processing = "filter same base in itd"
        self.reason_str = 'filter_by_step2: duplicate by same base in itd'

    def _run(self, df):
        df[self.rank] = 10000  # 随便给一个较大的初始值
        data_slice = (df[neighbor_base] == '-') & (df[filter_reason] == "")
        df.loc[data_slice, self.rank] = df[data_slice].groupby(
            [intended_barcode, neighbor_barcode]
        )[position].rank(method='first')
        update_slice = (df[filter_reason] == "") & (df[self.rank] != 1)
        df.loc[update_slice, filter_reason] = self.reason_str


class FilterDuplicateBySubstitutionWithIdl(FilterBase):
    def __init__(self, pbar: tqdm):
        super(FilterDuplicateBySubstitutionWithIdl, self).__init__(pbar)
        self.processing = "filter substitution with idl"
        self.reason_str = 'filter_by_step3: duplicate by substitution with idl'

    def _run(self, df):
        df[self.rank] = 10000  # 重新使用rank列按新的条件排序
        data_slice = df[filter_reason] == ""
        df.loc[data_slice, self.rank] = df[data_slice].groupby(
            [neighbor_barcode, intended_barcode]
        )[neighbor_base].rank(method='first')
        update_slice = (df[filter_reason] == "") & (df[self.rank] != 1)
        df.loc[update_slice, filter_reason] = self.reason_str


# filter itd with mult nbr
class FilterItdWithMultNbr(FilterBase):
    def __init__(self, pbar: tqdm):
        super(FilterItdWithMultNbr, self).__init__(pbar)
        self.processing = "filter itd with mult nbr"
        self.reason_str = 'filter_by_step4: itd with mult nbr'

    def _run(self, df):
        df[self.rank] = 10000  # 重新使用rank列按新的条件排序
        data_slice = df[filter_reason] == ""
        df.loc[data_slice, self.rank] = df[data_slice].groupby(
            [neighbor_barcode]
        )[intended_size].rank(method='max', ascending=True)
        update_slice = (df[filter_reason] == "") & (df[self.rank] != 1)
        df.loc[update_slice, filter_reason] = self.reason_str


class FilterItdInNbr(FilterBase):
    def __init__(self, pbar: tqdm):
        super(FilterItdInNbr, self).__init__(pbar)
        self.processing = "filter itd in nbr"
        self.reason_str = 'filter_by_step5: itd in nbr'

    def _run(self, df):
        update_slice = (df[filter_reason] == "") & (
            df[intended_barcode].isin(df[neighbor_barcode]))
        df.loc[update_slice, filter_reason] = self.reason_str


class AggItdAllNbr(FilterBase):
    def __init__(self, pbar: tqdm):
        super(AggItdAllNbr, self).__init__(pbar)
        self.processing = "aggregate itd all nbr"
        self.reason_str = 'filter_by_step5: itd in nbr'

    def _run(self, df):
        # 将满足条件的itd-nbr pair按itd分组聚合
        data_slice = df[filter_reason] == ""
        sub_df = df[data_slice]
        data = []
        for (itd_b, itd_s), ddf in sub_df.groupby([intended_barcode, intended_size]):
            meta = []
            nbr_size = 0
            nbr_n = 0
            for d in ddf[[intended_barcode,  # 0
                          neighbor_size,  # 1
                          intended_base,  # 2
                          neighbor_base,  # 3
                          position  # 4
                          ]].values.tolist():
                nbr_size += d[1]
                nbr_n += 1
                meta.append(f'{d[0]}:{d[1]}@POS[{d[4]}]{d[2]}->{d[3]}')
            nbr = ";".join(meta)
            data.append([itd_b, itd_s, nbr, nbr_n, nbr_size, round(nbr_size / itd_s, 2)])

        grouped_df = pd.DataFrame(data=data,
                                  columns=[
                                      intended_barcode,
                                      intended_size,
                                      neighbors,
                                      n_neighbors,
                                      neighbor_size_total,
                                      percent,
                                  ])
        return grouped_df

    def get_filter_slice(self, grouped_df: pd.DataFrame, limit):
        le_percent_slice = grouped_df.loc[grouped_df[percent] > limit, intended_barcode]
        return le_percent_slice


class FilterBarcodeNbrUMIsLtLimit(FilterBase):
    def __init__(self, pbar: tqdm, data_slice, limit=0.01):
        super(FilterBarcodeNbrUMIsLtLimit, self).__init__(pbar)
        self.processing = f"filter barcode nbr umis > {limit}"
        self.reason_str = f'filter_by_step6: percent of nbr umis < {limit}'
        self.filter_limit = limit
        self.data_slice = data_slice

    def _run(self, df):
        # less than percent，这些数据是需要校正的
        data_slice = (df[filter_reason] == "") & (df[intended_barcode].isin(self.data_slice))
        df.loc[data_slice, filter_reason] = self.reason_str


class NeedCorrectBarcode(FilterBase):

    def __init__(self, pbar: tqdm, limit=0.01):
        super(NeedCorrectBarcode, self).__init__(pbar)
        self.processing = "obtain barcode correct dict"
        self.reason_str = "Passed:NeedCorrectBarcode"
        self.filter_limit = limit

    def _run(self, df):
        data_slice = df[filter_reason] == ""
        itd_nbr_df = df[data_slice].copy()
        df.loc[data_slice, filter_reason] = "Passed:NeedCorrectBarcode"

        # 去除一些不需要的列
        remove_cols = list(set(itd_nbr_df.columns.tolist()) &
                           {"id", filter_reason, self.rank})

        itd_nbr_df.drop(remove_cols, axis=1, inplace=True)
        return itd_nbr_df


def get_correct_dict(df, out_util: BaseOut, filter_limit=0.01) -> dict:
    """

    :param df: 输入数据
    :param out_util: 输出帮助类
    :param filter_limit: nbr/itd > filter_limit will be filtered
    :return:
    """

    rank = "rank"
    df[filter_reason] = ""
    df[rank] = 100000

    pbar = tqdm(total=8)
    pbar.set_description('step2: filter invalid itd-nbr pairs in barcode')

    FilterItdSizeLeNbrSize(pbar).run(df)
    FilterDuplicateSameBaseInItd(pbar).run(df)
    FilterDuplicateBySubstitutionWithIdl(pbar).run(df)
    FilterItdWithMultNbr(pbar).run(df)
    FilterItdInNbr(pbar).run(df)

    agg_itd_all_nbr = AggItdAllNbr(pbar)
    grouped_df = agg_itd_all_nbr.run(df)
    # todo grouped_df 记录到数据表
    # DBUtil.ItdNbrGroup.df2db(db_util, grouped_df)
    le_percent_slice = agg_itd_all_nbr.get_filter_slice(grouped_df, filter_limit)

    FilterBarcodeNbrUMIsLtLimit(pbar, le_percent_slice).run(df)
    need_correct_barcode = NeedCorrectBarcode(pbar).run(df)

    # 去除一些不需要的列
    remove_cols = list(set(df.columns.tolist()) & {"id", rank})
    df.drop(remove_cols, axis=1, inplace=True)
    # df 添加原因后记录到文件
    out_util.to_csv("barcode_itd_nbr_pairs.raw", df)
    out_util.to_csv("barcode_itd_nbr_pairs", need_correct_barcode)
    correct_dict = dict(need_correct_barcode[[neighbor_barcode, intended_barcode]].values.tolist())

    pbar.close()
    return correct_dict


if __name__ == '__main__':
    _out = "/mnt/data/cellcosmo/xxx"
    _df = pd.read_csv(f"{_out}/.s1_barcode_itd_nbr.tsv", sep='\t')

    get_correct_dict(_df, BaseOut(_out, "s1"))
