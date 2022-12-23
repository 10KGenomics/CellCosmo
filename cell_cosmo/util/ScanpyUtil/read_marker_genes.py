#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : read_marker_genes.py
@Time       : 2022/07/05
@Version    : 1.0
@Desc       : None
"""
import os
import pandas as pd
from cell_cosmo.util import runtime, PathUtil
from .constant import (DiffExpSubDir, DiffExpFile,
                       NAME_DICT, FEATURE_ID, FEATURE_NAME,
                       FEATURE, L2FC_FN, ADJ_P_VAL_FN, TopDirFmt)

N = 50


def trans_col_name(col_name: str):
    # Feature
    # 'Cluster 15 Log2 fold change', 'Cluster 15 Adjusted p value
    if col_name.startswith(FEATURE):
        return [FEATURE, col_name.replace(FEATURE, "").strip()]
    elif col_name.endswith(L2FC_FN):
        return [col_name.replace(L2FC_FN, "").strip(), "L2FC"]
    elif col_name.endswith(ADJ_P_VAL_FN):
        return [col_name.replace(ADJ_P_VAL_FN, "").strip(), "p-value"]
    else:
        raise Exception(f"not support this col name trans,{col_name}")


def _get_display_data(df):
    # 过滤基因,只留 mean count >1.0 ,l2fc top 50 的基因
    n_cluster = len([c for c in df.columns.tolist() if c.endswith("Mean Counts")])
    cols = [FEATURE_ID, FEATURE_NAME]
    dfs = set()
    p_value_cols = []
    for g in range(1, n_cluster + 1):
        mean_counts_col = f"Cluster {g} Mean Counts"
        l2fc_col = f"Cluster {g} {L2FC_FN}"
        pVal_col = f"Cluster {g} {ADJ_P_VAL_FN}"
        cols.extend([l2fc_col, pVal_col])
        ndf = df[[FEATURE_NAME, mean_counts_col, l2fc_col]]
        ndf = ndf[ndf[mean_counts_col] > 1.0]
        ndf = ndf.sort_values(by=l2fc_col, axis=0, ascending=False)

        p_value_cols.append(pVal_col)
        dfs.update(ndf.head(N)[FEATURE_NAME].tolist())
    sub_df = df[df[FEATURE_NAME].isin(list(dfs))][cols]
    # 改变 p-value 的输出格式
    for pv_col in p_value_cols:
        # d = d.round({'L2FC': 3, 'p-value': 3})
        # d['L2FC'] = d['L2FC'].apply(lambda x: format(x, '.2%s'))
        sub_df[pv_col] = sub_df[pv_col].apply(lambda x: f"%.0e" % x)

    cur_dict = sub_df.to_dict(orient="split")
    del cur_dict['index']
    cur_dict["columns"] = [trans_col_name(col_name) for col_name in cur_dict["columns"]]
    return cur_dict


def read_marker_genes(sample, outdir):
    table_data = {"umap": {}, "tsne": {}}
    workdir = os.path.join(outdir, TopDirFmt % sample)
    with PathUtil.chdir(workdir):
        for d in os.listdir(DiffExpSubDir):
            # 过滤无关的目录
            if d[:4] in table_data or d in NAME_DICT:
                df = pd.read_csv(f"{DiffExpSubDir}/{d}/{DiffExpFile}")
                cur_dict = _get_display_data(df)
                if d[:4] in table_data:
                    # todo 需要更优雅的设计实现,最好不要使用名称进行组合
                    table_data[d[:4]][NAME_DICT[d[4:]]] = cur_dict
                else:
                    table_data['umap'][NAME_DICT[d]] = cur_dict
                    table_data['tsne'][NAME_DICT[d]] = cur_dict
    return table_data
