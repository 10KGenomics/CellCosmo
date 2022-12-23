#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : db_util.py
@Time       : 2022/08/04
@Version    : 1.0
@Desc       : None
"""
import abc
import logging
import os.path
import pandas as pd
from sqlalchemy import Column, String, Integer, Float, Text, create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

logger = logging.getLogger(__name__)

Base = declarative_base()


class DBUtilBase(metaclass=abc.ABCMeta):
    def __init__(self, outdir, sample):
        self.db_file = f"{outdir}/.temp.db"
        self.out_prefix = f"{outdir}/{sample}"
        self._name2model = {}
        self._init_tables()

    def _get_engine(self):
        # 多进程中engine无法 pickle
        return create_engine(f'sqlite:///{self.db_file}')

    def _get_session(self):
        engine = self._get_engine()
        return sessionmaker(bind=engine)

    def _init_tables(self):
        engine = self._get_engine()
        Base.metadata.create_all(engine)

    def model_name_register(self, name, clazz):
        self._name2model[name] = clazz

    def df2db(self, df: pd.DataFrame, table_name: str):
        if df.shape[0] == 0:
            return
        engine = self._get_engine()
        df.to_sql(name=table_name, con=engine, if_exists='append', index=False)

    def get_table_columns(self, table_name):
        model_clazz = self._name2model.get(table_name, None)
        if model_clazz is None:
            raise Exception(f"{table_name} not exists!")
        cols = [col.name for col in model_clazz.__table__.columns if col.name != "id"]
        return cols

    def db2df(self, table_name):
        cols = self.get_table_columns(table_name)
        return self.query(sql=f"SELECT * FROM {table_name}")

    def query(self, sql):
        # 将数字类型的字段转换成数字
        engine = self._get_engine()
        df = pd.read_sql_query(sql, engine)
        return df

    def to_csv(self, filename, df: pd.DataFrame):
        df.to_csv(
            f"{self.out_prefix}_{filename}", sep='\t', index=False
        )


class DBUtil(DBUtilBase):
    class BaseModel(Base):
        __abstract__ = True
        __tablename__ = None
        id = Column(Integer, autoincrement=True, primary_key=True)

        @classmethod
        def db2df(cls, db: DBUtilBase):
            return db.db2df(cls.__tablename__)

        @classmethod
        def df2db(cls, db: DBUtilBase, df: pd.DataFrame):
            db.df2db(df, cls.__tablename__)

        @classmethod
        def query(cls, db: DBUtilBase, sql):
            return db.query(sql=sql)

    class FilterStat(BaseModel):
        __tablename__ = 'filter_log'

        id = Column(Integer, autoincrement=True, primary_key=True)
        desc = Column(Text, nullable=False)

        @classmethod
        def write_log(cls, db: DBUtilBase, msg: str):
            df = pd.DataFrame(data=[msg], columns=["desc"])
            cls.df2db(db, df)

    class ErrorLog(BaseModel):
        __tablename__ = 'error_log'

        id = Column(Integer, autoincrement=True, primary_key=True)
        desc = Column(Text, nullable=False)

        @classmethod
        def write_log(cls, db: DBUtilBase, msg: str):
            df = pd.DataFrame(data=[msg], columns=["desc"])
            cls.df2db(db, df)

    class CountDetailAll(BaseModel):
        __tablename__ = 'full_data'  # 'count_detail_all'

        id = Column(Integer, autoincrement=True, primary_key=True)
        barcode = Column(Text, nullable=False)
        gene_id = Column(Text, nullable=False)
        umi = Column(Text, nullable=False)
        count = Column(Integer, nullable=False)
        miss_pos = "miss%d"

        @classmethod
        def check_df(cls, db: DBUtilBase, df: pd.DataFrame) -> (int, int):
            # 对于full_data数据的入库对列名进行检查，以便后续操作都在预期的schema下进行
            columns_check_set = set(db.get_table_columns(cls.__tablename__))
            for col in df.columns.tolist():
                assert col in columns_check_set, f"the column {col} not in file header!"
            # 对 barcode 和 umi length 进行检查
            res = db.query(sql=f"SELECT distinct length(barcode) "
                               f"as barcode_len FROM {cls.__tablename__};")
            assert res.shape[0] == 1, f"the barcode has inconsistent length,please check it!!"
            barcode_len = int(res.iloc[0, 0])  # todo 5<=len<=30,是否需要验证一下？

            res = db.query(sql=f"SELECT distinct length(umi) "
                               f"as umi_len FROM {cls.__tablename__};")
            assert res.shape[0] == 1, f"the umi has inconsistent length,please check it!!"
            umi_len = int(res.iloc[0, 0])
            return barcode_len, umi_len

        @classmethod
        def get_barcode_matrix(cls, db: DBUtilBase, barcode_len, n_umi_filter=20) -> pd.DataFrame:
            barcode, count = "barcode", "count"  # 将列名提取出来，避免后续修改出错

            stat_df = db.query(sql=f"SELECT {barcode},SUM(count) as {count} "
                                   f"FROM {cls.__tablename__} GROUP BY {barcode};")
            # 按需求逻辑， barcode umi count <= 20,这些数据将被去除
            le_df = stat_df[stat_df[count] <= n_umi_filter]
            df = stat_df[stat_df[count] > n_umi_filter].copy()
            db.to_csv(f"barcode_umi_less_than_{n_umi_filter}_filter_{le_df.shape[0]}_data.tsv", le_df)
            logger.info(f"there is {le_df.shape[0]} barcode umi less than {n_umi_filter}")

            # 构建缺失列，以便后续分析
            remove_cols = [0, barcode_len + 1]
            # 过滤0
            ndf = pd.concat([df, df[barcode].str.split('', expand=True)], axis=1)
            ndf.drop(remove_cols, axis=1, inplace=True)
            # 将缺失序列提前构建
            for i in range(1, barcode_len + 1):
                cols = [c for c in range(1, barcode_len + 1) if c != i]
                ndf[cls.miss_pos % i] = ndf[list(cols)].apply("".join, axis=1)
            return ndf

    class CountDetailAll2(BaseModel):
        __tablename__ = 'full_data2'  # 'count_detail_all'

        id = Column(Integer, autoincrement=True, primary_key=True)
        barcode = Column(Text, nullable=False)
        gene_id = Column(Text, nullable=False)
        umi = Column(Text, nullable=False)
        count = Column(Integer, nullable=False)

    class ItdNbr(BaseModel):
        __tablename__ = "barcode_itd_nbr"  # 'intended_neighbor'

        id = Column(Integer, autoincrement=True, primary_key=True)
        intended_barcode = Column(Text, nullable=False)
        neighbor_barcode = Column(Text, nullable=False)
        intended_size = Column(Integer, nullable=False)
        neighbor_size = Column(Integer, nullable=False)
        position = Column(Integer, nullable=False)
        intended_base = Column(Text, nullable=False)
        neighbor_base = Column(Text, nullable=False)

    class ItdNbrReason(BaseModel):
        __tablename__ = "barcode_itd_nbr_reason"  # 'intended_neighbor'

        id = Column(Integer, autoincrement=True, primary_key=True)
        intended_barcode = Column(Text, nullable=False)
        neighbor_barcode = Column(Text, nullable=False)
        intended_size = Column(Integer, nullable=False)
        neighbor_size = Column(Integer, nullable=False)
        position = Column(Integer, nullable=False)
        intended_base = Column(Text, nullable=False)
        neighbor_base = Column(Text, nullable=False)
        filter_reason = Column(Text, default="")

        @classmethod
        def update_filter_data(cls, db: DBUtilBase, save_df: pd.DataFrame, filter_df: pd.DataFrame, reason: str):
            """df1: 需要保留的df； df2: 按reason条件过滤的数据"""

            n1, uniq_n1 = save_df.shape[0], save_df[DBUtil.intended_barcode].unique().size
            n2, uniq_n2 = filter_df.shape[0], filter_df[DBUtil.intended_barcode].unique().size

            filter_log = f'当前共用{n1 + n2}条数据,{uniq_n1 + uniq_n2}个唯一barcode出现。' \
                         f'通过`{reason}`过滤后：' \
                         f'过滤{n2}条数据,包含{uniq_n2}个barcode；' \
                         f'剩余{n1}条数据,包含{uniq_n1}个barcode。'

            DBUtil.FilterStat.write_log(db, msg=filter_log)
            filter_df[DBUtil.filter_reason] = reason
            cls.df2db(db, filter_df)  # tbl_itd_nbr_filter

        @classmethod
        def update_others(cls, db: DBUtilBase, save_df: pd.DataFrame, reason="Passed:NeedCorrectBarcode"):
            save_df[DBUtil.filter_reason] = reason
            cls.df2db(db, save_df)

    class ItdNbrGroup(BaseModel):
        __tablename__ = 'intended_neighbor_group'

        id = Column(Integer, autoincrement=True, primary_key=True)
        intended_barcode = Column(Text, nullable=False)
        intended_size = Column(Integer, nullable=False)
        neighbors = Column(Text, nullable=False)
        n_neighbors = Column(Integer, nullable=False)
        neighbor_size_total = Column(Integer, nullable=False)
        percent = Column(Float, nullable=False)

    # 处理过程中使用到的一些字段名
    # count_detail 类型结果中使用的表名
    barcode = CountDetailAll.barcode.name
    gene_id = CountDetailAll.gene_id.name
    umi = CountDetailAll.umi.name
    count = CountDetailAll.count.name
    # itd_nbr 类型结果中使用的表名
    intended_barcode = ItdNbr.intended_barcode.name
    intended_size = ItdNbr.intended_size.name
    intended_base = ItdNbr.intended_base.name
    neighbor_barcode = ItdNbr.neighbor_barcode.name
    neighbor_size = ItdNbr.neighbor_size.name
    neighbor_base = ItdNbr.neighbor_base.name
    position = ItdNbr.position.name
    # 额外增加的字段
    neighbors = ItdNbrGroup.neighbors.name
    neighbor_size_total = ItdNbrGroup.neighbor_size_total.name
    n_neighbors = ItdNbrGroup.n_neighbors.name
    percent = ItdNbrGroup.percent.name
    filter_reason = ItdNbrReason.filter_reason.name
    # 其他表中使用到的列名
    miss_pos = "miss%d"
    key_col = "temp_key"
    error_desc = ErrorLog.desc.name

    COUNT_DETAIL_COLS = [barcode, gene_id, umi, count]

    def __init__(self, outdir, sample):
        super(DBUtil, self).__init__(outdir, sample)
        self.model_name_register(DBUtil.FilterStat.__tablename__, DBUtil.FilterStat)
        self.model_name_register(DBUtil.ErrorLog.__tablename__, DBUtil.ErrorLog)
        self.model_name_register(DBUtil.CountDetailAll.__tablename__, DBUtil.CountDetailAll)
        self.model_name_register(DBUtil.CountDetailAll2.__tablename__, DBUtil.CountDetailAll2)
        self.model_name_register(DBUtil.ItdNbr.__tablename__, DBUtil.ItdNbr)
        self.model_name_register(DBUtil.ItdNbrReason.__tablename__, DBUtil.ItdNbrReason)
        self.model_name_register(DBUtil.ItdNbrGroup.__tablename__, DBUtil.ItdNbrGroup)
        self.clean_db()

    def clean_db(self):
        if os.path.exists(self.db_file):
            try:
                os.remove(self.db_file)
            except Exception as e:
                logger.error(f"db 文件删除失败，请手动删除。{str(e)}")
