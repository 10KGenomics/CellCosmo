#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Author     : ice-melt@outlook.com
@File       : base_report_runner.py
@Time       : 2022/06/07
@Version    : 1.0
@Desc       : None
"""
import abc
import sys
import io
import json
import numbers
import os
import subprocess
from jinja2 import Environment, FileSystemLoader, select_autoescape
from cell_cosmo.output_runner.base_runner import BaseRunner
from cell_cosmo.util import runtime
import logging

logger = logging.getLogger(__name__)


def cap_str_except_preposition(my_string):
    prepositions = {"and", "or", "the", "a", "of", "in", "per", "after", 'with'}
    lowercase_words = my_string.split(" ")

    final_words = [word if word in prepositions else word[0].upper() + word[1:] for word in lowercase_words]
    final_words = " ".join(final_words)
    return final_words


class TemplateReader:
    def __init__(self, output):
        # self.report_html = f"{self.outdir}/../{self.sample}_report.html"
        self.report_html = output
        # jinja env
        self._env = Environment(
            loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "templates")),
            autoescape=select_autoescape(['html', 'xml'])
        )

    def render_html(self, data: dict, assay: str, sample: str):
        template = self._env.get_template(f'html/{assay}/base.html')
        with io.open(self.report_html, 'w', encoding='utf8') as f:
            from cell_cosmo.__main__ import __VERSION__
            report_data = data.copy()
            # 页面需要展示版本号和样本号
            report_data["version"] = __VERSION__
            report_data["sample"] = sample
            sequence_summary = "sequence_summary"
            count_summary = "count_summary"
            analysis_summary = "analysis_summary"

            sequence_data = report_data.get(sequence_summary, None)
            count_data = report_data.get(count_summary, None)
            analysis_data = report_data.get(analysis_summary, None)

            # 去除Sequencing模块中以下指标的展示
            # Total Base Pairs Processed
            # Base Pairs Quality-Trimmed
            # Base Pairs Written(filtered)
            if sequence_data is not None:
                new_data = []
                cur_data = sequence_data.get('metric_list', [])
                for m in cur_data:
                    if m["name"] not in {
                        "Total Base Pairs Processed",
                        "Base Pairs Quality-Trimmed",
                        "Base Pairs Written(filtered)",
                    }:
                        new_data.append(m)
                report_data[sequence_summary]['metric_list'] = new_data
            # 将 Sequencing Saturation图和Median Genes per Cell图
            # 放入 Analysis 中
            if count_data is not None and analysis_data is not None:
                count_help, analysis_help = [], []
                line_saturation = count_data.get("line_saturation", None)
                line_median = count_data.get("line_median", None)
                if line_saturation is not None:
                    del count_data['line_saturation']
                    analysis_data['line_saturation'] = line_saturation
                if line_median is not None:
                    del count_data['line_median']
                    analysis_data['line_median'] = line_median
                for m in count_data['help_content']:
                    if m['name'] in {
                        "Sequencing Saturation Plot",
                        "Median Genes per Cell Plot",
                    }:
                        analysis_help.append(m)
                    else:
                        count_help.append(m)
                analysis_help.extend(analysis_data.get('help_content', []))
                count_data["help_content"] = count_help
                analysis_data["help_content"] = analysis_help

            html = template.render(report_data)
            f.write(html)


class BaseReportRunner(BaseRunner):
    _STEP_NAME = None
    _DISPLAY_TITLE = None

    def __init__(self, outdir, sample, thread=4, debug=False, **kwargs):
        super(BaseReportRunner, self).__init__()
        if self._STEP_NAME is None:
            raise Exception(f"{self.__class__.__name__} must set _STEP_NAME value")
        if self._DISPLAY_TITLE is None:
            # raise Exception(f"{self.__class__.__name__} must set _DISPLAY_TITLE value")
            self._DISPLAY_TITLE = self.__class__.__name__
        self.display_title = self._DISPLAY_TITLE
        self._display_title = self._DISPLAY_TITLE
        self._step_name = self._STEP_NAME

        # class_name = self.__class__.__name__
        # self._step_name = class_name[0].lower() + class_name[1:]

        # outdir, sample, subparser_assay,
        self.outdir = outdir  # kwargs.get("output_dir", "./")
        self.sample = sample  # kwargs.get("sample_name", "ttt")
        # TODO subparser_assay 参数的处理,暂时由command模块自动填充
        self.assay = kwargs.get("subparser_assay")
        self.thread = int(thread)
        self.debug = debug

        self.out_prefix = f'{self.outdir}/{self.sample}'

        self.__slots = ['data', 'metrics']
        self._step_summary_name = f'{self._step_name}_summary'

        self.__metric_list = []
        self.__help_content = []
        self._path_dict = {}
        for slot in self.__slots:
            self._path_dict[slot] = f'{self.outdir}/../.{slot}.json'

        self.__content_dict = {}
        for slot, path in self._path_dict.items():
            if not os.path.exists(path):
                self.__content_dict[slot] = {}
            else:
                with open(path) as f:
                    self.__content_dict[slot] = json.load(f)
            # clear step_summary
            self.__content_dict[slot][self._step_summary_name] = {}

        # jinja env
        # self.env = Environment(
        #     loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), "templates")),
        #     autoescape=select_autoescape(['html', 'xml'])
        # )
        self.env = TemplateReader(
            output=f"{self.outdir}/../{self.sample}_report.html"
        )

        # out file
        self.__stat_file = f'{self.outdir}/stat.txt'

    @abc.abstractmethod
    def collect_matrix(self):
        pass

    @property
    def assay_text(self):
        return 'Single-cell ' + self.assay

    def add_metric(self, name, value, total=None, help_info=None, display=None, show=True):
        """
        add metric to metric_list
        :param name:
        :param value:
        :param total:   int or float, used to calculate fraction
        :param help_info: str, help info for metric in html output_runner
        :param display: str, controls how to display the metric in HTML output_runner.
        :param show:    bool, whether to add to `.data.json` and `stat.txt`.
                        `.data.json` is used for HTML output_runner. `stat.txt` is used in house.
        :return:
        """

        name = cap_str_except_preposition(name)
        if help_info:
            help_info = help_info[0].upper() + help_info[1:]
            if help_info[-1] != '.':
                help_info += '.'
        if not display:
            if isinstance(value, numbers.Number):
                display = str(format(value, ','))
            else:
                display = value
        fraction = None
        if total:
            fraction = round(value / total * 100, 2)
            display += f'({fraction}%)'
        self.__metric_list.append(
            {
                "name": name,
                "value": value,
                "total": total,
                "fraction": fraction,
                "display": display,
                "help_info": help_info,
                "show": show,
            }
        )

    def _write_stat(self):
        with open(self.__stat_file, 'w') as writer:
            for metric in self.__metric_list:
                if metric['show']:
                    name = metric['name']
                    display = metric['display']

                    line = f'{name}: {display}'
                    writer.write(line + '\n')

    def _dump_content(self):
        """dump content to json file"""
        for slot, path in self._path_dict.items():
            if self.__content_dict[slot]:
                with open(path, 'w') as f:
                    json.dump(self.__content_dict[slot], f, indent=4)

    @runtime(f"{__name__}.render_html")
    def _render_html(self):
        self.env.render_html(
            self.__content_dict['data'],
            self.assay,
            self.sample
        )

    def _add_content_data(self):
        step_summary = {'display_title': self._display_title}
        metric_list = []
        for metric in self.__metric_list:
            if metric['show']:
                metric_list.append(metric)
        step_summary['metric_list'] = metric_list
        step_summary['help_content'] = self.__help_content
        self.__content_dict['data'][self._step_summary_name].update(step_summary)

    def _add_content_metric(self):
        metric_dict = dict()
        for metric in self.__metric_list:
            name = metric['name']
            value = metric['value']
            fraction = metric['fraction']
            metric_dict[name] = value
            if fraction:
                metric_dict[f'{name} Fraction'] = fraction

        self.__content_dict['metrics'][self._step_summary_name].update(metric_dict)

    def add_data(self, **kwargs):
        """
        add data(other than metrics) to self.content_dict['data']
        for example: add plots and tables
        """
        for key, value in kwargs.items():
            self.__content_dict['data'][self._step_summary_name][key] = value

    def add_help_content(self, name, content):
        """add help info before metrics' help_info"""
        if content and content[-1] != ".":
            content = content + "."
        self.__help_content.append(
            {
                'name': name,
                'content': content
            }
        )

    # @utils.add_log # TODO 是否打印用時
    def get_slot_key(self, slot, step_name, key):
        """read slot from json file"""
        try:
            # print("sssssssssssssssssssss", step_name)
            return self.__content_dict[slot][step_name + '_summary'][key]
        except KeyError:
            # self.get_slot_key.logger.warning(f'{key} not found in {step_name}_summary.{slot}')
            logger.warning(f'{key} not found in {step_name}_summary.{slot}')
            raise

    def get_table_dict(self, title, table_id, df_table):
        """
        table_dict {title: '', table_id: '', df_table: pd.DataFrame}
        """
        table_dict = {'title': title, 'table': df_table.to_html(
            escape=False,
            index=False,
            table_id=table_id,
            justify="center"), 'id': table_id}
        return table_dict

    def merge_step(self, key1, key2, new_key):
        key1, key2 = [f"{key1}_summary", f"{key2}_summary"]
        data = self.__content_dict['data']
        metrics = self.__content_dict['metrics']
        if all([n in data for n in [key1, key2]]):
            data[f'{new_key}_summary'] = {
                'display_title': 'Sequencing',
                'metric_list': data[key1]['metric_list'] + data[key2]['metric_list']
            }
            del data[key1]
            del data[key2]
        if all([n in metrics for n in [key1, key2]]):
            metrics[f'{new_key}_summary'] = dict(metrics[key1], **metrics[key2])
            del metrics[key1]
            del metrics[key2]
        # 为了保证合并后数据的顺序，这里重新对字典进行排序
        sort_list = [
            'sample_summary',
            'sequence_summary',
            'star_summary',
            'featureCounts_summary',
            'count_summary',
            'analysis_summary',
            'summary_summary',
        ]
        self.__content_dict['data'] = dict(sorted(data.items(), key=lambda x: sort_list.index(x[0])))
        self.__content_dict['metrics'] = dict(sorted(metrics.items(), key=lambda x: sort_list.index(x[0])))

    def process_starsolo_demultiplexing_step(self):
        """starsolo中转换barcode报告指标结果"""
        # key1, key2 = [f"{key1}_summary", f"{key2}_summary"]
        data = self.__content_dict['data']
        metrics = self.__content_dict['metrics']
        if "barcode_summary" in data:
            data['sequence_summary'] = {
                'display_title': 'Sequencing',
                'metric_list': data['barcode_summary'].get('metric_list', [])
            }
            del data['barcode_summary']
        if "barcode_summary" in metrics:
            metrics['sequence_summary'] = metrics['barcode_summary']
            del metrics['barcode_summary']

        # 为了保证合并后数据的顺序，这里重新对字典进行排序
        sort_list = [
            'sample_summary',
            'sequence_summary',
            'star_summary',
            'featureCounts_summary',
            'count_summary',
            'analysis_summary',
            'summary_summary',
        ]
        self.__content_dict['data'] = dict(sorted(data.items(), key=lambda x: sort_list.index(x[0])))
        self.__content_dict['metrics'] = dict(sorted(metrics.items(), key=lambda x: sort_list.index(x[0])))

    def set_summary_step(self):
        data = self.__content_dict['data']
        metrics = self.__content_dict['metrics']
        # extract info for summary display,step is count_summary
        if 'count_summary' in data:
            # 这里的名称需要和 Count 类中一致
            extract_list = ['Estimated Number of Cells',
                            'Median UMI Counts per Cell',
                            'Median Genes per Cell',
                            'Mean Reads per Cell']
            # Sample Name
            # Estimated Number of Cells
            # Mean Reads per Cell
            # Median Genes per Cell
            # Median UMI Counts per Cell

            metric_list = [metric for metric in data['count_summary']['metric_list']
                           if metric['name'] in extract_list]
            data['summary_summary'] = {
                'display_title': 'Summary',
                'metric_list': metric_list,
            }

            metrics['summary_summary'] = dict([
                (key, val) for key, val in metrics['count_summary'].items()
                if key in extract_list])

    def _mtx2tsv(self):
        # 新增需求,用户想可以打开文件复制统计结果,
        # 这里直接将.metrics.json的文件内容转成tsv输出到当前目录
        os.path.join(os.path.dirname(self.outdir))
        slot = self.__slots[1]  # ['data', 'metrics'],取metrics
        path = self._path_dict[slot]
        if not os.path.exists(path):
            return
        try:
            # {Sample}_summary.tsv
            outfile = os.path.join(self.outdir, f"{self.sample}_summary.tsv")
            with open(path) as f, open(outfile, "w") as o:
                metrics = json.load(f)
                data = []
                for k, v in metrics.items():
                    if k == "summary_summary":
                        # summary 下是一些重复的指标
                        continue
                    for a, b in v.items():
                        # 20220810,客户委托去除一些数据的展示
                        if str(a).startswith("Fraction of Cells Have Mito Gene Percent>"):
                            continue
                        if str(a).startswith("Read Fraction "):
                            continue
                        data.append(f"{a}\t{b}")
                o.write("\n".join(data))
        except Exception as e:
            logger.error(str(e))

    @runtime(f"{__name__}.clean_up")
    def _clean_up(self):
        self._add_content_data()
        self._add_content_metric()
        self._write_stat()
        self._dump_content()
        self._render_html()
        self._mtx2tsv()

    @runtime(f"{__name__}.debug_subprocess_call")
    def debug_subprocess_call(self, cmd):
        """debug subprocess call"""
        logger.debug(cmd)
        subprocess.check_call(cmd, shell=True)

    def get_metric_list(self):
        return self.__metric_list

    def set_metric_list(self, metric_list):
        self.__metric_list = metric_list

    def __exit__(self, *args, **kwargs):
        self.collect_matrix()
        self._clean_up()
