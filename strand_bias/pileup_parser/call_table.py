
from itertools import product

from numpy import log2
import pandas as pd
from pandas.api.extensions import register_dataframe_accessor
from pandas import IndexSlice as idx
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"


def _filter_dups(table, duplicate_keys):
    table = (table
             .reset_index(level="project_id")
             .drop(duplicate_keys)
             .reset_index().set_index(["project_id", "case_id", "reference"]))
    return table



def _simplify_for_reference_asymmetry(table, change_numerator, change_denominator,
                                      add_ratio_column=True):
    ratio_label = f"{change_numerator}{change_denominator}_ratio".lower()
    changes = pd.DataFrame([list(change) for change in [change_numerator, change_denominator]],
                           columns=["Ref", "Alt"])
    dups = (table.groupby(["case_id", "file_id"]).agg(sum)
            .reset_index()
            .case_id.value_counts())
    dups = set(dups.loc[dups > 1].to_dict().keys())
    table = table.sort_index()
    table = (table.groupby(["alt", "filtered"], axis=1).agg(sum)
             .droplevel("file_id"))
    table = _filter_dups(table, dups)
    table = table.loc[idx[:, :, list(changes.Ref)], idx[list(changes.Alt), :]]
    keepers = [(alt, thing, ref) for thing in ["PASS", "FAIL"]
               for _, (ref, alt) in changes.iterrows()]
    table = table.unstack(level=-1).loc[:, keepers]
    cols = list(zip(*table.columns))
    cols = [cols[1], ["".join(x) for x in zip(cols[2], cols[0])]]
    table.columns = pd.MultiIndex.from_arrays(cols, names=["filtered", "mutation"])
    totals = table.groupby("mutation", axis=1).agg(sum)
    for col in totals.columns:
        table[("TOTAL", col)] = totals[col]
    if add_ratio_column:
        for col in table.columns.levels[0]:
            table[(col, ratio_label)] = (table[(col, change_numerator)]
                                         / table[(col, change_denominator)])
    table.sort_index(axis=1, inplace=True)
    return table


@register_dataframe_accessor("call_tools")
class CallTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", header=[0, 1, 2], index_col=[0, 1, 2, 3])
        # table.axes[1].names = ["coding_region", "mutation", "snv_filter"]
        return table

    def simplify_for_reference_asymmetry(self, change_numerator, change_denominator,
                                         add_ratio_column=True, normalize_by_nt_content=False,
                                         normalization_counts=None, log_transform_ratio=False):
        changes = list(zip(change_numerator, change_denominator))

        table = (self._obj
                 .loc[idx[:, :, :, list(changes[0]), :, :], idx[list(changes[1]), :]]
                 .groupby("alt", axis=1).agg(sum)
                 .groupby(["project_id", "case_id", "reference"]).agg(sum)
                 )
        table = table.unstack(level=-1)[[tuple(reversed(change_numerator)),
                                         tuple(reversed(change_denominator))]]
        table.columns = [change_numerator, change_denominator]
        if normalize_by_nt_content:
            for col in table.columns:
                table[col] = table[col] / normalization_counts[col[0]]
        if add_ratio_column:
            ratio_label = f"{change_numerator}{change_denominator}_ratio".lower()
            if log_transform_ratio:
                table[ratio_label] = log2(table[change_numerator] / table[change_denominator])
            else:
                table[ratio_label] = table[change_numerator] / table[change_denominator]
        return table

    # def simplify_for_reference_asymmetry(self, change_numerator, change_denominator):
    #     table = _simplify_for_reference_asymmetry(self._obj, change_numerator, change_denominator)
    #     return table

    def simplify_reference_asymmetry(self):
        dups = (self._obj.groupby(["case_id", "file_id"]).agg(sum)
                .reset_index()
                .case_id.value_counts())
        dups = set(dups.loc[dups > 1].to_dict().keys())
        table = self._obj.sort_index()
        table = (table.groupby(["alt", "filtered"], axis=1).agg(sum)
                 .droplevel("file_id"))
        table = _filter_dups(table, dups)
        table = table.unstack(level=-1).reorder_levels([1, 2, 0], axis=1).sort_index(axis=1)

        totals = table.groupby(["reference", "alt"], axis=1).agg(sum)
        totals.columns = pd.MultiIndex.from_tuples([("TOTAL",) + col
                                                    for col in totals.columns])
        table = table.join(totals)
        return table

    def aggregate_total_by_coding_region(self):
        table = self._obj.groupby("coding_region", axis=1).agg(sum)
        return table

    def aggregate_total_by_mutation(self):
        table = self._obj.groupby("mutation", axis=1).agg(sum)
        return table

    def collapse_snv_filter(self):
        table = self._obj.groupby(["coding_region", "mutation"], axis=1).agg(sum)
        return table

    def collapse_regions(self):
        table = self._obj.groupby(["mutation", "snv_filter"], axis=1).agg(sum)
        return table

    def plot_boxplot_by_coding_region(self):
        table = self.aggregate_total_by_coding_region()
        plot = go.Figure()
        plot.add_trace(go.Box(y=table.forward, name="Forward Coding"))
        plot.add_trace(go.Box(y=table.reverse, name="Reverse Coding"))
        return plot

    def plot_scatter_by_coding_region(self):
        table = self.aggregate_total_by_coding_region()
        plot = go.Figure(go.Scatter(x=table.forward, y=table.reverse,
                                    mode="markers"))
        plot.update_layout(xaxis_title="Forward Coding",
                           yaxis_title="Reverse Coding")
        return plot

    def calculate_gtca_ratio(self):
        table = self._obj.asym_var_tools.aggregate_total_by_mutation()
        table["gtca_ratio"] = table.GT / table.CA
        return table

    def plot_gtca_scatter(self):
        plot = go.Figure(go.Scatter(x=self._obj.GT, y=self._obj.CA, mode="markers"))
        plot.add_trace(go.Scatter(x=[0, self._obj.GT.max()], y=[0, self._obj.GT.max()], mode="lines"))
        return plot

    def plot_gtca_ratio_boxplot(self):
        plot = go.Figure(go.Box(y=self._obj.gtca_ratio))
        return plot


@register_dataframe_accessor("change_tools")
class ChangeTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    def subset_change_pair(self, change_numerator, change_denominator):
        keepers = [(x, y, z) for x in ["FAIL", "PASS", "TOTAL"]
                   for y, z in [tuple(change_numerator), tuple(change_denominator)]]
        table = self._obj[keepers]
        table.columns = pd.MultiIndex.from_tuples([(x, y+z) for x, y, z in keepers],
                                                  names=["filtered", "mutation"])
        for filtered in table.columns.levels[0]:
            label = f"{change_numerator}{change_denominator}_ratio".lower()
            ratio = pd.DataFrame(table.loc[:, (filtered, change_numerator)]
                                 / table.loc[:, (filtered, change_denominator)])
            ratio.columns = pd.MultiIndex.from_tuples([(filtered, label)], names=["filtered", "mutation"])
            table = table.join(ratio)
        table = table.sort_index(axis=1)
        return table


@register_dataframe_accessor("call_tally_tools")
class CallTallyTable:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    @staticmethod
    def read_csv(file_path):
        table = pd.read_csv(file_path, sep="\t", index_col=list(range(8)))
        return table

    def simplify_for_reference_asymmetry(self, change_numerator, change_denominator,
                                         add_ratio_column=True, add_fraction_column=True,
                                         log_transform_ratio=False,
                                         normalization_counts=None, normalization_factor=1):
        table = self._obj.groupby(["project_id", "case_id", "file_id",
                                   "ref", "alt"]).agg(sum)
        table = table[["PASS", "FAIL"]]
        table = table.unstack(level=[-2, -1])
        keepers = [(x,) + y for x in ["PASS", "FAIL"]
                   for y in [tuple(change_numerator), tuple(change_denominator)]]
        table = table[keepers]
        table.columns = pd.MultiIndex.from_tuples([(x, y) for x in ["PASS", "FAIL"]
                                                   for y in [change_numerator, change_denominator]])

        if normalization_counts is not None:
            for col in table.columns:
                table[col] = table[col] / normalization_counts[col[1][0]] * normalization_factor

        if add_ratio_column:
            ratio_label = f"{change_numerator}{change_denominator}_ratio".lower()
            for filter_status in ["PASS", "FAIL"]:
                denom = (table[(filter_status, change_numerator)]
                         + table[(filter_status, change_denominator)])
                ratios = table[(filter_status, change_numerator)] / denom
                if log_transform_ratio:
                    table[(filter_status, ratio_label)] = log2(ratios)
                else:
                    table[(filter_status, ratio_label)] = ratios
            table = table[sorted(table.columns)]
        return table
